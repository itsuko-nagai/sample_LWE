#include <stdint.h>
#include <stddef.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>

#define KB_Q          3329
#define KB_N          256
#define KB_K          4
#define KB_ETA        2
#define KB_HALF_Q     1665
#define KB_SEED_BYTES 32
#define KB_MAX_BITS   (KB_K * KB_N)   /* 1024 bits available per generate() call */

typedef struct { int16_t c[KB_N]; } kb_poly;
typedef struct { kb_poly v[KB_K]; } kb_polyvec;

/* fixed key material -- set once by kb_setup(), never touched again */
static kb_poly    kb_A[KB_K][KB_K];
static kb_polyvec kb_s;
static kb_polyvec kb_As;          /* A.s, precomputed since both are fixed */
static int        kb_ready = 0;

static inline int16_t kb_mod_q(int32_t x)
{
    int32_t r = x % KB_Q;
    if (r < 0) r += KB_Q;
    return (int16_t)r;
}

static void kb_poly_add(const kb_poly *a, const kb_poly *b, kb_poly *out)
{
    for (int i = 0; i < KB_N; i++) out->c[i] = kb_mod_q((int32_t)a->c[i] + b->c[i]);
}

static void kb_poly_mul(const kb_poly *a, const kb_poly *b, kb_poly *out)
{
    int32_t acc[2 * KB_N - 1];
    memset(acc, 0, sizeof(acc));
    for (int i = 0; i < KB_N; i++) {
        if (a->c[i] == 0) continue;
        for (int j = 0; j < KB_N; j++)
            acc[i + j] += (int32_t)a->c[i] * (int32_t)b->c[j];
    }
    int32_t red[KB_N];
    memset(red, 0, sizeof(red));
    for (int i = 0; i < 2 * KB_N - 1; i++) {
        if (i < KB_N) red[i] += acc[i];
        else          red[i - KB_N] -= acc[i];
    }
    for (int i = 0; i < KB_N; i++) out->c[i] = kb_mod_q(red[i]);
}

static void kb_poly_zero(kb_poly *p) { memset(p->c, 0, sizeof(p->c)); }

static void kb_poly_uniform_from_seed(const uint8_t seed[KB_SEED_BYTES], uint8_t i, uint8_t j, kb_poly *out)
{
    uint8_t buf[KB_SEED_BYTES + 2];
    memcpy(buf, seed, KB_SEED_BYTES);
    buf[KB_SEED_BYTES]     = i;
    buf[KB_SEED_BYTES + 1] = j;

    uint8_t stream[1024];
    shake128(buf, KB_SEED_BYTES + 2, stream, sizeof(stream));

    int pos = 0, count = 0, extra = 0;
    while (count < KB_N) {
        if (pos + 3 > (int)sizeof(stream)) {
            uint8_t buf2[KB_SEED_BYTES + 3];
            memcpy(buf2, buf, KB_SEED_BYTES + 2);
            buf2[KB_SEED_BYTES + 2] = (uint8_t)(++extra);
            shake128(buf2, KB_SEED_BYTES + 3, stream, sizeof(stream));
            pos = 0;
        }
        uint16_t d1 = stream[pos] | ((uint16_t)(stream[pos + 1] & 0x0F) << 8);
        uint16_t d2 = (stream[pos + 1] >> 4) | ((uint16_t)stream[pos + 2] << 4);
        pos += 3;
        if (d1 < KB_Q && count < KB_N) out->c[count++] = (int16_t)d1;
        if (d2 < KB_Q && count < KB_N) out->c[count++] = (int16_t)d2;
    }
}

static void kb_gen_matrix_A(const uint8_t seed[KB_SEED_BYTES])
{
    for (uint8_t i = 0; i < KB_K; i++)
        for (uint8_t j = 0; j < KB_K; j++)
            kb_poly_uniform_from_seed(seed, i, j, &kb_A[i][j]);
}

static void kb_cbd_sample_poly(const uint8_t *stream, int eta, kb_poly *out)
{
    int bit_idx = 0;
    for (int n = 0; n < KB_N; n++) {
        int a = 0, b = 0;
        for (int t = 0; t < eta; t++) {
            int byte_i = bit_idx / 8, bit_i = bit_idx % 8;
            a += (stream[byte_i] >> bit_i) & 1;
            bit_idx++;
        }
        for (int t = 0; t < eta; t++) {
            int byte_i = bit_idx / 8, bit_i = bit_idx % 8;
            b += (stream[byte_i] >> bit_i) & 1;
            bit_idx++;
        }
        out->c[n] = kb_mod_q(a - b);
    }
}

static void kb_poly_noise(const uint8_t seed[KB_SEED_BYTES], uint8_t nonce, int eta, kb_poly *out)
{
    uint8_t buf[KB_SEED_BYTES + 1];
    memcpy(buf, seed, KB_SEED_BYTES);
    buf[KB_SEED_BYTES] = nonce;
    uint8_t stream[64 * 3];
    shake256(buf, KB_SEED_BYTES + 1, stream, (size_t)eta * 64);
    kb_cbd_sample_poly(stream, eta, out);
}

static void kb_polyvec_noise(const uint8_t seed[KB_SEED_BYTES], uint8_t nonce_base, int eta, kb_polyvec *out)
{
    for (int i = 0; i < KB_K; i++)
        kb_poly_noise(seed, (uint8_t)(nonce_base + i), eta, &out->v[i]);
}

static void kb_matrix_vec_mul(void)
{
    kb_poly tmp;
    for (int i = 0; i < KB_K; i++) {
        kb_poly_zero(&kb_As.v[i]);
        for (int j = 0; j < KB_K; j++) {
            kb_poly_mul(&kb_A[i][j], &kb_s.v[j], &tmp);
            kb_poly_add(&kb_As.v[i], &tmp, &kb_As.v[i]);
        }
    }
}

static void kb_get_random_bytes(uint8_t *out, size_t len)
{
    for (size_t i = 0; i < len; i++) out[i] = (uint8_t)(rand() & 0xFF);
}

/* ------------------------------------------------------------------
 * kb_setup()
 * Call ONCE, before any kb_generate_bits() calls. Generates A and s
 * from random seeds. Fixed for the rest of the program's life.
 * If your code already calls srand() elsewhere, delete the srand()
 * line here so you don't reseed rand() twice.
 * ------------------------------------------------------------------ */
static void kb_setup(void)
{
    srand((unsigned)time(NULL));

    uint8_t matrix_seed[KB_SEED_BYTES];
    kb_get_random_bytes(matrix_seed, KB_SEED_BYTES);
    kb_gen_matrix_A(matrix_seed);

    uint8_t secret_seed[KB_SEED_BYTES];
    kb_get_random_bytes(secret_seed, KB_SEED_BYTES);
    kb_polyvec_noise(secret_seed, 0, KB_ETA, &kb_s);

    kb_matrix_vec_mul();   /* precompute A.s */

    kb_ready = 1;
}

/* ------------------------------------------------------------------
 * kb_generate_bits(mu, bitbuf, nbits)
 * Computes c = A.s + e + mu*floor(q/2) (mod q) with a FRESH e every
 * call, and writes `nbits` output bits (0/1, one per byte) into
 * bitbuf, taken as the LSB of each coefficient of c, in order
 * c[0][0..255], c[1][0..255], c[2][0..255], c[3][0..255].
 * nbits must be <= KB_MAX_BITS (1024).
 * ------------------------------------------------------------------ */
static void kb_generate_bits(int mu, uint8_t *bitbuf, size_t nbits)
{
    if (!kb_ready) kb_setup();
    if (nbits > KB_MAX_BITS) nbits = KB_MAX_BITS;

    uint8_t noise_seed[KB_SEED_BYTES];
    kb_get_random_bytes(noise_seed, KB_SEED_BYTES);
    kb_polyvec e;
    kb_polyvec_noise(noise_seed, 0, KB_ETA, &e);

    kb_polyvec c;
    for (int i = 0; i < KB_K; i++) {
        kb_poly withNoise;
        kb_poly_add(&kb_As.v[i], &e.v[i], &withNoise);
        for (int n = 0; n < KB_N; n++)
            withNoise.c[n] = kb_mod_q(withNoise.c[n] + mu * KB_HALF_Q);
        c.v[i] = withNoise;
    }

    size_t written = 0;
    for (int i = 0; i < KB_K && written < nbits; i++)
        for (int n = 0; n < KB_N && written < nbits; n++)
            bitbuf[written++] = (uint8_t)(c.v[i].c[n] & 1);
}

/* ============================================================
 * mu ARRAY + "reseed every 2N bits" GLUE
 * (also just plain functions/globals, paste alongside the above)
 * ============================================================ */

#define KB_LFSR_N   128     /* <-- SET to your LFSR's actual N */
#define KB_SEED_BITS KB_LFSR_N
#define KB_MU_COUNT  4096   /* how many mu values to pre-generate */

static int  kb_mu_array[KB_MU_COUNT];
static int  kb_mu_index = 0;
static long kb_bits_since_reseed = 0;

static void kb_init_mu_array(void)
{
    /* rand() is already seeded by kb_setup(); no need to srand() again
     * here unless you call this before kb_setup(). */
    for (int i = 0; i < KB_MU_COUNT; i++)
        kb_mu_array[i] = rand() % 2;
}

/*
 * Call this once per output bit, right before you emit the bit,
 * e.g.:
 *     kb_maybe_reseed();
 *     int bit = your_lfsr_step();
 *     fputc(bit ? '1' : '0', fp);
 *     kb_bits_since_reseed++;
 *
 * It checks whether 2*KB_LFSR_N bits have been produced since the
 * last reseed, and if so, pulls the next mu, generates fresh seed
 * bits, and reseeds your LFSR -- you fill in the marked line below
 * with your team's real LFSR seed function.
 */
static void kb_maybe_reseed(void)
{
    if (kb_bits_since_reseed < 2 * KB_LFSR_N && kb_bits_since_reseed != 0)
        return; /* not time yet -- but always reseed on the very first call */

    if (kb_mu_index >= KB_MU_COUNT) kb_mu_index = 0; /* wrap, or re-init if you prefer */

    uint8_t seed_bits[KB_SEED_BITS];
    kb_generate_bits(kb_mu_array[kb_mu_index], seed_bits, KB_SEED_BITS);
    kb_mu_index++;

    /* ---- TODO: replace with your team's real LFSR seed call ---- */
    /* your_lfsr_load_state(seed_bits, KB_SEED_BITS); */

    kb_bits_since_reseed = 0;
}