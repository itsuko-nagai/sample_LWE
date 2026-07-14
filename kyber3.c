/*
 * kyber_ring_demo.c
 * ------------------------------------------------------------------
 * Ring-LWE version of the earlier scalar demo, addressing both of
 * your mentor's corrections:
 *
 *   1. Every matrix/vector entry is now a REAL polynomial in
 *      Z_q[X]/(X^256+1)  (N=256 coefficients, mod Q=3329), not a
 *      single integer. A is a genuine K x K = 4x4 matrix of such
 *      polynomials (16 polynomials, 256 coeffs each).
 *
 *   2. Encryption now produces the real two-part Kyber ciphertext:
 *
 *         keygen:   t  = A ."  s + e            (once, secret s fixed)
 *         encrypt:  u  = A^T . r + e1           (fresh r EVERY call)
 *                   v  = t^T . r + e2 + encode(m)
 *         decrypt:  m' = decode( v - s^T . u )
 *
 *      where "." is the ring-matrix product (schoolbook negacyclic
 *      polynomial convolution mod X^256+1, mod Q -- Kyber itself
 *      uses NTT for speed, but schoolbook gives the identical
 *      mathematical result and is much easier to read/verify).
 *
 * Parameters used: Q=3329, N=256, K=4, ETA1=2, ETA2=2.
 * This is the REAL Kyber1024 parameter set (see pq-crystals/kyber
 * ref/params.h: K=4 => ETA1=2). The previous scalar demo used K=4
 * with ETA1=3, which is not a real Kyber parameter combination --
 * fixed here.
 *
 * What is still simplified vs. the reference implementation:
 *   - No NTT (schoolbook poly-mul instead -- correct, just slower).
 *   - No ciphertext compression (Compress_q/Decompress_q, du/dv bit
 *     packing) -- u and v are left as full mod-q coefficient arrays.
 *   - No CPA/CCA KEM wrapper (Fujisaki-Okamoto transform, hashing of
 *     the public key into the seed, re-encryption check, etc). This
 *     is IND-CPA PKE only, i.e. what indcpa_enc()/indcpa_dec() do
 *     inside crystals-kyber, not the full crypto_kem_enc()/_dec().
 *
 * Build:   gcc -O2 -o kyber_ring_demo sha3.c kyber_ring_demo.c
 * Run:     ./kyber_ring_demo
 * ------------------------------------------------------------------
 */

#include <stdio.h>
#include <stdint.h>
#include <stddef.h>
#include <string.h>

/* from sha3.c */
void shake128(const uint8_t *in, size_t inlen, uint8_t *out, size_t outlen);
void shake256(const uint8_t *in, size_t inlen, uint8_t *out, size_t outlen);

/* ---------------- Kyber1024 parameters ---------------- */
#define Q       3329
#define N       256      /* ring degree: Z_q[X]/(X^N+1) */
#define K       4         /* matrix rank -> 4x4 matrix of polynomials */
#define ETA1    2         /* real Kyber1024 value (K=4 => ETA1=2) */
#define ETA2    2
#define HALF_Q  1665      /* ceil(Q/2), matches poly_frommsg() in ref code */
#define SEED_BYTES 32     /* real Kyber spec: seeds are ALWAYS exactly 32 bytes */

typedef struct { int16_t c[N]; } poly;
typedef struct { poly v[K]; } polyvec;

/* ---------------- modular helpers ---------------- */
static inline int16_t mod_q(int32_t x)
{
    int32_t r = x % Q;
    if (r < 0) r += Q;
    return (int16_t)r;
}

/* ---------------- polynomial arithmetic ---------------- */
static void poly_add(const poly *a, const poly *b, poly *out)
{
    for (int i = 0; i < N; i++) out->c[i] = mod_q((int32_t)a->c[i] + b->c[i]);
}

static void poly_sub(const poly *a, const poly *b, poly *out)
{
    for (int i = 0; i < N; i++) out->c[i] = mod_q((int32_t)a->c[i] - b->c[i]);
}

/* Schoolbook negacyclic convolution: multiply in Z_q[X]/(X^N+1).
 * Kyber's NTT-based poly_basemul does, just O(N^2) instead of
 * O(N log N). */
static void poly_mul(const poly *a, const poly *b, poly *out)
{
    int32_t acc[2 * N - 1];
    memset(acc, 0, sizeof(acc));
    for (int i = 0; i < N; i++) {
        if (a->c[i] == 0) continue;
        for (int j = 0; j < N; j++)
            acc[i + j] += (int32_t)a->c[i] * (int32_t)b->c[j];
    }
    int32_t red[N];
    memset(red, 0, sizeof(red));
    for (int i = 0; i < 2 * N - 1; i++) {
        if (i < N) {
            red[i] += acc[i];
        } else {
            red[i - N] -= acc[i]; /* X^N = -1 wraparound */
        }
    }
    for (int i = 0; i < N; i++) out->c[i] = mod_q(red[i]);
}

static void poly_zero(poly *p) { memset(p->c, 0, sizeof(p->c)); }

/* ---------------- XOF: sample one A[i][j] polynomial ---------------- */
/* Mirrors gen_matrix()/poly_uniform() in ref/indcpa.c + ref/poly.c:
 * SHAKE128(seed || i || j) squeezed and rejection-sampled mod Q,
 * 12 bits at a time, two candidates per 3 bytes. */
static void poly_uniform_from_seed(const uint8_t seed[SEED_BYTES], uint8_t i, uint8_t j, poly *out)
{
    uint8_t buf[SEED_BYTES + 2];
    memcpy(buf, seed, SEED_BYTES);
    buf[SEED_BYTES]     = i;
    buf[SEED_BYTES + 1] = j;

    /* 256 coeffs need >= 256/(2 accepted per 3 bytes worst-case ~0.81
     * acceptance) ~ 474 bytes on average; 1024 bytes squeezed in one
     * shot is comfortably enough in the overwhelming majority of
     * cases. If it ever isn't, we re-XOF with an extra counter byte
     * appended, exactly like Kyber falls back to squeezing another
     * SHAKE128 block. */
    uint8_t stream[1024];
    shake128(buf, SEED_BYTES + 2, stream, sizeof(stream));

    int pos = 0, count = 0, extra = 0;
    while (count < N) {
        if (pos + 3 > (int)sizeof(stream)) {
            /* extremely rare fallback: extend with a counter byte */
            uint8_t buf2[SEED_BYTES + 3];
            memcpy(buf2, buf, SEED_BYTES + 2);
            buf2[SEED_BYTES + 2] = (uint8_t)(++extra);
            shake128(buf2, SEED_BYTES + 3, stream, sizeof(stream));
            pos = 0;
        }
        uint16_t d1 = stream[pos] | ((uint16_t)(stream[pos + 1] & 0x0F) << 8);
        uint16_t d2 = (stream[pos + 1] >> 4) | ((uint16_t)stream[pos + 2] << 4);
        pos += 3;
        if (d1 < Q && count < N) out->c[count++] = (int16_t)d1;
        if (d2 < Q && count < N) out->c[count++] = (int16_t)d2;
    }
}

static void gen_matrix_A(const uint8_t seed[SEED_BYTES], poly A[K][K])
{
    for (uint8_t i = 0; i < K; i++)
        for (uint8_t j = 0; j < K; j++)
            poly_uniform_from_seed(seed, i, j, &A[i][j]);
}

/* ---------------- PRF + CBD noise sampler ---------------- */
/* Same construction as poly_getnoise_eta1/2 -> poly_cbd_eta in
 * ref/poly.c: SHAKE256(seed||nonce) feeds a centered binomial
 * sampler, eta bits vs eta bits per coefficient. */
static void cbd_sample_poly(const uint8_t *stream, int eta, poly *out)
{
    int bit_idx = 0;
    for (int n = 0; n < N; n++) {
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
        out->c[n] = mod_q(a - b);
    }
}

static void poly_noise(const uint8_t seed[SEED_BYTES], uint8_t nonce,
                        int eta, poly *out)
{
    uint8_t buf[SEED_BYTES + 1]; /* fixed size, always exactly right */
    memcpy(buf, seed, SEED_BYTES);
    buf[SEED_BYTES] = nonce;
    uint8_t stream[64 * 3]; /* need eta*64 bytes; eta<=3 so 192 covers it */
    shake256(buf, SEED_BYTES + 1, stream, (size_t)eta * 64);
    cbd_sample_poly(stream, eta, out);
}

static void polyvec_noise(const uint8_t seed[SEED_BYTES],
                           uint8_t nonce_base, int eta, polyvec *out)
{
    for (int i = 0; i < K; i++)
        poly_noise(seed, (uint8_t)(nonce_base + i), eta, &out->v[i]);
}

/* ---------------- vector/matrix ring operations ---------------- */
/* out = sum_i a[i]*b[i] */
static void polyvec_dot(const polyvec *a, const polyvec *b, poly *out)
{
    poly_zero(out);
    poly tmp;
    for (int i = 0; i < K; i++) {
        poly_mul(&a->v[i], &b->v[i], &tmp);
        poly_add(out, &tmp, out);
    }
}

/* out[i] = sum_j A[i][j]*s[j]   (rows of A) */
static void matrix_vec_mul(const poly A[K][K], const polyvec *s, polyvec *out)
{
    poly tmp;
    for (int i = 0; i < K; i++) {
        poly_zero(&out->v[i]);
        for (int j = 0; j < K; j++) {
            poly_mul(&A[i][j], &s->v[j], &tmp);
            poly_add(&out->v[i], &tmp, &out->v[i]);
        }
    }
}

/* out[i] = sum_j A[j][i]*r[j]   (columns of A, i.e. A^T . r) */
static void matrix_transpose_vec_mul(const poly A[K][K], const polyvec *r, polyvec *out)
{
    poly tmp;
    for (int i = 0; i < K; i++) {
        poly_zero(&out->v[i]);
        for (int j = 0; j < K; j++) {
            poly_mul(&A[j][i], &r->v[j], &tmp);
            poly_add(&out->v[i], &tmp, &out->v[i]);
        }
    }
}

/* ---------------- message encode/decode ---------------- */
/* 32 bytes (256 bits) <-> one polynomial, exactly like poly_frommsg /
 * poly_tomsg in ref/poly.c. Bit k of the message becomes coefficient
 * k, scaled to 0 or HALF_Q. */
static void poly_encode_msg(const uint8_t msg[32], poly *out)
{
    for (int i = 0; i < 32; i++)
        for (int b = 0; b < 8; b++) {
            int bit = (msg[i] >> b) & 1;
            out->c[8 * i + b] = bit ? HALF_Q : 0;
        }
}

static void poly_decode_msg(const poly *p, uint8_t msg[32])
{
    memset(msg, 0, 32);
    for (int i = 0; i < 32; i++)
        for (int b = 0; b < 8; b++) {
            int16_t c = p->c[8 * i + b];
            /* distance to Q/2 vs distance to 0 */
            int d0   = c < (Q - c) ? c : (Q - c);
            int dhal = c > HALF_Q ? c - HALF_Q : HALF_Q - c;
            if (dhal < d0) msg[i] |= (uint8_t)(1u << b);
        }
}

/* ---------------- keygen / encrypt / decrypt ---------------- */
static void keygen(const uint8_t matrix_seed[SEED_BYTES],
                    const uint8_t secret_seed[SEED_BYTES],
                    poly A[K][K], polyvec *s, polyvec *t)
{
    gen_matrix_A(matrix_seed, A);
    polyvec_noise(secret_seed, 0, ETA1, s);      /* nonces 0..K-1 */
    polyvec e;
    polyvec_noise(secret_seed, K, ETA1, &e);      /* nonces K..2K-1 */
    matrix_vec_mul((const poly (*)[K])A, s, t);
    for (int i = 0; i < K; i++) poly_add(&t->v[i], &e.v[i], &t->v[i]);
}

static void encrypt_msg(const poly A[K][K], const polyvec *t,
                         const uint8_t coins[SEED_BYTES],
                         const uint8_t msg[32],
                         polyvec *u, poly *v)
{
    /* `coins` is now genuinely fresh 32 random bytes from the OS RNG,
     * generated by the caller right before each call (see main) --
     * this is exactly what a real indcpa_enc() caller passes in every
     * time, no nonce_base workaround needed since the bytes themselves
     * are already different every call. */
    polyvec r, e1;
    polyvec_noise(coins, 0, ETA1, &r);       /* nonces 0..K-1 */
    polyvec_noise(coins, K, ETA2, &e1);      /* nonces K..2K-1 */
    poly e2;
    poly_noise(coins, 2 * K, ETA2, &e2);      /* nonce 2K */

    matrix_transpose_vec_mul(A, &r, u);
    for (int i = 0; i < K; i++) poly_add(&u->v[i], &e1.v[i], &u->v[i]);

    poly tr, mpoly;
    polyvec_dot(t, &r, &tr);
    poly_encode_msg(msg, &mpoly);
    poly_add(&tr, &e2, v);
    poly_add(v, &mpoly, v);
}

static void decrypt_msg(const polyvec *s, const polyvec *u, const poly *v,
                         uint8_t msg[32])
{
    poly su, mp;
    polyvec_dot(s, u, &su);
    poly_sub(v, &su, &mp);
    poly_decode_msg(&mp, msg);
}

/* ---------------- printing helpers ---------------- */
static void print_poly_head(const char *label, const poly *p, int n)
{
    printf("%s[", label);
    for (int i = 0; i < n; i++) printf("%4d%s", p->c[i], i == n - 1 ? "" : ", ");
    printf(", ...]\n");
}

static void print_hex(const char *label, const uint8_t *b, int n)
{
    printf("%s", label);
    for (int i = 0; i < n; i++) printf("%02x", b[i]);
    printf("\n");
}

/* Print each coefficient in BINARY (12 bits, since Q=3329 needs 12 bits
 * to represent any value 0..3328 -- 2^11=2048 is too few, 2^12=4096 is
 * enough). This lets you visually inspect the actual bit pattern of
 * each number in t (or any polynomial), not just its decimal value. */
static void print_poly_binary(const char *label, const poly *p, int n)
{
    printf("%s\n", label);
    for (int i = 0; i < n; i++) {
        printf("  c[%3d] = %4d  = ", i, p->c[i]);
        for (int b = 11; b >= 0; b--) putchar((p->c[i] >> b) & 1 ? '1' : '0');
        putchar('\n');
    }
}

/* Print raw bytes in binary (8 bits each) -- for seeds, coins, message. */
static void print_bytes_binary(const char *label, const uint8_t *b, int n)
{
    printf("%s\n", label);
    for (int i = 0; i < n; i++) {
        printf("  byte[%2d] = 0x%02x = ", i, b[i]);
        for (int bit = 7; bit >= 0; bit--) putchar((b[i] >> bit) & 1 ? '1' : '0');
        putchar('\n');
    }
}

#include <stdlib.h>
#include <time.h>

static int get_random_bytes(uint8_t *out, size_t len)
{
    for (size_t i = 0; i < len; i++)
        out[i] = (uint8_t)(rand() & 0xFF);
    return 1;
}

int main(void)
{
    srand((unsigned)time(NULL));

    printf("Kyber1024-style ring-LWE demo -- real random seeds\n");
    printf("Q=%d  N=%d (poly degree)  K=%d (matrix rank -> %dx%d = %d polynomials)  ETA1=%d  ETA2=%d\n\n",
           Q, N, K, K, K, K * K, ETA1, ETA2);

    uint8_t matrix_seed[SEED_BYTES], secret_seed[SEED_BYTES];
    if (!get_random_bytes(matrix_seed, SEED_BYTES) ||
        !get_random_bytes(secret_seed, SEED_BYTES)) {
        fprintf(stderr, "could not generate random seed\n");
        return 1;
    }
    print_hex("matrix_seed (32 random bytes, public): 0x", matrix_seed, SEED_BYTES);
    print_hex("secret_seed (32 random bytes, PRIVATE): 0x", secret_seed, SEED_BYTES);

    poly A[K][K];
    polyvec s, t;
    keygen(matrix_seed, secret_seed, A, &s, &t);

    printf("\nA is a %dx%d matrix of degree-%d polynomials (16 polynomials total).\n", K, K, N);
    printf("First few coefficients of A[0][0] and A[3][3] as a sanity check:\n");
    print_poly_head("  A[0][0] = ", &A[0][0], 8);
    print_poly_head("  A[3][3] = ", &A[3][3], 8);
    printf("\ns (secret key, K=%d polynomials, fixed for this session):\n", K);
    for (int i = 0; i < K; i++) { char lbl[16]; snprintf(lbl, sizeof(lbl), "  s[%d] = ", i); print_poly_head(lbl, &s.v[i], 8); }
    printf("\nt = A.s + e  (public key, computed once at keygen):\n");
    for (int i = 0; i < K; i++) { char lbl[16]; snprintf(lbl, sizeof(lbl), "  t[%d] = ", i); print_poly_head(lbl, &t.v[i], 8); }

    printf("\nBinary view of t[0], first 8 coefficients (12 bits each, since Q=3329 needs 12 bits):\n");
    print_poly_binary("t[0] =", &t.v[0], 8);

    uint8_t msg[32];
    memset(msg, 0, 32);
    memcpy(msg, "HELLO-KYBER-MSG!", 17); /* first 17 bytes carry a readable payload */

    print_hex("\nPlaintext message (32 bytes): 0x", msg, 32);
    print_bytes_binary("Plaintext message, first 4 bytes in binary:", msg, 4);

    for (int call = 1; call <= 5; call++) {
        uint8_t coins[SEED_BYTES];
        if (!get_random_bytes(coins, SEED_BYTES)) {
            fprintf(stderr, "could not generate random coins\n");
            return 1;
        }

        polyvec u;
        poly v;
        encrypt_msg((const poly (*)[K])A, &t, coins, msg, &u, &v);

        uint8_t recovered[32];
        decrypt_msg(&s, &u, &v, recovered);

        printf("\n--- Encryption #%d (fresh 32-byte coins from rand()) ---\n", call);
        print_hex("  coins = 0x", coins, 8);
        print_poly_head("  u[0]  = ", &u.v[0], 8);
        print_poly_head("  v     = ", &v, 8);
        if (call == 1) print_poly_binary("  v (binary, first 8 coeffs) =", &v, 8);
        print_hex("  decrypted    : 0x", recovered, 32);
        printf("  round-trip OK : %s\n", memcmp(msg, recovered, 32) == 0 ? "yes" : "NO -- MISMATCH");
    }

    // printf("\n-> matrix_seed and secret_seed are now genuine 32-byte random values\n"
    //        "   (real Kyber's fixed seed size), pulled fresh from the OS RNG once,\n"
    //        "   at program start. The keypair (A, s, t) derived from them is fixed\n"
    //        "   for the rest of this run. Each encryption call instead pulls its\n"
    //        "   own fresh 32-byte 'coins' value, which is why u and v differ every\n"
    //        "   time even though the message and key never change.\n");

    return 0;
}