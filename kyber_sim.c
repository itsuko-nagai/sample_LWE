/*
 * kyber_ring_bit_demo.c
 * ------------------------------------------------------------------
 * Exactly ONE equation, exactly as your mentor wrote it:
 *
 *      c = A . s + e + m * floor(q/2)     (mod q)
 *
 * but keeping the REAL ring structure underneath it:
 *   - A is a genuine 4x4 matrix of polynomials (16 polys total),
 *     each poly has 256 coefficients in Z_q[X]/(X^256+1).
 *   - A is generated from a seed via SHAKE128 + rejection sampling,
 *     exactly like gen_matrix_A() in the earlier ring demo.
 *   - s (the "key", K=4 polynomials) is generated from a seed via
 *     SHAKE256 + a centered-binomial-distribution (CBD) sampler,
 *     exactly like poly_noise()/cbd_sample_poly() before.
 *   - e (K=4 polynomials, the NOISE) is generated the same way --
 *     SHAKE256 + CBD -- but from a FRESH random seed every single
 *     call, which is the only thing that changes between runs.
 *
 * What is deliberately removed vs. the earlier ring demo:
 *   - No keygen()/encrypt_msg()/decrypt_msg() split.
 *   - No ciphertext pair (u, v). There is just ONE output vector c.
 *   - No 32-byte message. m is a single bit, 0 or 1, exactly what
 *     your mentor asked for ("he will give either 0 or 1").
 *   - m*floor(q/2) is added as a constant offset to every
 *     coefficient of every component of c -- the simplest possible
 *     way to fold a single bit into every polynomial in the vector
 *     while keeping the equation to that one line.
 *
 * A and s are computed ONCE at program start and never touch rand()
 * again -- they are "the key". e is resampled fresh on every run of
 * the loop below. Feed the same m in every time and watch c change
 * anyway, purely because of e.
 *
 * Build:  gcc -O2 -o kyber_ring_bit_demo sha3.c kyber_ring_bit_demo.c
 * Run:    ./kyber_ring_bit_demo 0      (bit 0, five times)
 *         ./kyber_ring_bit_demo 1      (bit 1, five times)
 * ------------------------------------------------------------------
 */

#include <stdio.h>
#include <stdint.h>
#include <stddef.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>

/* from sha3.c */
void shake128(const uint8_t *in, size_t inlen, uint8_t *out, size_t outlen);
void shake256(const uint8_t *in, size_t inlen, uint8_t *out, size_t outlen);

/* ---------------- Kyber1024-style parameters ---------------- */
#define Q          3329
#define N          256      /* ring degree: Z_q[X]/(X^N+1)        */
#define K          4        /* matrix rank -> 4x4 matrix of polys */
#define ETA        2        /* CBD noise width (Kyber1024's eta1) */
#define HALF_Q     1665     /* ceil(Q/2)                          */
#define SEED_BYTES 32        /* real Kyber spec: seeds are 32 bytes */

typedef struct { int16_t c[N]; } poly;
typedef struct { poly v[K]; } polyvec;

/* ---------------- modular helper ---------------- */
static inline int16_t mod_q(int32_t x)
{
    int32_t r = x % Q;
    if (r < 0) r += Q;
    return (int16_t)r;
}

/* ---------------- polynomial add ---------------- */
static void poly_add(const poly *a, const poly *b, poly *out)
{
    for (int i = 0; i < N; i++) out->c[i] = mod_q((int32_t)a->c[i] + b->c[i]);
}

/* Schoolbook negacyclic convolution in Z_q[X]/(X^N+1). */
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
        if (i < N) red[i] += acc[i];
        else       red[i - N] -= acc[i]; /* X^N = -1 wraparound */
    }
    for (int i = 0; i < N; i++) out->c[i] = mod_q(red[i]);
}

static void poly_zero(poly *p) { memset(p->c, 0, sizeof(p->c)); }

/* ---------------- SHAKE128 + rejection sampling -> one A[i][j] ---------------- */
static void poly_uniform_from_seed(const uint8_t seed[SEED_BYTES], uint8_t i, uint8_t j, poly *out)
{
    uint8_t buf[SEED_BYTES + 2];
    memcpy(buf, seed, SEED_BYTES);
    buf[SEED_BYTES]     = i;
    buf[SEED_BYTES + 1] = j;

    uint8_t stream[1024];
    shake128(buf, SEED_BYTES + 2, stream, sizeof(stream));

    int pos = 0, count = 0, extra = 0;
    while (count < N) {
        if (pos + 3 > (int)sizeof(stream)) {
            uint8_t buf2[SEED_BYTES + 3];
            memcpy(buf2, buf, SEED_BYTES + 2);
            buf2[SEED_BYTES + 2] = (uint8_t)(++extra);
            shake128(buf2, SEED_BYTES + 3, stream, sizeof(stream));
            pos = 0;
        }
        /* rejection sampling: two 12-bit candidates per 3 bytes,
         * keep only those < Q, exactly like Kyber's poly_uniform() */
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

/* ---------------- SHAKE256 + CBD -> one noise polynomial ---------------- */
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

static void poly_noise(const uint8_t seed[SEED_BYTES], uint8_t nonce, int eta, poly *out)
{
    uint8_t buf[SEED_BYTES + 1];
    memcpy(buf, seed, SEED_BYTES);
    buf[SEED_BYTES] = nonce;
    uint8_t stream[64 * 3]; /* eta*64 bytes needed; eta<=3 so 192 is plenty */
    shake256(buf, SEED_BYTES + 1, stream, (size_t)eta * 64);
    cbd_sample_poly(stream, eta, out);
}

static void polyvec_noise(const uint8_t seed[SEED_BYTES], uint8_t nonce_base, int eta, polyvec *out)
{
    for (int i = 0; i < K; i++)
        poly_noise(seed, (uint8_t)(nonce_base + i), eta, &out->v[i]);
}

/* out[i] = sum_j A[i][j]*s[j]  -- the "A . s" matrix-vector product */
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

/* ---------------- random seed helper ---------------- */
static void get_random_bytes(uint8_t *out, size_t len)
{
    for (size_t i = 0; i < len; i++) out[i] = (uint8_t)(rand() & 0xFF);
}

static void print_poly_head(const char *label, const poly *p, int n)
{
    printf("%s[", label);
    for (int i = 0; i < n; i++) printf("%4d%s", p->c[i], i == n - 1 ? "" : ", ");
    printf(", ...]\n");
}

/* Print every coefficient of a polynomial in BINARY, 12 bits each
 * (Q=3329 needs 12 bits: 2^11=2048 is too few, 2^12=4096 is enough
 * to cover 0..3328). */
static void print_poly_binary(const char *label, const poly *p)
{
    printf("%s\n", label);
    for (int i = 0; i < N; i++) {
        printf("  c[%3d] = %4d = ", i, p->c[i]);
        for (int b = 11; b >= 0; b--) putchar((p->c[i] >> b) & 1 ? '1' : '0');
        putchar('\n');
    }
}

int main(int argc, char **argv)
{
    srand((unsigned)time(NULL));

    int m = (argc > 1) ? atoi(argv[1]) : 1;
    if (m != 0 && m != 1) {
        fprintf(stderr, "usage: %s [0|1]\n", argv[0]);
        return 1;
    }

    printf("c = A . s + e + m*floor(q/2)   (mod q)\n");
    printf("Q=%d  N=%d (coeffs/poly)  K=%d (matrix rank -> %dx%d = %d polys)  ETA=%d\n",
           Q, N, K, K, K, K * K, ETA);
    printf("m = %d, fixed for every run below\n\n", m);

    /* ---- A: 4x4 matrix, 16 polynomials, from SHAKE128 + rejection
     * sampling. Generated ONCE, from a random seed. Fixed key material. ---- */
    uint8_t matrix_seed[SEED_BYTES];
    get_random_bytes(matrix_seed, SEED_BYTES);
    poly A[K][K];
    gen_matrix_A(matrix_seed, A);

    /* ---- s: 4 polynomials, from SHAKE256 + CBD. Generated ONCE,
     * from a random seed. Fixed key material, same as A. ---- */
    uint8_t secret_seed[SEED_BYTES];
    get_random_bytes(secret_seed, SEED_BYTES);
    polyvec s;
    polyvec_noise(secret_seed, 0, ETA, &s);

    /* A . s is fixed too, since neither A nor s ever change again */
    polyvec As;
    matrix_vec_mul((const poly (*)[K])A, &s, &As);

    printf("A[0][0] (first 6 of 256 coeffs, fixed for the whole run):\n");
    print_poly_head("  ", &A[0][0], 6);
    printf("s[0]    (first 6 of 256 coeffs, fixed for the whole run):\n");
    print_poly_head("  ", &s.v[0], 6);
    printf("A.s[0]  (first 6 of 256 coeffs, fixed for the whole run):\n");
    print_poly_head("  ", &As.v[0], 6);

    //printf("\n only thing that changes: e, resampled every run.\n");
    printf("------------------------------------------------------------\n");

    polyvec c; /* keeps the LAST run's result for the full binary dump below */

    for (int run = 1; run <= 5; run++) {
        /* fresh random seed EVERY run -> fresh e EVERY run */
        uint8_t noise_seed[SEED_BYTES];
        get_random_bytes(noise_seed, SEED_BYTES);
        polyvec e;
        polyvec_noise(noise_seed, 0, ETA, &e);

        /* c = A.s + e + m*floor(q/2), broadcast onto every component */
        for (int i = 0; i < K; i++) {
            poly withNoise;
            poly_add(&As.v[i], &e.v[i], &withNoise);
            for (int n = 0; n < N; n++)
                withNoise.c[n] = mod_q(withNoise.c[n] + m * HALF_Q);
            c.v[i] = withNoise;
        }

        printf("\nrun %d:\n", run);
        print_poly_head("  e[0] = ", &e.v[0], 6);
        print_poly_head("  c[0] = ", &c.v[0], 6);
    }

    // printf("\nSame m=%d, same A, same s every run -- only e changed --\n"
    //        "and c[0] came out different every single time.\n", m);

    /* ---- full binary dump of the FINAL matrix output c ----
     * c is a polyvec: K=4 polynomials, N=256 coefficients each,
     * every coefficient shown as 12 bits (Q=3329 needs 12 bits). */
    // printf("\n============================================================\n");
    // printf("FINAL MATRIX OUTPUT c (last run, m=%d), full binary dump:\n", m);
    // printf("c has K=%d polynomials, N=%d coefficients each, 12 bits/coeff\n", K, N);
    // printf("============================================================\n");
    // for (int i = 0; i < K; i++) {
    //     char label[32];
    //     snprintf(label, sizeof(label), "\nc[%d]:", i);
    //     print_poly_binary(label, &c.v[i]);
    // }

    return 0;
}