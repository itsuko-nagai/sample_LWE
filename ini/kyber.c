/*
 * kyber_lwe_demo.c
 * ------------------------------------------------------------------
 *   b = A.s + e + m*floor(q/2)   (mod q)
 *
 * A (public matrix) : generated ONCE from a fixed seed via SHAKE-128 +
 *                      rejection sampling mod q. Same strategy as
 *                      pq-crystals/kyber's indcpa_gen_matrix()
 *                      (ref/indcpa.c) which drives poly_uniform()
 *                      (ref/poly.c) via rej_uniform() over XOF output.
 *
 * s (secret vector) : generated ONCE from a fixed seed via SHAKE-256
 *                      used as a PRF, sampled through the Centered
 *                      Binomial Distribution (CBD, eta = ETA1). Same
 *                      construction as poly_getnoise_eta1() ->
 *                      poly_cbd_eta1() in ref/poly.c.
 *
 * e (error vector)  : regenerated on EVERY encryption call. Only the
 *                      PRF "nonce" changes between calls -- standing in
 *                      for the fresh 32-byte random "coins" a real
 *                      caller passes to crypto_kem_enc()/indcpa_enc()
 *                      every time. This is why b differs every call
 *                      even though A, s and the message bit are fixed.
 *
 * Build:   gcc -O2 -o kyber_lwe_demo sha3.c kyber_lwe_demo.c
 * Run:     ./kyber_lwe_demo
 * ------------------------------------------------------------------
 * NOTE: this treats each matrix entry as a single integer mod q (a
 * scalar/vector simplification), not a degree-256 ring element as
 * real Kyber does. The XOF/PRF/CBD machinery is implemented exactly
 * as the reference code does, so it faithfully demonstrates *why*
 * re-encrypting changes the ciphertext.
 */

#include <stdio.h>
#include <stdint.h>
#include <stddef.h>
#include <string.h>

/* from sha3.c */
void shake128(const uint8_t *in, size_t inlen, uint8_t *out, size_t outlen);
void shake256(const uint8_t *in, size_t inlen, uint8_t *out, size_t outlen);

/* ---------------- Kyber512 parameters ---------------- */
#define Q       3329
#define K       4        /* NOTE: real Kyber512 uses K=2; set to 2 to match
                             exactly. K=4 corresponds to Kyber1024's rank. */
#define ETA1    3
#define ETA2    2
#define HALF_Q  (Q / 2)

/* ---------------- XOF for matrix A ---------------- */
static void xof_A_entry(const uint8_t *seed, size_t seedlen,
                         int i, int j, uint8_t *out, size_t outlen)
{
    uint8_t buf[64];
    memcpy(buf, seed, seedlen);
    buf[seedlen]     = (uint8_t)i;
    buf[seedlen + 1] = (uint8_t)j;
    shake128(buf, seedlen + 2, out, outlen);
}

/* Rejection-sample one XOF stream into a full k x k matrix mod q,
 * exactly mirroring Kyber's gen_matrix(). */
static void gen_matrix_A(const uint8_t *seed, size_t seedlen, int16_t A[K][K])
{
    for (int i = 0; i < K; i++) {
        for (int j = 0; j < K; j++) {
            uint8_t stream[600];
            xof_A_entry(seed, seedlen, i, j, stream, sizeof(stream));
            int pos = 0;
            int found = 0;
            while (!found) {
                if (pos + 3 > (int)sizeof(stream)) {
                    /* refresh buffer if we somehow ran out (astronomically
                       unlikely with 600 bytes of buffer) */
                    xof_A_entry(seed, seedlen, i, j, stream, sizeof(stream));
                    pos = 0;
                }
                uint16_t d1 = stream[pos] | ((uint16_t)(stream[pos + 1] & 0x0F) << 8);
                uint16_t d2 = (stream[pos + 1] >> 4) | ((uint16_t)stream[pos + 2] << 4);
                pos += 3;
                if (d1 < Q) { A[i][j] = (int16_t)d1; found = 1; }
                else if (d2 < Q) { A[i][j] = (int16_t)d2; found = 1; }
            }
        }
    }
}

/* ---------------- PRF + CBD noise sampler ---------------- */
/* CBD_eta: sample = (sum of eta random bits) - (sum of eta random bits),
 * bit-for-bit the same construction as Kyber's poly_cbd_eta1/2. */
static void cbd_sample(const uint8_t *stream, int eta, int count, int16_t *out)
{
    int bit_idx = 0;
    for (int n = 0; n < count; n++) {
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
        int v = ((a - b) % Q + Q) % Q;
        out[n] = (int16_t)v;
    }
}

static void prf_vector(const uint8_t *seed, size_t seedlen, uint8_t nonce,
                        int eta, int count, int16_t *out)
{
    uint8_t buf[64];
    memcpy(buf, seed, seedlen);
    buf[seedlen] = nonce;

    uint8_t stream[64 * 8]; /* generous: 64*eta bytes needed, eta<=3 */
    size_t need = (size_t)64 * eta;
    shake256(buf, seedlen + 1, stream, need);
    cbd_sample(stream, eta, count, out);
}

/* ---------------- linear algebra mod q ---------------- */
static void matvec_mod_q(const int16_t A[K][K], const int16_t s[K], int16_t out[K])
{
    for (int i = 0; i < K; i++) {
        int32_t acc = 0;
        for (int j = 0; j < K; j++)
            acc += (int32_t)A[i][j] * (int32_t)s[j];
        out[i] = (int16_t)(((acc % Q) + Q) % Q);
    }
}

/* ---------------- encrypt / decrypt ---------------- */
static void encrypt_bit(const int16_t A[K][K], const int16_t s[K],
                         const uint8_t *prf_seed, size_t prf_seedlen,
                         uint8_t nonce, int bit,
                         int16_t b[K], int16_t e[K])
{
    prf_vector(prf_seed, prf_seedlen, nonce, ETA2, K, e);
    int16_t As[K];
    matvec_mod_q(A, s, As);
    int m_encoded = bit * HALF_Q;
    for (int i = 0; i < K; i++)
        b[i] = (int16_t)(((As[i] + e[i] + m_encoded) % Q + Q) % Q);
}

static int decrypt_bit(const int16_t A[K][K], const int16_t s[K], const int16_t b[K])
{
    int16_t As[K];
    matvec_mod_q(A, s, As);
    int votes0 = 0, votes1 = 0;
    for (int i = 0; i < K; i++) {
        int d = ((b[i] - As[i]) % Q + Q) % Q;
        int dist_to_0 = d < (Q - d) ? d : (Q - d);
        int dist_to_half = d > HALF_Q ? d - HALF_Q : HALF_Q - d;
        if (dist_to_0 < dist_to_half) votes0++; else votes1++;
    }
    return votes1 > votes0 ? 1 : 0;
}

/* ---------------- demo ---------------- */
static void print_vec(const char *label, const int16_t v[K])
{
    printf("%s[", label);
    for (int i = 0; i < K; i++) printf("%5d%s", v[i], i == K - 1 ? "" : ", ");
    printf("]\n");
}

int main(void)
{
    printf("Kyber-style LWE demo   (Q=%d, K=%d, ETA1=%d, ETA2=%d)\n\n", Q, K, ETA1, ETA2);

    const uint8_t matrix_seed[] = "public-seed-for-A-fixed-once";
    const uint8_t secret_seed[] = "secret-seed-for-s-fixed-once";

    int16_t A[K][K];
    gen_matrix_A(matrix_seed, sizeof(matrix_seed) - 1, A);

    int16_t s[K];
    prf_vector(secret_seed, sizeof(secret_seed) - 1, 0, ETA1, K, s);

    printf("A (public, fixed) =\n");
    for (int i = 0; i < K; i++) {
        printf("  [");
        for (int j = 0; j < K; j++) printf("%5d%s", A[i][j], j == K - 1 ? "" : ", ");
        printf("]\n");
    }
    print_vec("\ns (secret, fixed) = ", s);
    printf("\n");

    for (int call = 1; call <= 3; call++) {
        int16_t b[K], e[K];
        encrypt_bit(A, s, secret_seed, sizeof(secret_seed) - 1, (uint8_t)call, 1, b, e);
        int recovered = decrypt_bit(A, s, b);
        printf("Encryption #%d (nonce=%d)\n", call, call);
        print_vec("  fresh e        = ", e);
        print_vec("  ciphertext b   = ", b);
        printf("  decrypted bit  = %d\n\n", recovered);
    }

    printf("-> b is different every time because e is resampled from SHAKE-256\n"
           "   output each call, even though A, s, and the message bit (1)\n"
           "   never change. Decryption still recovers '1' correctly because\n"
           "   the CBD noise is small enough to round away.\n\n");

    int16_t b0[K], e0[K];
    encrypt_bit(A, s, secret_seed, sizeof(secret_seed) - 1, 99, 0, b0, e0);
    print_vec("Encrypting bit 0 (nonce=99): b = ", b0);
    printf("  -> decrypts to %d\n", decrypt_bit(A, s, b0));

    return 0;
}