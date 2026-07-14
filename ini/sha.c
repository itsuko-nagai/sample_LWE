/* sha3.c / sha3.h combined -- public-domain style compact Keccak/SHA-3
 * implementation (same sponge construction Kyber's fips202.c uses).
 * Provides SHAKE128 / SHAKE256 extendable-output functions.
 */
#include <stdint.h>
#include <stddef.h>
#include <string.h>

#define KECCAKF_ROUNDS 24
#define ROTL64(x, y) (((x) << (y)) | ((x) >> (64 - (y))))

typedef struct {
    union {
        uint8_t b[200];
        uint64_t q[25];
    } st;
    int pt, rsiz, mdlen;
} sha3_ctx_t;

static const uint64_t keccakf_rndc[24] = {
    0x0000000000000001ULL, 0x0000000000008082ULL, 0x800000000000808aULL,
    0x8000000080008000ULL, 0x000000000000808bULL, 0x0000000080000001ULL,
    0x8000000080008081ULL, 0x8000000000008009ULL, 0x000000000000008aULL,
    0x0000000000000088ULL, 0x0000000080008009ULL, 0x000000008000000aULL,
    0x000000008000808bULL, 0x800000000000008bULL, 0x8000000000008089ULL,
    0x8000000000008003ULL, 0x8000000000008002ULL, 0x8000000000000080ULL,
    0x000000000000800aULL, 0x800000008000000aULL, 0x8000000080008081ULL,
    0x8000000000008080ULL, 0x0000000080000001ULL, 0x8000000080008008ULL
};
static const int keccakf_rotc[24] = {
    1,  3,  6,  10, 15, 21, 28, 36, 45, 55, 2,  14,
    27, 41, 56, 8,  25, 43, 62, 18, 39, 61, 20, 44
};
static const int keccakf_piln[24] = {
    10, 7,  11, 17, 18, 3, 5,  16, 8,  21, 24, 4,
    15, 23, 19, 13, 12, 2, 20, 14, 22, 9,  6,  1
};

static void sha3_keccakf(uint64_t st[25])
{
    int i, j, r;
    uint64_t t, bc[5];

    for (r = 0; r < KECCAKF_ROUNDS; r++) {
        /* Theta */
        for (i = 0; i < 5; i++)
            bc[i] = st[i] ^ st[i + 5] ^ st[i + 10] ^ st[i + 15] ^ st[i + 20];
        for (i = 0; i < 5; i++) {
            t = bc[(i + 4) % 5] ^ ROTL64(bc[(i + 1) % 5], 1);
            for (j = 0; j < 25; j += 5)
                st[j + i] ^= t;
        }
        /* Rho + Pi */
        t = st[1];
        for (i = 0; i < 24; i++) {
            j = keccakf_piln[i];
            bc[0] = st[j];
            st[j] = ROTL64(t, keccakf_rotc[i]);
            t = bc[0];
        }
        /* Chi */
        for (j = 0; j < 25; j += 5) {
            for (i = 0; i < 5; i++)
                bc[i] = st[j + i];
            for (i = 0; i < 5; i++)
                st[j + i] ^= (~bc[(i + 1) % 5]) & bc[(i + 2) % 5];
        }
        /* Iota */
        st[0] ^= keccakf_rndc[r];
    }
}

static void sha3_init(sha3_ctx_t *c, int mdlen)
{
    memset(c, 0, sizeof(*c));
    c->mdlen = mdlen;
    c->rsiz = 200 - 2 * mdlen;
}

static void sha3_update(sha3_ctx_t *c, const void *data, size_t len)
{
    size_t i;
    int j = c->pt;
    for (i = 0; i < len; i++) {
        c->st.b[j++] ^= ((const uint8_t *)data)[i];
        if (j >= c->rsiz) {
            sha3_keccakf(c->st.q);
            j = 0;
        }
    }
    c->pt = j;
}

/* switch context from absorbing to squeezing (SHAKE domain sep = 0x1F) */
static void shake_xof(sha3_ctx_t *c)
{
    c->st.b[c->pt] ^= 0x1F;
    c->st.b[c->rsiz - 1] ^= 0x80;
    sha3_keccakf(c->st.q);
    c->pt = 0;
}

static void shake_out(sha3_ctx_t *c, void *out, size_t len)
{
    size_t i;
    int j = c->pt;
    for (i = 0; i < len; i++) {
        if (j >= c->rsiz) {
            sha3_keccakf(c->st.q);
            j = 0;
        }
        ((uint8_t *)out)[i] = c->st.b[j++];
    }
    c->pt = j;
}

void shake128(const uint8_t *in, size_t inlen, uint8_t *out, size_t outlen)
{
    sha3_ctx_t c;
    sha3_init(&c, 16);          /* rate = 200-32 = 168 bytes, SHAKE128 rate */
    sha3_update(&c, in, inlen);
    shake_xof(&c);
    shake_out(&c, out, outlen);
}

void shake256(const uint8_t *in, size_t inlen, uint8_t *out, size_t outlen)
{
    sha3_ctx_t c;
    sha3_init(&c, 32);          /* rate = 200-64 = 136 bytes, SHAKE256 rate */
    sha3_update(&c, in, inlen);
    shake_xof(&c);
    shake_out(&c, out, outlen);
}