static void kb_generate_bits(int mu, uint8_t *bitbuf, size_t nbits)
{
    if (!kb_ready) kb_setup();
    if (nbits > KB_MAX_BITS) nbits = KB_MAX_BITS;

    uint8_t noise_seed[KB_SEED_BYTES];
    kb_get_random_bytes(noise_seed, KB_SEED_BYTES);
    kb_polyvec e;
    kb_polyvec_noise(noise_seed, 0, KB_ETA, &e);   /* internally uses shake256 */

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