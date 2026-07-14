int main(void)
{
    kb_setup();
    kb_init_mu_array();

    FILE *fp = fopen("stream.txt", "w");
    if (!fp) { perror("fopen"); return 1; }

    for (kb_mu_index = 0; kb_mu_index < KB_MU_COUNT; kb_mu_index++) {
        uint8_t seed_bits[KB_SEED_BITS];
        kb_generate_bits(kb_mu_array[kb_mu_index], seed_bits, KB_SEED_BITS);
        LFSR_seed(seed_bits, KB_SEED_BITS);   /* your team's real LFSR seed call */

        g_total_bits_written = 0;
        g_stop_now = 0;

        LFSR_call(fp);
    }

    fclose(fp);
    printf("done\n");
    return 0;
}