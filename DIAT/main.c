int main(void)
{
    kb_setup();          /* generates kb_A / kb_s once */
    kb_init_mu_array();  /* fills the mu array with random 0/1 */

    FILE *fp = fopen("prng_output.txt", "w");
    if (!fp) { perror("fopen"); return 1; }

    long total_bits_target = 1000000;  /* set to however many bits you need */
    long bits_written = 0;

    while (bits_written < total_bits_target) {

        kb_maybe_reseed();          /* reseeds on the very first iteration,
                                        then again every 2*KB_LFSR_N bits */

        int bit = your_lfsr_step(); /* <-- replace with your team's real
                                        LFSR step/clock function */
        fputc(bit ? '1' : '0', fp);

        bits_written++;
        kb_bits_since_reseed++;
    }

    fclose(fp);
    printf("done: wrote %ld bits\n", bits_written);
    return 0;
}