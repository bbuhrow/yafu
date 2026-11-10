/* tinyrho tests and timings (factor/gmp-ecm/tinyrho.c).
   usage: rho_test mode [num [args]]     (run from the yafu root for files)
     mode f64  : semiprimes_32..64bit.dat through tinyrho64 (num per file)
     mode f52  : semiprimes_32..52bit.dat through tinyrho52_x8_list and
                 tinyrho64 (same inputs; time per input, speedup)
     mode u    : unbalanced n = p q, p of pb bits, n of nb bits (GMP random
                 primes), through tinyrho128 and, for nb <= 104,
                 tinyrho104_x8_list; args: pb0 pb1 nb (default 16 30 100)
     mode f128 : semiprimes_32..64bit.dat as 128-bit inputs through
                 tinyrho104_x8_list (n < 2^60 goes to the scalar code, so
                 this checks the x8 104 path only for 60..64 bits)
   Every returned factor is checked (divides n, nontrivial).
   build: arith/tinyrho/build.sh */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <time.h>
#include "gmp.h"
#include "tinyrho.h"

static double now(void)
{
    struct timespec t;
    clock_gettime(CLOCK_MONOTONIC, &t);
    return t.tv_sec + 1e-9 * t.tv_nsec;
}

static int read_file(int bits, uint64_t *comp, uint32_t *f1, uint32_t *f2, int max)
{
    char fn[64];
    int i = 0;
    FILE *in;
    sprintf(fn, "semiprimes_%dbit.dat", bits);
    in = fopen(fn, "r");
    if (!in)
        return 0;
    while (i < max && fscanf(in, "%lu,%u,%u", &comp[i], &f1[i], &f2[i]) == 3)
        i++;
    fclose(in);
    return i;
}

static int ok64(uint64_t n, uint64_t f)
{
    return f > 1 && f < n && n % f == 0;
}

static gmp_randstate_t grs;

static void rprime(mpz_t p, int bits)
{
    do {
        mpz_urandomb(p, grs, bits);
        mpz_setbit(p, bits - 1);
        mpz_nextprime(p, p);
    } while ((int)mpz_sizeinbase(p, 2) != bits);
}

static void to2(const mpz_t x, uint64_t *lo, uint64_t *hi)
{
    mpz_t t;
    mpz_init(t);
    mpz_tdiv_r_2exp(t, x, 64);
    *lo = mpz_get_ui(t);
    mpz_tdiv_q_2exp(t, x, 64);
    *hi = mpz_get_ui(t);
    mpz_clear(t);
}

static int ok128(const mpz_t n, uint64_t lo, uint64_t hi)
{
    mpz_t f;
    int r;
    mpz_init(f);
    mpz_set_ui(f, hi);
    mpz_mul_2exp(f, f, 64);
    mpz_add_ui(f, f, lo);
    r = mpz_cmp_ui(f, 1) > 0 && mpz_cmp(f, n) < 0 && mpz_divisible_p(n, f);
    mpz_clear(f);
    return r;
}

int main(int argc, char **argv)
{
    const char *mode = argc > 1 ? argv[1] : "f64";
    int num = argc > 2 ? atoi(argv[2]) : 10000;
    static uint64_t comp[200000], fo[400000], q104[400000];
    static uint32_t f1[200000], f2[200000];
    int bits, i;
    gmp_randinit_default(grs);
    gmp_randseed_ui(grs, 1);

    if (!strcmp(mode, "f64") || !strcmp(mode, "f52")) {
        int top = strcmp(mode, "f52") ? 64 : 52;
        printf("bits |  n    | rho64 us  ok  | %s\n", top == 52 ? "x8 us  ok  its/input | speedup" : "");
        for (bits = 32; bits <= top; bits += 2) {
            int n = read_file(bits, comp, f1, f2, num), good = 0, good8 = 0;
            double t0, t1, t2;
            if (n == 0) { printf("no file for %d bits\n", bits); continue; }
            t0 = now();
            for (i = 0; i < n; i++)
                good += ok64(comp[i], tinyrho64(comp[i], 1, 0));
            t1 = now();
            printf("%4d | %5d | %8.2f %5d", bits, n, 1e6 * (t1 - t0) / n, good);
            if (top == 52) {
                tinyrho_list_iterations();
                t1 = now();
                tinyrho52_x8_list(comp, fo, n, 0);
                t2 = now();
                for (i = 0; i < n; i++)
                    good8 += ok64(comp[i], fo[i]);
                printf(" | %7.2f %5d %8.0f | %.2fx", 1e6 * (t2 - t1) / n, good8,
                    (double)tinyrho_list_iterations() / n, (t1 - t0) / (t2 - t1));
            }
            printf("\n");
        }
        return 0;
    }

    if (!strcmp(mode, "f128")) {
        for (bits = 56; bits <= 64; bits += 2) {
            int n = read_file(bits, comp, f1, f2, num), good = 0;
            double t0, t1;
            for (i = 0; i < n; i++) {
                q104[2 * i] = comp[i] & 0xfffffffffffffULL;
                q104[2 * i + 1] = comp[i] >> 52;
            }
            t0 = now();
            tinyrho104_x8_list(q104, fo, n, 0);
            t1 = now();
            for (i = 0; i < n; i++)
                good += fo[2 * i + 1] == 0 && ok64(comp[i], fo[2 * i]);
            printf("%d bits: x8 104 %.2f us/input, %d of %d\n", bits, 1e6 * (t1 - t0) / n, good, n);
        }
        return 0;
    }

    if (!strcmp(mode, "u")) {
        int pb0 = argc > 3 ? atoi(argv[3]) : 16, pb1 = argc > 4 ? atoi(argv[4]) : 30;
        int nb = argc > 5 ? atoi(argv[5]) : 100, pb;
        mpz_t p, q, n, *ns = malloc(num * sizeof(mpz_t));
        mpz_inits(p, q, n, NULL);
        for (i = 0; i < num; i++)
            mpz_init(ns[i]);
        printf("n %d bits | p bits | rho128 us  ok | x8 104 us  ok | speedup\n", nb);
        for (pb = pb0; pb <= pb1; pb += 2) {
            int good = 0, good8 = 0;
            double t0, t1, t2;
            uint64_t nn[2], ff[2];
            for (i = 0; i < num; i++) {
                do {
                    rprime(p, pb);
                    rprime(q, nb - pb);
                    mpz_mul(ns[i], p, q);
                } while ((int)mpz_sizeinbase(ns[i], 2) != nb);
            }
            t0 = now();
            for (i = 0; i < num; i++) {
                to2(ns[i], &nn[0], &nn[1]);
                tinyrho128(nn, ff, 1, 0);
                good += ok128(ns[i], ff[0], ff[1]);
            }
            t1 = now();
            printf("%9d | %6d | %9.2f %5d", nb, pb, 1e6 * (t1 - t0) / num, good);
            if (nb <= 104) {
                for (i = 0; i < num; i++) {
                    to2(ns[i], &nn[0], &nn[1]);
                    q104[2 * i] = nn[0] & 0xfffffffffffffULL;
                    q104[2 * i + 1] = (nn[0] >> 52) | (nn[1] << 12);
                }
                t1 = now();
                tinyrho104_x8_list(q104, fo, num, 0);
                t2 = now();
                for (i = 0; i < num; i++)
                    good8 += ok128(ns[i], fo[2 * i], fo[2 * i + 1]);
                printf(" | %9.2f %5d | %.2fx", 1e6 * (t2 - t1) / num, good8, (t1 - t0) / (t2 - t1));
            }
            printf("\n");
        }
        return 0;
    }
    printf("unknown mode %s\n", mode);
    return 1;
}
