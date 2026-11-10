/* P-1 reach on the semiprimes files: the fraction of inputs n = p q where
   p-1 or q-1 is (B1, B2)-smooth (every prime power q^e || p-1 has q^e <= B1
   except at most one prime q in (B1, B2]; micropm1's stage 1 uses prime
   powers <= B1... this is the textbook criterion), for a grid of B1 and
   B2 / B1.  Compare with pm1_eval's measured rates at B2 = 25 B1.
   usage: pm1_smooth bits [num]     (run from the yafu root)
   gcc -O2 pm1_smooth.c -o pm1_smooth */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>

/* largest prime power factor (as q^e) other than the largest prime, and
   the largest prime */
static void split(uint64_t x, uint64_t *big, uint64_t *rest)
{
    uint64_t p, b = 1, r = 1;
    for (p = 2; p * p <= x; p += (p == 2 ? 1 : 2)) {
        if (x % p)
            continue;
        uint64_t pe = 1;
        while (x % p == 0) { x /= p; pe *= p; }
        // a prime power p^e: smooth part needs p^e <= B1 unless p is the
        // single large prime (then e must be 1)
        if (pe == p && p > b) { if (b > r) r = b; b = p; }
        else if (pe > r) r = pe;
    }
    if (x > 1) {
        if (x > b) { if (b > r) r = b; b = x; }
        else if (x > r) r = x;
    }
    *big = b;
    *rest = r;
}

int main(int argc, char **argv)
{
    int bits = argc > 1 ? atoi(argv[1]) : 64;
    int num = argc > 2 ? atoi(argv[2]) : 20000;
    static uint64_t comp[200000], bg[2][200000], rs[2][200000];
    static uint32_t f1[200000], f2[200000];
    static const int b1s[] = { 100, 200, 333, 500, 666, 1000 };
    static const int ms[] = { 25, 50, 100, 200, 400 };
    char fn[64];
    int n = 0, i, a, b;
    sprintf(fn, "semiprimes_%dbit.dat", bits);
    FILE *in = fopen(fn, "r");
    while (n < num && fscanf(in, "%lu,%u,%u", &comp[n], &f1[n], &f2[n]) == 3)
        n++;
    fclose(in);
    for (i = 0; i < n; i++) {
        split(f1[i] - 1, &bg[0][i], &rs[0][i]);
        split(f2[i] - 1, &bg[1][i], &rs[1][i]);
    }
    printf("%d bits: fraction with p-1 or q-1 (B1, m B1)-smooth\n  B1  ", bits);
    for (b = 0; b < 5; b++)
        printf("   m=%-4d", ms[b]);
    printf("\n");
    for (a = 0; a < 6; a++) {
        printf("%5d ", b1s[a]);
        for (b = 0; b < 5; b++) {
            uint64_t B1 = b1s[a], B2 = B1 * ms[b];
            int c = 0;
            for (i = 0; i < n; i++) {
                int k, hit = 0;
                for (k = 0; k < 2; k++)
                    hit |= rs[k][i] <= B1 && bg[k][i] <= B2;
                c += hit;
            }
            printf("  %.4f ", (double)c / n);
        }
        printf("\n");
    }
    return 0;
}
