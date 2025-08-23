#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include <sys/time.h>
#include "../ytools/ytools.h"
#include "batch_factor.h"
#include "gmp.h"

// build line:
// clang -O2 -I. -I../../ytools/ -I../../include -I../../ms_include -I../../top/aprcl -I../../../gmp-install/6.2.0-aocc/include -L../../../gmp-install/6.2.0-aocc/lib ../../ytools/util.c ../../ytools/threadpool.c batch_factor.c process_batch.c tinyecm.c microecm.c prime_sieve.c micropm1.c -o bfact -lm -lgmp -pthread
// 
// usage: bfact lpb fbmax relsfilein relsfileout

int main(int argc, char** argv)
{
    relation_batch_t rb;

    uint64_t lpb = atoi(argv[1]);
    uint64_t pmax = strtoull(argv[2], NULL, 10);
    char infile[80], outfile[80];
    char buf[1024], str1[1024], str2[1024];
    uint32_t fr[32], fa[32], numr = 0, numa = 0;
    mpz_t res1, res2;
    struct timeval start;
    struct timeval stop;
    double ttime;
    uint64_t lcg_state = 42;
    int i;
    uint32_t line = 0;
    uint32_t numfull = 0;

    strcpy(infile, argv[3]);
    strcpy(outfile, argv[4]);

    if (lpb > 40)
    {
        printf("expected large prime bound < 40 bits, saw %lu\n", lpb);
        exit(0);
    }

    lpb = 1ull << lpb;

    rb.num_uecm[0] = 0;
    rb.num_uecm[1] = 0;
    rb.num_uecm[2] = 0;
    rb.num_uecm[3] = 0;
    rb.num_tecm = 0;
    rb.num_tecm2 = 0;
    rb.num_qs = 0;
    rb.num_attempt = 0;
    rb.num_success = 0;
    for (i = 0; i < 8; i++)
    {
        rb.num_abort[i] = 0;
    }
    mpz_init(res1);
    mpz_init(res2);

    printf("initializing relation batch...\n");

    gettimeofday(&start, NULL);
    relation_batch_init(stdout, &rb, pmax, lpb / 2, lpb, lpb, NULL, 1);
    gettimeofday(&stop, NULL);
    ttime = ytools_difftime(&start, &stop);

    printf("init took %1.2f sec, reading input file...\n", ttime);

    FILE* fid = fopen(infile, "r");
    if (fid == NULL)
    {
        printf("could not open %s to read\n", infile);
        exit(0);
    }

    gettimeofday(&start, NULL);

    while (~feof(fid))
    {
        int64_t a;
        uint32_t b;
        char* tok;

        line++;
        char *ptr = fgets(buf, 1024, fid);
        if (ptr == NULL)
            break;

        // printf("buffer read: %s\n", buf);

        tok = strtok(buf, ":");
        if (tok == NULL)
        {
            printf("could not read relation %u, no lfactors token\n", line);
            continue;
        }

        ptr = strchr(tok, ',');
        *ptr = '\0';
        mpz_set_str(res1, tok, 10);
        mpz_set_str(res2, ptr+1, 10);

        // gmp_printf("large factors: %Zd, %Zd\n", res1, res2);

        tok = strtok(NULL, ":");
        if (tok == NULL)
        {
            printf("could not read relation %u, no a/b token\n", line);
            continue;
        }

        sscanf(tok, "%ld,%u", &a, &b);

        tok = strtok(NULL, ":");
        if (tok == NULL)
        {
            printf("could not read relation %u, no rfactors token\n", line);
            continue;
        }
        numr = 0;
        while (strlen(tok) > 0)
        {
            fr[numr++] = strtoul(tok, NULL, 16);
            ptr = strchr(tok, ',');
            if (ptr == NULL)
                break;
            tok = ptr + 1;
        }

        tok = strtok(NULL, ":");
        if (tok == NULL)
        {
            printf("could not read relation %u, no afactors token\n", line);
            continue;
        }

        numa = 0;
        while (strlen(tok) > 0)
        {
            fa[numa++] = strtoul(tok, NULL, 16);
            ptr = strchr(tok, ',');
            if (ptr == NULL)
                break;
            tok = ptr + 1;
        }

        if ((mpz_sgn(res1) > 0) && (mpz_sgn(res2) > 0))
        {
            numfull++;
        }
        else
        {
            relation_batch_add(a, b, fr, numr, res1, fa, numa, res2, &rb);
        }
    }

    gettimeofday(&stop, NULL);
    ttime = ytools_difftime(&start, &stop);
    printf("file parsing took %1.2f sec, found %d fulls, batched %u rels, running batch solve...\n", 
        ttime, numfull, rb.num_relations);

    gettimeofday(&start, NULL);
    relation_batch_run(&rb, &lcg_state);
    gettimeofday(&stop, NULL);

    ttime = ytools_difftime(&start, &stop);
    printf("relation_batch_run took %1.4f sec producing %u relations\n",
        ttime, rb.num_success);

    mpz_clear(res1);
    mpz_clear(res2);
    relation_batch_free(&rb);

    return 0;
}