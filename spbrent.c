#include <stdlib.h>
#include <stdint.h>
#include <stdio.h>
#include <math.h>
#include <sys/types.h>
#include <sys/time.h>
#include <time.h>

typedef int32_t int32;
typedef int64_t int64;
typedef uint16_t uint16;
typedef uint32_t uint32;
typedef uint64_t uint64;

#define MIN(a,b) ((a) < (b)? (a) : (b))

uint64 bingcd(uint64 x, uint64 y) 
{
    if (y) {
        y >>= __builtin_ctzl(y);
        while (x != y)
            if (x < y)
                y -= x, y >>= __builtin_ctzl(y);
            else
                x -= y, x >>= __builtin_ctzl(x);
    }
    return x;
}

__inline uint64 addmod(uint64 x, uint64 y, uint64 n)
{
    uint64_t t = x-n;
    x += y;
    asm ("add %2, %1\n\t"
         "cmovc %1, %0\n\t"
         :"+r" (x), "+&r" (t)
         :"r" (y)
         :"cc"
        );
    return x;
}

__inline uint64 u64div(uint64 c, uint64 n)
{
    __asm__("divq %4"
        : "=a"(c), "=d"(n)
        : "1"(c), "0"(0), "r"(n));

    return n;
}

__inline uint64 mulredc(uint64 x, uint64 y, uint64 n, uint64 nhat)
{
    if (n & 0x8000000000000000)
    {
        __asm__(
          "mulq %2	\n\t"
          "movq %%rax, %%r10		\n\t"
          "movq %%rdx, %%r11		\n\t"
          "movq $0, %%r12 \n\t"
          "mulq %3 \n\t"
          "mulq %4 \n\t"
          "addq %%r10, %%rax \n\t"        
          "adcq %%r11, %%rdx \n\t"
          "cmovae %4, %%r12 \n\t"  
          "xorq %%rax, %%rax \n\t"        
          "subq %4, %%rdx \n\t"
          "cmovc %%r12, %%rax \n\t"
          "addq %%rdx, %%rax \n\t"
          : "=&a"(x)
          : "0"(x), "r"(y), "r"(nhat), "r"(n)
          : "rdx", "r10", "r11", "r12", "cc");
    }
    else
    {
        __asm__(
          "mulq %2	\n\t"
          "movq %%rax, %%r10		\n\t"
          "movq %%rdx, %%r11		\n\t"
          "mulq %3 \n\t"
          "mulq %4 \n\t"
          "addq %%r10, %%rax \n\t"
          "adcq %%r11, %%rdx \n\t"
          "movq $0, %%rax \n\t"
          "subq %4, %%rdx \n\t"
          "cmovc %4, %%rax \n\t"
          "addq %%rdx, %%rax \n\t"
          : "=&a"(x)
          : "0"(x), "r"(y), "r"(nhat), "r"(n)
          : "rdx", "r10", "r11", "cc");
      
    }
    return x;
}

__inline uint64 mulredc63(uint64 x, uint64 y, uint64 n, uint64 nhat)
{
    __asm__(
        "mulq %2	\n\t"
        "movq %%rax, %%r10		\n\t"
        "movq %%rdx, %%r11		\n\t"
        "mulq %3 \n\t"
        "mulq %4 \n\t"
        "addq %%r10, %%rax \n\t"
        "adcq %%r11, %%rdx \n\t"
        "xorq %%rax, %%rax \n\t"
        "subq %4, %%rdx \n\t"
        "cmovc %4, %%rax \n\t"
        "addq %%rdx, %%rax \n\t"
        : "=a"(x)
        : "0"(x), "r"(y), "r"(nhat), "r"(n)
        : "rdx", "r10", "r11", "cc");

    return x;
}

// for comparison... arbooker's C-version mulredc from http://www.mersenneforum.org/showthread.php?t=22525
// #define mulredc(x,y,n,nbar) (t.d=(u128)(x)*(y),t.d+=(t.l[0]*nbar)*(u128)n,\
// t.l[1]-=n,t.l[1]+(n&((i64)t.l[1]>>63)))

uint64 spbrent(uint64 N, uint64 c, int imax)
{
    /*
    run pollard's rho algorithm on n with Brent's modification,
    returning the first factor found in f, or else 0 for failure.
    use f(x) = x^2 + c
    see, for example, bressoud's book.
    */
    uint64 x, y, q, g, ys, t1, f, nhat;
    uint32 i = 0, k, r, m;
    int it;

    // start out checking gcd fairly often
    r = 1;

    // urder 48 bits, don't defer gcd quite as long
    i = __builtin_clzl(N);
    if (i > 20)
        m = 32;
    else if (i > 16)
        m = 160;
    else if (i > 3)
        m = 256;   
    else
        m = 384;

    it = 0;
    q = 1;    
    g = 1;

    x = (((N + 2) & 4) << 1) + N; // here x*a==1 mod 2**4
    x *= 2 - N * x;               // here x*a==1 mod 2**8
    x *= 2 - N * x;               // here x*a==1 mod 2**16
    x *= 2 - N * x;               // here x*a==1 mod 2**32         
    x *= 2 - N * x;               // here x*a==1 mod 2**64
    nhat = (uint64)0 - x;

    // Montgomery representation of c
    c = u64div(c, N);
    y = c;    
    
    do
    {
        x = y;
        for (i = 0; i <= r; i++)
        {
            y = mulredc63(y, y+c, N, nhat);
        }

        k = 0;
        do
        {
            ys = y;
            for (i = 1; i <= MIN(m, r - k); i++)
            {
                y = mulredc63(y, y+c, N, nhat);
                t1 = x > y ? y - x + N : y - x;
                q = mulredc63(q, t1, N, nhat);
            }

            g = bingcd(N,q);
            k += m;
            it++;

            if (it>imax)
            {
                f = 0;
                goto done;
            }

        } while ((k<r) && (g == 1));
        r *= 2;
    } while (g == 1);

    if (g == N)
    {
        //back track
        do
        {
            ys = mulredc63(ys, ys+c, N, nhat);
            t1 = x > ys ? ys - x + N : ys - x;
            g = bingcd(N,t1);
        } while (g == 1);

        if (g == N)
        {
            f = 0;
        }
        else
        {
            f = g;
        }
    }
    else
    {
        f = g;
    }

    done:
    return f;
}

double my_difftime (struct timeval * start, struct timeval * end)
{
    double secs;
    double usecs;

	if (start->tv_sec == end->tv_sec) {
		secs = 0;
		usecs = end->tv_usec - start->tv_usec;
	}
	else {
		usecs = 1000000 - start->tv_usec;
		secs = end->tv_sec - (start->tv_sec + 1);
		usecs += end->tv_usec;
		if (usecs >= 1000000) {
			usecs -= 1000000;
			secs += 1;
		}
	}
	
	return secs + usecs / 1000000.;
}

uint64 lr_pmod64(uint64 b, uint64 e, uint64 m)
{
    int sz_b;
    int i, j;
    uint64 x;
    uint64 s;
    uint64 one;
    uint64 mhat;
    uint64 ee = e;
    uint64 t = 1ULL << 63;
    
    x = (((m + 2) & 4) << 1) + m; // here x*a==1 mod 2**4
    x *= 2 - m * x;               // here x*a==1 mod 2**8
    x *= 2 - m * x;               // here x*a==1 mod 2**16
    x *= 2 - m * x;               // here x*a==1 mod 2**32         
    x *= 2 - m * x;               // here x*a==1 mod 2**64
    mhat = (uint64)0 - x;
    
    x = b;
    s = 1;
    
    // Montgomery representation
    s = one = u64div(s, m);
    x = u64div(x, m);
    
    while ((t & e) == 0)
      t >>= 1;
    
    // L-R exponentiation
    while (t)
    {
        s = mulredc(s, s, m, mhat);
        if (e & t)
            s = mulredc(x, s, m, mhat);
        t >>= 1;
    }

    s = mulredc(s, 1, m, mhat);
    //printf("%lu ^ %lu mod %lu = %lu (one = %lu, mhat = %lu)\n", b, e, m, s, one, mhat);
    return s;
}

uint64 rl_pmod64(uint64 b, uint64 e, uint64 m)
{
    int sz_b;
    int i, j;
    uint64 x;
    uint64 s;
    uint64 one;
    uint64 mhat;
    uint64 ee = e;
    
    x = (((m + 2) & 4) << 1) + m; // here x*a==1 mod 2**4
    x *= 2 - m * x;               // here x*a==1 mod 2**8
    x *= 2 - m * x;               // here x*a==1 mod 2**16
    x *= 2 - m * x;               // here x*a==1 mod 2**32         
    x *= 2 - m * x;               // here x*a==1 mod 2**64
    mhat = (uint64)0 - x;
    
    x = 1;
    s = b;
    
    // Montgomery representation
    x = one = u64div(x, m);
    s = u64div(s, m);
    
    // R-L exponentiation
    while (e)
    {
        if (e & 1)
          x = mulredc(x, s, m, mhat);
        s = mulredc(s, s, m, mhat);
        e >>= 1;
    }

    x = mulredc(x, 1, m, mhat);
    //printf("%lu ^ %lu mod %lu = %lu\n", b, ee, m, x);
    return x;
}

uint64 lrwin_pmod64(uint64 b, uint64 e, uint64 m)
{
    int sz_b;
    int i, j;
    int bit; // = mpz_sizeinbase(e, 2) - 1;
    int mask;
    int nsqr, nmul;
    int bstr;
    int k; // = get_winsize();
    //mpz_t g[256];
    //mpz_t acc;

    mask = 0;
    for (j = 0; j < k; j++)
    {
        mask = (mask << 1) | 1;
    }

    //for (i = 0; i < 256; i++)
    //    mpz_init(g[i]);

    // setup accumulator: 1
    //mpz_set_ui(x, 1);

    // precomputations, b^i for 0 <= i < 2^k
    //mpz_set(g[1], b);
    //for (i = 2; i < (1 << k); i++)
    //{
    //    mpz_mul(g[i], g[i - 1], b);
    //    mpz_mod(g[i], g[i], m);
    //}

    // L-R windowed exponentiation.  Scan the exponent bit-vector
    // backward instead of flipping and shifting it.
    while (bit >= 0)
    {
        if (bit < k)
        {
            // grab the last bits of the exponent.
            // accommodates exponent lengths not divisible
            // by the window size
            mask = 0x0;
            for (j = 0; j < (bit + 1); j++)
            {
                //mpz_mul(x, x, x);
                //mpz_mod(x, x, m);
                mask = (mask << 1) | 1;
            }

            //bstr = e->_mp_d[0] & mask;
        }
        else
        {
            // grab the next k bits of the exponent.
            //bstr = get_bitwin(e, bit, k, mask);
            for (j = 0; j < k; j++)
            {
                //mpz_mul(x, x, x);
                //mpz_mod(x, x, m);
            }
        }

        if (bstr > 0)
        {
            //mpz_mul(x, g[bstr], x);
            //mpz_mod(x, x, m);
        }

        bit -= k;

    }

    //for (i = 0; i < 256; i++)
    //    mpz_clear(g[i]);

    return;
}

int main(int argc, char **argv)
{
    FILE *in;
    uint64 *comp, f64;
    uint32 *f1;
    uint32 *f2;
    double t_time;
    clock_t start, stop;
    int i,j,k,num,correct;
    struct timeval gstart;	
    struct timeval gstop;
    int nf;
    int num_files;
    char filenames[15][80];
    
    comp = (uint64 *)malloc(200000 * sizeof(uint64));
    f1 = (uint32 *)malloc(200000 * sizeof(uint32));
    f2 = (uint32 *)malloc(200000 * sizeof(uint32));

    i = 0;
    //strcpy(filenames[i++], "pseudoprimes_36bit.dat");
    //strcpy(filenames[i++], "pseudoprimes_38bit.dat");
    //strcpy(filenames[i++], "pseudoprimes_40bit.dat");
    //strcpy(filenames[i++], "pseudoprimes_42bit.dat");
    //strcpy(filenames[i++], "pseudoprimes_44bit.dat");
    //strcpy(filenames[i++], "pseudoprimes_46bit.dat");
    //strcpy(filenames[i++], "pseudoprimes_48bit.dat");
    //strcpy(filenames[i++], "pseudoprimes_50bit.dat");
    //strcpy(filenames[i++], "pseudoprimes_52bit.dat");
    //strcpy(filenames[i++], "pseudoprimes_54bit.dat");
    //strcpy(filenames[i++], "pseudoprimes_56bit.dat");
    //strcpy(filenames[i++], "pseudoprimes_58bit.dat");
    strcpy(filenames[i++], "pseudoprimes_60bit.dat");
    strcpy(filenames[i++], "pseudoprimes_62bit.dat");
    strcpy(filenames[i++], "pseudoprimes_63bit.dat");
    strcpy(filenames[i++], "pseudoprimes_64bit.dat");
    num_files = 4;
    
    
    for (nf = 0; nf < num_files; nf++)
    {
        uint64 c[3] = {1, 2, 3};
        in = fopen(filenames[nf], "r");

        start = clock();
        i = 0;

        //read in everything
        while (!feof(in))
        {
            fscanf(in, "%lu,%u,%u", comp + i, f1 + i, f2 + i);
            i++;
        }
        num = i;
        
        fclose(in);
        stop = clock();
        t_time = (double)(stop - start) / (double)CLOCKS_PER_SEC;
        printf("%d inputs read in %2.4f sec\n", num, t_time);

        num = 100000;
        gettimeofday(&gstart, NULL);

        correct = 0;
        k = 0;
        for (i = 0; i < num; i++)
        {
          
            #if 1
            int k;
            uint64 x = 0;
            for (k = 0; k < 100; k++)
            {
                //f64 = lr_pmod64(k+2, comp[i]-1, comp[i]);
                f64 = rl_pmod64(k+2, comp[i]-1, comp[i]);
                if (f64 == 1)
                {
                  correct++;
                  printf("%lu is a base-%d fermat pseudoprime\n", comp[i], k+2);
                }
            }
            #else
            int p;
            for (p = 0; p < 3; p++)
            {
                f64 = spbrent(comp[i], c[p], 8192);
                if ((f64 == f1[i]) || (f64 == f2[i]))
                {
                    correct++;
                    //printf("%lu = %lu * %lu\n", comp[i], f64, comp[i] / f64);
                    if (p > 0)
                      k++;
                    break;
                }
            }    
            #endif
        }

        gettimeofday(&gstop, NULL);
        t_time = my_difftime(&gstart, &gstop);

        printf("rho got %d of %d correct in file %s in %2.2f sec (%d with backup polynomial)\n", 
          correct, num, filenames[nf], t_time, k);
        printf("percent correct = %.2f\n", 100.0*(double)correct / (double)num);
        printf("average time per input = %1.4f ms\n", 1000 * t_time / (double)num);
    }
    
    free(comp);
    free(f1);
    free(f2);

    return 0;
}
