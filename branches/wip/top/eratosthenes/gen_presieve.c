#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>

#define NUMP 22
#define BOUND 512

int main(int argc, char **argv)
{
  uint64_t u64;
  int i, j, k, bk;
  uint32_t primes[NUMP] = {5, 7, 11, 13, 17, 19, 23, 29, 
    31, 37, 41, 43, 47, 53, 59, 61,
    67, 71, 73, 79, 83, 89};

  // for each prime
  for (i = 0; i < NUMP; i++)
  {    
    // for each possible starting position
    printf(" _%u_MASKS_%d[%u] = {\n", primes[i], BOUND, primes[i]);
    for (j = 0; j < primes[i]; j++)
    {
      uint64_t interval[BOUND/64];
      for (bk = 0; bk < BOUND/64; bk++) {
        interval[bk] = 0xffffffffffffffffULL;
      }
        
      // sieve up to the bound, printing each 64-bit word as we fill it
      k = j;

      for (; k < BOUND; k += primes[i])
      {
        interval[k>>6] &= ~(1ULL << k);          
      }

      for (k = 0; k < BOUND/64; k++)
      {
        printf("0x%lxULL,", interval[k]);
      }
      printf("\n");
    }
    printf("\b};\n\n");
  }

  return 0;
}