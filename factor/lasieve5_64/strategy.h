#include <stdarg.h>
#include <stdio.h>

typedef struct {
  u32_t **stindex; /* nr of strategy to use for (n0,n1) */
  u32_t **stlist;  /* list of strategies */
  u16_t mc[2];     /* maximal number of bits for composites */
  unsigned char **bit;
} strat_t;

void print_strategy(strat_t s);
void read_strategy(strat_t *s, u16_t *maxcomp, char *basename, u16_t *maxpr);
int cofactorisation(strat_t *st, mpz_t **large_primes, mpz_t *large_factors,
                      u16_t *max_primebits, u32_t *nlp, mpz_t *FBb_sq,
                      mpz_t *FBb_cu);
void print_strategy_stat();
