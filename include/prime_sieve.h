/*--------------PRIME SIEVE RELATED DECLARATIONS ---------------------*/

/* many separate places in the code need a list
   of primes. Include such a list pregenerated, with
   the differences between primes stored */
#include "yafu.h"
#include "qs.h"
#include "factor.h"
#include "ytools.h"
#include <stdint.h>

#define PRECOMPUTED_PRIME_BOUND 100000
#define PRECOMPUTED_NUM_PRIMES 9592
extern const uint8_t_t prime_delta[PRECOMPUTED_NUM_PRIMES];

typedef struct {
    uint32_t_t p;
    uint32_t_t r;
} prime_aux_t;

typedef struct {
    uint32_t_t num_aux;
    prime_aux_t *aux;
    uint8_t_t *sieve;
    uint32_t_t curr_block;
    uint32_t_t curr_off;
} prime_sieve_t;

void init_prime_sieve(prime_sieve_t *sieve,
    uint32_t_t min_prime,
    uint32_t_t max_prime);

void free_prime_sieve(prime_sieve_t *sieve);

uint32_t get_next_prime(prime_sieve_t *sieve);

typedef struct {
    uint32_t_t *list;
    uint32_t_t num_primes;
} prime_list_t;

void fill_prime_list(prime_list_t *prime_list,
    uint32_t_t max_size,
    uint32_t_t max_prime);
