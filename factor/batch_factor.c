/*--------------------------------------------------------------------
This source distribution is placed in the public domain by its author,
Jason Papadopoulos. You may use it for any purpose, free of charge,
without having to notify anyone. I disclaim any responsibility for any
errors.

Optionally, please be nice and tell me if you find this source to be
useful. Again optionally, if you add to the functionality present here
please consider making those additions public too, so that others may 
benefit from your work.

$Id: batch_factor.c 638 2011-09-11 15:31:19Z jasonp_sf $
--------------------------------------------------------------------*/

#include <batch_factor.h>
#include <stdint.h>
#include "monty.h"
#include "prime_sieve.h"
#include "yafu_ecm.h"

/*------------------------------------------------------------------

Jasonp's batch GCD for relation factoring, using R. Gerbicz's 
balanced remainder tree idea. See below for the thread where
this idea first came up for reference.

*/



/*
R. Gerbicz: https://www.mersenneforum.org/showpost.php?p=521621&postcount=144

Nice code, now understand more of this, at least for me it was new, so:
Set S = p[0] * p[1] * ... * p[m - 1] where p[i] is the i - th prime and also
using a product tree get Z = r[0] * r[1] * ... * r[n - 1], where you want to
get only the smooth parts of r[].
Then using a remainder tree get v[i] = S mod r[i].
Trivially the i - th number's smooth part is gcd(r[i],v[i]) !

What I don't understand is that you could make a (balanced) remainder tree
in an explicit way, instead of your recursion. With that you'd compute the
product tree for Z only once, not computing multiple times the same subproducts
in the tree while you are doing the remainder tree algorithm.
But you need more memory for that, by a factor of log2(n).

https://www.mersenneforum.org/showpost.php?p=521656&postcount=147
Then a "hybrid" method would work here, build up the whole subtree in memory
if you can hold it (otherwise do recursion), surely at
depth=tree height-6 you should be able to do that:
(where size is the size of integer in the root)
to hold the subtree in memory you need
size/64*(height-6)<size (because even height<64 is true), and that's not much,
you need more memory to do the multiplication/division in the root of the tree.

*/
typedef struct 
{
    mpz_t prod;
    uint32_t low;
    uint32_t high;
    uint32_t complete;
    int left_id;
    int right_id;
} bintree_element_t;

typedef struct
{
    bintree_element_t* nodes;
    uint32_t size;
    uint32_t alloc;
} bintree_t;

uint32_t getNode(bintree_t* tree, uint32_t low, uint32_t high)
{
    uint32_t nodeid = 0;
    bintree_element_t* node = &tree->nodes[nodeid];

    //printf("commencing tree search for (%u:%u)\n", low, high);
    while (node != NULL)
    {
        //printf("node is now %d: (%u:%u)\n", nodeid, node->low, node->high);
        if ((node->low == low) && (node->high == high))
        {
            return nodeid;
        }
        else 
        {
            if (node->left_id == -1)
                return -1;

            if (high <= tree->nodes[node->left_id].high)
            {
                //printf("taking left path\n");
                nodeid = node->left_id;
                node = &tree->nodes[node->left_id];
            }
            else
            {
                if (node->right_id == -1)
                    return -1;

                //printf("taking right path\n");
                nodeid = node->right_id;
                node = &tree->nodes[node->right_id];
            }
        }
    }

    return -1;
}

void addNode(bintree_t* tree, int id, int side, uint32_t low, uint32_t high, mpz_t prod)
{
    uint32_t nodeid = 0;
    bintree_element_t* node = &tree->nodes[id];

    if (tree->size == tree->alloc)
    {
        //printf("growing the tree to size %u\n", tree->alloc + 16);
        tree->alloc *= 2;
        tree->nodes = (bintree_element_t*)xrealloc(tree->nodes,
            tree->alloc * sizeof(bintree_element_t));
    }

    if ((side == 0) && (node->left_id != -1))
    {
        printf("node %d already has a left child\n", id);
        exit(1);
    }
    if ((side == 1) && (node->right_id != -1))
    {
        printf("node %d already has a right child\n", id);
        exit(1);
    }

    //printf("adding node %d as %d-child of %d: (%u:%u)\n", tree->size, side, id, low, high);
    tree->nodes[tree->size].low = low;
    tree->nodes[tree->size].high = high;
    mpz_init(tree->nodes[tree->size].prod);
    if (prod != NULL)
        mpz_set(tree->nodes[tree->size].prod, prod);
    tree->nodes[tree->size].left_id = -1;
    tree->nodes[tree->size].right_id = -1;
    tree->nodes[tree->size].complete = 0;
    if (side == 0)
        tree->nodes[id].left_id = tree->size;
    else
        tree->nodes[id].right_id = tree->size;
    
    //printf("added new node %d and assigned to node %d side %d\n", tree->size, id, side);
    tree->size++;

    return;
}

#define USE_TREE
#define TREE_CUTOFF 16


#define BREAKOVER_WORDS 50
void multiply_primes(uint32_t first, uint32_t last,
    prime_sieve_t *sieve, mpz_t prod) {

    /* recursive routine to multiply all the elements of a
       list of consecutive primes. The current invocation
       multiplies elements first to last of the list, inclusive.
       The primes are created on demand */

    mpz_t half_prod;
    uint32_t mid = (last + first) / 2;

    /* base case; accumulate a few primes */

    if (last < first + BREAKOVER_WORDS) {
        mpz_set_ui(prod, (unsigned long)get_next_prime(sieve));
        while (++first <= last) {
            mpz_mul_ui(prod, prod,
                (unsigned long)get_next_prime(sieve));
        }
        return;
    }

    /* recursively handle the left and right halves of the list */

    mpz_init(half_prod);
    multiply_primes(first, mid, sieve, prod);
    multiply_primes(mid + 1, last, sieve, half_prod);

    /* multiply them together. We can take advantage of
       fast multiplication since in general the two half-
       products are about the same size*/

    mpz_mul(prod, prod, half_prod);
    mpz_clear(half_prod);
}



/*------------------------------------------------------------------*/
void relation_to_gmp(relation_batch_t *rb,
				uint32_t index, mpz_t out) {

	/* multiply together the unfactored parts of an
	   NFS relation, then convert to an mpz_t */

    int j;
	uint32_t i, nwords;
	cofactor_t *c = rb->relations + index;
	uint32_t *f = rb->factors + c->factor_list_word + 
			c->num_factors_r + c->num_factors_a;
	
    mpz_ptr num1 = rb->f1a, num2 = rb->f2a, num3 = rb->n;

	if (c->lp_r_num_words > 1 && c->lp_a_num_words > 1) {

		/* rational and algebraic parts need to be
		   multiplied together first */

		//mp_clear(&num1);
		//mp_clear(&num2);
		//num1.nwords = c->lp_r_num_words;
		//for (i = 0; i < c->lp_r_num_words; i++)
		//	num1.val[i] = f[i];

        mpz_set_ui(num1, f[0]);
        for (i = 1; i < c->lp_r_num_words; i++)
        {
            mpz_mul_2exp(num1, num1, 32);
            mpz_add_ui(num1, num1, f[i]);
        }

		//num2.nwords = c->lp_a_num_words;
		f += c->lp_r_num_words;
		//for (i = 0; i < c->lp_a_num_words; i++)
		//	num2.val[i] = f[i];
        mpz_set_ui(num2, f[0]);
        for (i = 1; i < c->lp_a_num_words; i++)
        {
            mpz_mul_2exp(num2, num2, 32);
            mpz_add_ui(num2, num2, f[i]);
        }

		mpz_mul(num3, num2, num1);
	}
	else {
		/* no multiply necesary */

		//mp_clear(&num3);
		nwords = c->lp_r_num_words;

		if (c->lp_a_num_words > 1) {
			nwords = c->lp_a_num_words;
			f += c->lp_r_num_words;
		}

		//num3.nwords = nwords;
		//for (i = 0; i < nwords; i++)
		//	num3.val[i] = f[i];

        mpz_set_ui(num3, f[nwords - 1]);
        for (j = nwords - 2; j >= 0; j--)
        {
            mpz_mul_2exp(num3, num3, 32);
            mpz_add_ui(num3, num3, f[j]);
        }
	}

    mpz_set(out, num3);

    return;
}

/*------------------------------------------------------------------*/
void multiply_relations(bintree_t* tree, uint32_t first, uint32_t last,
				relation_batch_t *rb,
				mpz_t prod) {
	mpz_t half_prod;

	/* recursive routine to multiply together a list of
	   relations. The current invocation multiplies relations
	   first to last of the list, inclusive */

	/* base case of recursion */

	if (first == last) {
        //printf("leaf relation %u\n", first);
		relation_to_gmp(rb, first, prod);
		return;
	}

	/* recurse on the left and right halves */

	mpz_init(half_prod);
	if (last == first + 1) {
        //printf("forming product (prod(%u) * prod(%u))\n", first, last);
		relation_to_gmp(rb, first, half_prod);
		relation_to_gmp(rb, last, prod);
	}
	else {
        uint32_t mid = (last + first) / 2;

#ifdef USE_TREE
        uint32_t pnode;
        pnode = getNode(tree, first, last);

        // Paul K's proposed fix, which for me caused segfaults... 
        // need to look closer at the patch.
        //if ((pnode >= 0) && ((mid - first) > TREE_CUTOFF))
        //{
        //    uint32_t node;
        //    node = getNode(tree, first, mid);
        //    if ((node != -1) && (tree->nodes[node].complete))
        //    {
        //        //printf("getting half-product (prod(%u:%u)\n", first, mid);
        //        mpz_set(prod, tree->nodes[node].prod);
        //    }
        //    else
        //    {
        //        //printf("forming half-product (prod(%u:%u)\n", first, mid);
        //        addNode(tree, pnode, 0, first, mid, NULL);
        //        multiply_relations(tree, first, mid, rb, prod);
        //    }
        //}
        //
        //if ((pnode >= 0) && ((last - mid - 1) > TREE_CUTOFF))
        //{
        //    uint32_t node;
        //    node = getNode(tree, mid + 1, last);
        //    if ((node != -1) && (tree->nodes[node].complete))
        //    {
        //        //printf("getting half-product (prod(%u:%u)\n", mid + 1, last);
        //        mpz_set(prod, tree->nodes[node].prod);
        //    }
        //    else
        //    {
        //        //printf("forming half-product (prod(%u:%u)\n", mid + 1, last);
        //        addNode(tree, pnode, 1, mid + 1, last, NULL);
        //        multiply_relations(tree, mid + 1, last, rb, half_prod);
        //    }
        //    //}
        //    //else
        //    //{
        //    //    multiply_relations(tree, first, mid, rb, prod);
        //}
        //else
        //{
        //    multiply_relations(tree, mid + 1, last, rb, half_prod);
        //}


        if ((pnode >= 0) && ((mid - first) > TREE_CUTOFF))
        {
            uint32_t node;
            node = getNode(tree, first, mid);
            if ((node != -1) && (tree->nodes[node].complete))
            {
                //printf("getting half-product (prod(%u:%u)\n", first, mid);
                mpz_set(prod, tree->nodes[node].prod);
            }
            else
            {
                //printf("forming half-product (prod(%u:%u)\n", first, mid);
                addNode(tree, pnode, 0, first, mid, NULL);
                multiply_relations(tree, first, mid, rb, prod);
            }
        
        
            node = getNode(tree, mid + 1, last);
            if ((node != -1) && (tree->nodes[node].complete))
            {
                //printf("getting half-product (prod(%u:%u)\n", mid + 1, last);
                mpz_set(prod, tree->nodes[node].prod);
            }
            else
            {
                //printf("forming half-product (prod(%u:%u)\n", mid + 1, last);
                addNode(tree, pnode, 1, mid + 1, last, NULL);
                multiply_relations(tree, mid + 1, last, rb, half_prod);
            }
        }
        else
        {
            multiply_relations(tree, first, mid, rb, prod);
            multiply_relations(tree, mid + 1, last, rb, half_prod);
        }
#else
        multiply_relations(tree, first, mid, rb, prod);
        multiply_relations(tree, mid + 1, last, rb, half_prod);
#endif
		
	}

	/* multiply the halves */
#ifdef USE_TREE
    if ((last - first) > TREE_CUTOFF)
    {
        uint32_t pnode;
        pnode = getNode(tree, first, last);
        if (pnode == -1)
        {
            printf("could not find parent node for full product (%u:%u)\n", first, last);
            exit(1);
        }

        if (tree->nodes[pnode].complete)
        {
            //printf("getting product (prod(%u:%u)\n", first, last);
            mpz_set(prod, tree->nodes[pnode].prod);
        }
        else
        {
            //printf("forming product (prod(%u:%u)\n", first, last);
            mpz_mul(prod, prod, half_prod);
            
            mpz_set(tree->nodes[pnode].prod, prod);
            tree->nodes[pnode].complete = 1;
        }
    }
    else
    {
        mpz_mul(prod, prod, half_prod);
    }
#else
    mpz_mul(prod, prod, half_prod);
#endif

	mpz_clear(half_prod);
}

/*------------------------------------------------------------------*/

uint64_t pow2m(uint64_t b, uint64_t n)
{
    // compute 2^b mod n
    uint64_t acc, x, rho, am, g[8], mask;
    int i;
    int j;
    int bstr;
    

#ifdef __INTEL_COMPILER
    int bits = 64 - _lead_zcnt64(b);
#elif defined(__GNUC__)
    int bits = 64 - __builtin_clzll(b);
#elif defined _MSC_VER
    int bits = 64 - _lead_zcnt64(b);
#endif

    x = (((n + 2) & 4) << 1) + n;   // here x*a==1 mod 2**4
    x *= 2 - n * x;                 // here x*a==1 mod 2**8
    x *= 2 - n * x;                 // here x*a==1 mod 2**16
    x *= 2 - n * x;                 // here x*a==1 mod 2**32         
    x *= 2 - n * x;                 // here x*a==1 mod 2**64

    rho = (uint64_t)0 - x;
    am = u64div(2, n);              // put 2 into Monty rep.

    // precomputations, b^i for 0 <= i < 2^k
    g[1] = am;
    g[2] = mulredc(g[1], am, n, rho);
    g[3] = mulredc(g[2], am, n, rho);
    g[4] = sqrredc(g[2], n, rho);
    g[5] = mulredc(g[4], am, n, rho);
    g[6] = sqrredc(g[3], n, rho);
    g[7] = mulredc(g[6], am, n, rho);

    // L-R windowed exponentiation.
    mask = 0x7ULL << (bits - 3);

    // the first 5 iterations can be done with no mod if BITS==64,
    // before we put the accumulator into monty rep.
    // iter1: sqr 1 and mul
    acc = 1;
    bstr = (b & mask) >> (bits - 3);
    if (bstr > 0)
        acc = acc * (1ULL << bstr);
    mask >>= 3;
    bits -= 3;
    acc *= acc;
    acc *= acc;
    bstr = (b & mask) >> (bits - 2);
    if (bstr > 0)
        acc = acc * (1ULL << bstr);
    mask >>= 2;
    bits -= 2;

    acc = u64div(acc, n);             // put acc into Monty rep.
    //acc = u64div(1, n);

    while (bits > 3)
    {
        bstr = (b & mask) >> (bits - 3);
        acc = sqrredc(acc, n, rho);
        acc = sqrredc(acc, n, rho);
        acc = sqrredc(acc, n, rho);
        mask >>= 3;
        if (bstr > 0)
            acc = mulredc(acc, g[bstr], n, rho);
        bits -= 3;
    }

    bstr = (b & mask);
    for (j = 0; mask > 0; j++)
    {
        acc = sqrredc(acc, n, rho);
        mask >>= 1;
    }

    if (bstr > 0)
        acc = mulredc(acc, g[bstr], n, rho);

    // take result out of monty rep.
    acc = mulredc(acc, 1, n, rho);

    // final check to ensure c < N
    if (acc >= n)
    {
        acc -= n;
    }

    return acc;
}

/* the following recursion base-case is specialized for relations
   containing <= 3 rational and/or algebraic large primes. The rest
   of the batch factoring handles an arbitrary number of large primes,
   so only this routine needs to change when more advanced factoring
   code becomes available. */

void check_batch_relation(relation_batch_t *rb,
			uint32_t index, uint64_t * lcg_state,
			mpz_t prime_product) {

	uint32_t i;
    int j;
	cofactor_t *c = rb->relations + index;
	uint32_t *f = rb->factors + c->factor_list_word;
	uint32_t *lp1 = f + c->num_factors_r + c->num_factors_a;
	uint32_t *lp2 = lp1 + c->lp_r_num_words;
    mpz_ptr f1r = rb->f1r;
    mpz_ptr f1a = rb->f1a;
    mpz_ptr f2r = rb->f2r;
    mpz_ptr f2a = rb->f2a;
    mpz_ptr _small = rb->_small;     // rpcndr.h defines "small" as "char", problem for msvc
    mpz_ptr _large = rb->_large;
    mpz_ptr n = rb->n;
	uint32_t lp_r[MAX_LARGE_PRIMES];
	uint32_t lp_a[MAX_LARGE_PRIMES];
	uint32_t num_r, num_a;

	/* first compute gcd(prime_product, rational large cofactor).
	   The rational part will split into a cofactor with all 
	   factors <= the largest prime in rb->prime_product (stored
	   in f1r) and a cofactor with all factors larger (stored in f2r) */

    // the larger we make lp_cutoff_r, the more TLP's we will find, at
    // the cost of having to split more TLP candidates in f1r.

    c->success = 0;
    c->lp_r[0] = c->lp_r[1] = c->lp_r[2] = 1;

	if (c->lp_r_num_words) {

        mpz_set_ui(n, lp1[c->lp_r_num_words - 1]);
        for (j = c->lp_r_num_words - 2; j >= 0; j--)
        {
            mpz_mul_2exp(n, n, 32);
            mpz_add_ui(n, n, lp1[j]);
        }

        //gmp_printf("processing n = %Zd: ", n); // , prime_product has %u bits\n", n, mpz_sizeinbase(prime_product, 2));

		if (mpz_sizeinbase(n, 2) < 32) {
            // this input is already small, assign it to the
            // "all factors <= the largest prime" side, and the other side 
            // has no factors.  
            // todo: this actually might not be true (it could be larger and
            // therefore technically belong in f2r), but does it matter?
            mpz_set(f1r, n);
            mpz_set_ui(f2r, 1);
		}
		else {
            mpz_gcd(f1r, prime_product, n);

            if (mpz_cmp_ui(f1r, 1) == 0)
            {
                // if the gcd is 1, that means all factors
                // of n are larger than what made up our 
                // pre-multiplied product, so assign everything to f2r.
                mpz_set(f2r, n);
            }
            else
            {
                // assign the leftover portion to f2r.
                mpz_tdiv_q(f2r, n, f1r);
            }
		}

		/* give up on this relation if
		     - f1r has a single factor, and that factor
		       exceeds the rational large prime cutoff
		     - f1r is more than 3 words long */

        if ((mpz_sizeinbase(f1r, 2) > 96) || 
            ((mpz_sizeinbase(f1r, 2) < 32) && (mpz_get_ui(f1r) > rb->lp_cutoff_r)))
        {
            return;
        }

		/* give up on this relation if
		     - f2r has a single factor, and that factor
		       exceeds the rational large prime cutoff
		     - f2r has two words and exceeds the square of the
		       cutoff (meaning at least one of the two factors in
		       f2r would exceed the large prime cutoff)
		     - f2r has two words and is smaller than the square
		       of the largest prime in rb->prime_product, meaning
		       that f2r is prime and way too large
		     - f2r has 3+ words. We don't know anything about
		       f2r in this case, except that all factors exceed the
		       largest prime in rb->prime_product, and it would be
		       much too expensive to find out anything more */

        uint32_t r_cutoff_bits = spBits(rb->lp_cutoff_r) * 2;

        if ((mpz_sizeinbase(f2r, 2) > r_cutoff_bits) || ((mpz_get_ui(f2r) > 1) &&
            (((mpz_sizeinbase(f2r, 2) < 32) && (mpz_get_ui(f2r) > rb->lp_cutoff_r)) ||
            ((mpz_sizeinbase(f2r, 2) < 64) &&
                ((mpz_cmp(f2r, rb->max_prime2) <= 0) ||
                (mpz_cmp(f2r, rb->lp_cutoff_r2) > 0))))))
        {
            return;
        }

        //if ((mpz_sizeinbase(f2r, 2) >= 64) && ((mpz_sizeinbase(f2r, 2) < 96)))
        if (0)
        {
            // experiment with the last case where f2r possibly splits
            // into a TLP but we know that all of the factors are
            // larger than primes in our GCD.  It is expensive to learn more
            // but for huge jobs the expense may be worth it.
            // We expend a small amount of ECM and hope to get lucky.
            if (mpz_probab_prime_p(f2r, 1))
            {
                //printf("abort 3\n");
                return;
            }

            for (i = num_r = num_a = 0; i < MAX_LARGE_PRIMES; i++)
                lp_r[i] = lp_a[i] = 1;

            tinyecm(f2r, _small, 70, 70 * 25, 16, lcg_state, 0);

            if ((mpz_sizeinbase(_small, 2) > 32) || (mpz_get_ui(_small) > rb->lp_cutoff_r))
            {
                //printf("abort 8\n");
                return;
            }

            // small is acceptable, get large part.
            mpz_tdiv_q(_large, f2r, _small);
            lp_r[num_r++] = mpz_get_ui(_small);

            if (mpz_sizeinbase(_large, 2) > 64)
            {
                //printf("abort 8\n");
                return;
            }

            uint64_t e = mpz_get_ui(_large);
            if (pow2m(e - 1, e) == 1)
            {
                //printf("abort 3\n");
                return;
            }

            // large still on track, try to split it
            uint64_t f64 = do_uecm(mpz_get_ui(_large));

            // see if any factors found are acceptable
            if (f64 <= 1 || f64 > rb->lp_cutoff_r)
            {
                //printf("abort 9\n");
                return;
            }

            lp_r[num_r++] = f64;
            mpz_tdiv_q_ui(_large, _large, f64);

            if ((mpz_sizeinbase(_large, 2) > 32) || (mpz_get_ui(_large) > rb->lp_cutoff_r))
            {
                //printf("abort 10\n");
                return;
            }

            lp_r[num_r++] = mpz_get_ui(_large);

            /* yay! Another relation found */

            rb->num_success++;

            for (i = 0; i < MIN(3, num_r); i++)
            {
                c->lp_r[i] = lp_r[i];
            }
            c->success = num_r;
            return;
        }

	}
	else {
        mpz_set_ui(f1r, 0);
        mpz_set_ui(f2r, 0);
	}

	/* repeat with the algebraic unfactored part, if any */

	if (c->lp_a_num_words) {
        printf("should not be processing any algebraic side relation in QS!\n");
        mpz_set_ui(n, lp2[c->lp_a_num_words - 1]);
        for (j = c->lp_a_num_words - 2; j >= 0; j--)
        {
            mpz_mul_2exp(n, n, 32);
            mpz_add_ui(n, n, lp2[j]);
        }

        if (mpz_sizeinbase(n, 2) < 32) {
            mpz_set(f1a, n);
            mpz_set_ui(f2a, 1);
        }
        else {
            mpz_gcd(f1a, prime_product, n);

            if (mpz_cmp_ui(f1a, 1) == 1)
            {
                mpz_set(f2a, n);
            }
            else
            {
                mpz_tdiv_q(f2a, n, f1a);
            }

        }

        if ((mpz_sizeinbase(f1a, 2) >= 96) ||
            ((mpz_sizeinbase(f1a, 2) < 32) && (mpz_get_ui(f1a) > rb->lp_cutoff_a)))
        {
            return;
        }

        if ((mpz_sizeinbase(f2a, 2) >= 64) ||
            ((mpz_sizeinbase(f2a, 2) < 32) && (mpz_get_ui(f2a) > rb->lp_cutoff_a)) ||
            ((mpz_sizeinbase(f2a, 2) < 64) &&
            ((mpz_cmp(f2a, rb->max_prime2) <= 0) ||
                (mpz_cmp(f2a, rb->lp_cutoff_a2) > 0))))
        {
            return;
        }
	}
	else {
        mpz_set_ui(f1a, 0);
        mpz_set_ui(f2a, 0);
	}

	/* the relation isn't obviously bad; do more work
	   trying to factor everything. Note that when relations 
	   are expected to have three large primes then ~98% of 
	   relations do not make it to this point
	
	   Begin by performing compositeness tests on f2r and f2a,
	   which are necessary if they are two words in size and
	   f1r or f1a is not one (the latter being true means
	   that the sieving already performed the compositeness test) */

	for (i = num_r = num_a = 0; i < MAX_LARGE_PRIMES; i++)
		lp_r[i] = lp_a[i] = 1;

    if ((mpz_sizeinbase(f2r, 2) > 32) && (mpz_sizeinbase(f2r, 2) <= 64))
    {
        uint64_t e = mpz_get_ui(f2r);
        if (pow2m(e - 1, e) == 1)
        {
            return;
        }
    }
    if ((mpz_sizeinbase(f2a, 2) > 32) && (mpz_sizeinbase(f2a, 2) <= 64))
    {
        uint64_t e = mpz_get_ui(f2a);
        printf("should not be processing any algebraic side relation in QS!\n");
        if (pow2m(e - 1, e) == 1)
        {
            return;
        }
    }

	/* now perform all the factorizations that
	   require SQUFOF, since it is much faster than the
	   QS code. We have to check all of f[1|2][r|a] but
	   for relations with three large primes then at most
	   two of the four choices need factoring */

    if (mpz_sizeinbase(f1r, 2) <= 32)
    {
        if (mpz_get_ui(f1r) > 1)
            lp_r[num_r++] = mpz_get_ui(f1r);
    }
    else if (mpz_sizeinbase(f1r, 2) <= 64)
    {
        uint64_t f64 = do_uecm(mpz_get_ui(f1r));

        if (f64 <= 1 || f64 > rb->lp_cutoff_r)
        {
            return;
        }
        lp_r[num_r++] = f64;
        mpz_tdiv_q_ui(f1r, f1r, f64);

        if ((mpz_sizeinbase(f1r, 2) > 32) || (mpz_get_ui(f1r) > rb->lp_cutoff_r))
        {
            return;
        }

        lp_r[num_r++] = mpz_get_ui(f1r);
    }

    if (mpz_sizeinbase(f2r, 2) <= 32)
    {
        if (mpz_get_ui(f2r) > 1)
            lp_r[num_r++] = mpz_get_ui(f2r);
    }
    else if (mpz_sizeinbase(f2r, 2) <= 64)
    {
        uint64_t f64 = do_uecm(mpz_get_ui(f2r));

        if (f64 <= 1 || f64 > rb->lp_cutoff_r)
        {
            return;
        }

        lp_r[num_r++] = f64;
        mpz_tdiv_q_ui(f2r, f2r, f64);

        if ((mpz_sizeinbase(f2r, 2) > 32) || (mpz_get_ui(f2r) > rb->lp_cutoff_r))
        {
            return;
        }

        lp_r[num_r++] = mpz_get_ui(f2r);
    }

    if (mpz_sizeinbase(f1a, 2) <= 32)
    {
        if (mpz_get_ui(f1a) > 1)
            lp_a[num_a++] = mpz_get_ui(f1a);
    }
    else if (mpz_sizeinbase(f1a, 2) <= 64)
    {
        uint64_t f64 = do_uecm(mpz_get_ui(f1a));
        printf("should not be processing any algebraic side relation in QS!\n");
        if (f64 <= 1 || f64 > rb->lp_cutoff_a)
            return;
        lp_a[num_a++] = f64;
        mpz_tdiv_q_ui(f1a, f1a, f64);

        if ((mpz_sizeinbase(f1a, 2) > 32) || (mpz_get_ui(f1a) > rb->lp_cutoff_a))
            return;

        lp_a[num_a++] = mpz_get_ui(f1a);
    }

    if (mpz_sizeinbase(f2a, 2) <= 32)
    {
        if (mpz_get_ui(f2a) > 1)
            lp_a[num_a++] = mpz_get_ui(f2a);
    }
    else if (mpz_sizeinbase(f2a, 2) <= 64)
    {
        uint64_t f64 = do_uecm(mpz_get_ui(f2a));
        printf("should not be processing any algebraic side relation in QS!\n");

        if (f64 <= 1 || f64 > rb->lp_cutoff_a)
            return;
        lp_a[num_a++] = f64;
        mpz_tdiv_q_ui(f2a, f2a, f64);

        if ((mpz_sizeinbase(f2a, 2) > 32) || (mpz_get_ui(f2a) > rb->lp_cutoff_a))
            return;

        lp_a[num_a++] = mpz_get_ui(f2a);
    }

	/* only use expensive factoring methods when we know 
	   f1r and/or f1a splits into three large primes and 
	   we know all three primes are smaller than the 
	   largest prime in rb->prime_product. When the latter 
	   is a good deal smaller than the large prime cutoff
	   this happens extremely rarely */

	//if (f1r.nwords == 3) {
	//	if (tinyqs(&f1r, &t0, &t1) == 0)
	//		return;
    //
	//	small = &t0;
	//	_large = &t1;
	//	if (mp_cmp(small, _large) > 0) {
	//		small = &t1;
	//		_large = &t0;
	//	}
    //
	//	if (small->nwords > 1 || small->val[0] > rb->lp_cutoff_r)
	//		return;
	//	lp_r[num_r++] = small->val[0];
	//	i = squfof(_large);
	//	if (i <= 1 || i > rb->lp_cutoff_r)
	//		return;
	//	lp_r[num_r++] = i;
	//	mp_divrem_1(_large, i, _large);
	//	if (_large->nwords > 1 || _large->val[0] > rb->lp_cutoff_r)
	//		return;
	//	lp_r[num_r++] = _large->val[0];
	//}

    if (mpz_sizeinbase(f1r, 2) > 64) {
        //gmp_printf("attempting to factor %u-bit n = %Zx\n", mpz_sizeinbase(f1r, 2), f1r);

        //if (tinyqs(qs_params, f1r, small, _large) == 0)
        //    goto done;
        //
        //if (mpz_cmp(small, _large) > 0) {
        //    mpz_set(n, small);
        //    mpz_set(small, _large);
        //    mpz_set(_large, n);
        //}

        int B1, B2, curves, bits = mpz_sizeinbase(f1r, 2);
        int targetBits = bits / 3 + 1;
        if (targetBits <= 25)
        {
            B1 = 70;
            curves = 16;
        }
        else if (targetBits <= 26)
        {
            B1 = 85;
            curves = 16;
        }
        else if (targetBits <= 29)
        {
            B1 = 125;
            curves = 16;
        }
        else if (targetBits <= 31)
        {
            B1 = 165;
            curves = 24;
        }
        else if (targetBits <= 32)
        {
            B1 = 205;
            curves = 24;
        }
        else
        {
            printf("something's wrong, bits = %u, targetBits = %u\n", bits, targetBits);
        }

        tinyecm(f1r, _small, B1, B1 * 25, curves, lcg_state, 0);

        if ((mpz_sizeinbase(_small, 2) > 32) || (mpz_get_ui(_small) > rb->lp_cutoff_r))
        {
            return;
        }

        mpz_tdiv_q(_large, f1r, _small);
        lp_r[num_r++] = mpz_get_ui(_small);

        uint64_t f64 = do_uecm(mpz_get_ui(_large));

        if (f64 <= 1 || f64 > rb->lp_cutoff_r)
        {
            return;
        }

        lp_r[num_r++] = f64;
        mpz_tdiv_q_ui(_large, _large, f64);

        if ((mpz_sizeinbase(_large, 2) > 32) || (mpz_get_ui(_large) > rb->lp_cutoff_r))
        {
            return;
        }

        lp_r[num_r++] = mpz_get_ui(_large);
    }


	//if (f1a.nwords == 3) {
	//	if (tinyqs(&f1a, &t0, &t1) == 0)
	//		return;
    //
	//	_small = &t0;
	//	_large = &t1;
	//	if (mp_cmp(_small, _large) > 0) {
	//		_small = &t1;
	//		_large = &t0;
	//	}
    //
	//	if (_small->nwords > 1 || _small->val[0] > rb->lp_cutoff_a)
	//		return;
	//	lp_a[num_a++] = _small->val[0];
	//	i = squfof(_large);
	//	if (i <= 1 || i > rb->lp_cutoff_a)
	//		return;
	//	lp_a[num_a++] = i;
	//	mp_divrem_1(_large, i, _large);
	//	if (_large->nwords > 1 || _large->val[0] > rb->lp_cutoff_a)
	//		return;
	//	lp_a[num_a++] = _large->val[0];
	//}

    if (mpz_sizeinbase(f1a, 2) > 64) {
        //if (tinyqs(qs_params, f1a, _small, _large) == 0)
        //    return;
        //
        //if (mpz_cmp(_small, _large) > 0) {
        //    mpz_set(n, _small);
        //    mpz_set(_small, _large);
        //    mpz_set(_large, n);
        //}

        int B1, B2, curves, bits = mpz_sizeinbase(f1a, 2);
        int targetBits = bits / 3 + 1;

        printf("should not be processing any algebraic side relation in QS!\n");
        if (targetBits <= 25)
        {
            B1 = 70;
            curves = 16;
        }
        else if (targetBits <= 26)
        {
            B1 = 85;
            curves = 16;
        }
        else if (targetBits <= 29)
        {
            B1 = 125;
            curves = 16;
        }
        else if (targetBits <= 31)
        {
            B1 = 165;
            curves = 24;
        }
        else if (targetBits <= 32)
        {
            B1 = 205;
            curves = 24;
        }
        else
        {
            printf("something's wrong, bits = %u, targetBits = %u\n", bits, targetBits);
        }
        
        tinyecm(f1a, _small, B1, B1 * 25, curves, lcg_state, 0);

        if ((mpz_sizeinbase(_small, 2) > 32) || (mpz_get_ui(_small) > rb->lp_cutoff_a))
            return;

        mpz_tdiv_q(_large, f1a, _small);
        lp_a[num_a++] = mpz_get_ui(_small);

        uint64_t f64 = do_uecm(mpz_get_ui(_large));

        if (f64 <= 1 || f64 > rb->lp_cutoff_a)
            return;

        lp_a[num_a++] = f64;
        mpz_tdiv_q_ui(_large, _large, f64);

        if ((mpz_sizeinbase(_large, 2) > 32) || (mpz_get_ui(_large) > rb->lp_cutoff_a))
            return;

        lp_a[num_a++] = mpz_get_ui(_large);
    }

	/* yay! Another relation found */

	rb->num_success++;
    // it is complicated to save into the siqs relation buffer
    // from here, so we just mark this cofactor as a success and
    // calling code will put it into the buffer.
    //printf("batch index %d success %u, %u factors: %u,%u,%u\n", 
    //    index, rb->num_success, num_r, lp_r[0], lp_r[1], lp_r[2]);

    for (i = 0; i < MIN(3, num_r); i++)
    {
        c->lp_r[i] = lp_r[i];
    }
    c->success = num_r;

	//rb->print_relation(rb->savefile, c->a, c->b,
	//		f, c->num_factors_r, lp_r,
	//		f + c->num_factors_r, c->num_factors_a, lp_a);
    //buffer_relation(c->offset, lp_r, c->num_factors_r, f, 
    //    c->a, c->b, c->parity, NULL, NULL, 0, 1);

    //printf("success\n");
    return;
}

/*------------------------------------------------------------------*/
void compute_remainder_tree(bintree_t* tree, uint32_t first, uint32_t last,
				relation_batch_t *rb, uint64_t *lcg_state,
				mpz_t numerator) {

	/* recursively compute numerator % (each relation in 
	   rb->relations) */

	uint32_t mid = (first + last) / 2;
	mpz_t relation_prod, remainder;

	/* recursion base case: numerator already fits in
	   an mp_t, so manually compute the remainder and
	   postprocess each relation */

	if (mpz_sizeinbase(numerator, 2) <= (MAX_MP_WORDS * 32)) {
		if (mpz_sgn(numerator) > 0) {
			while (first <= last)
                check_batch_relation(rb, first++, numerator, lcg_state);
		}
		return;
	}

	/* multiply together the unfactored parts of all the
	   relations from first to last */

	mpz_init(relation_prod);
	mpz_init(remainder);
	multiply_relations(tree, first, last, rb, relation_prod);

	/* use the remainder to deal with the left and right
	   halves of the relation list */

	if (mpz_cmp(numerator, relation_prod) < 0) {
		mpz_clear(relation_prod);
		compute_remainder_tree(tree, first, mid, rb, lcg_state, numerator);
		compute_remainder_tree(tree, mid + 1, last, rb, lcg_state, numerator);
	}
	else {
		mpz_tdiv_r(remainder, numerator, relation_prod);
		mpz_clear(relation_prod);
		compute_remainder_tree(tree, first, mid, rb, lcg_state, remainder);
		compute_remainder_tree(tree, mid + 1, last, rb, lcg_state, remainder);
	}
	mpz_clear(remainder);
}

/*------------------------------------------------------------------*/
void relation_batch_init(FILE *logfile, relation_batch_t *rb,
			uint32_t min_prime, uint32_t max_prime,
			uint32_t lp_cutoff_r, uint32_t lp_cutoff_a, 
			print_relation_t print_relation,
            int do_prime_product) {

    // qs_savefile_t *savefile,

    mpz_init(rb->prime_product);

    if (do_prime_product)
    {
        prime_sieve_t sieve;
        uint32_t num_primes, p;

        /* count the number of primes to multiply. Knowing this
           in advance makes the recursion a lot easier, at the cost
           of a small penalty in runtime */

        init_prime_sieve(&sieve, min_prime + 1, max_prime);
        p = min_prime;
        num_primes = 0;
        while (p < max_prime) {
        	p = get_next_prime(&sieve);
        	num_primes++;
        }
        free_prime_sieve(&sieve);

        ///* compute the product of primes */

        logprint(logfile, "multiplying %u primes from %u to %u\n",
        		num_primes, min_prime, max_prime);

        init_prime_sieve(&sieve, min_prime, max_prime);
        multiply_primes(0, num_primes - 2, &sieve, rb->prime_product);
        free_prime_sieve(&sieve);

        logprint(logfile, "multiply complete, product has %u bits\n",
        			(uint32_t)mpz_sizeinbase(rb->prime_product, 2));
    }
					
	//rb->savefile = savefile;
	rb->print_relation = print_relation;
    rb->conversion_ratio = 0.0;

	/* compute the cutoffs used by the recursion base-case. Large
	   primes have a maximum size specified as input arguments, 
	   but numbers that can be passed to the SQUFOF routine are
	   limited to size 2^62 */

	rb->lp_cutoff_r = lp_cutoff_r;
    // no longer an issue with new ecm routines.
	//lp_cutoff_r = MIN(lp_cutoff_r, 0x7fffffff);
	mpz_init(rb->lp_cutoff_r2);
    mpz_set_ui(rb->lp_cutoff_r2, lp_cutoff_r);
    mpz_mul_ui(rb->lp_cutoff_r2, rb->lp_cutoff_r2, lp_cutoff_r);

	rb->lp_cutoff_a = lp_cutoff_a;
    // no longer an issue with new ecm routines.
    //lp_cutoff_a = MIN(lp_cutoff_a, 0x7fffffff);
    mpz_init(rb->lp_cutoff_a2);
    mpz_set_ui(rb->lp_cutoff_a2, lp_cutoff_a);
    mpz_mul_ui(rb->lp_cutoff_a2, rb->lp_cutoff_a2, lp_cutoff_a);

    mpz_init(rb->max_prime2);
    mpz_set_ui(rb->max_prime2, max_prime);
    mpz_mul_ui(rb->max_prime2, rb->max_prime2, max_prime);

	/* allocate lists for relations and their factors */

	rb->target_relations = 500000;
	rb->num_relations = 0;
	rb->num_relations_alloc = 1000;
	rb->relations = (cofactor_t *)xmalloc(rb->num_relations_alloc *
						sizeof(cofactor_t));

	rb->num_factors = 0;
	rb->num_factors_alloc = 10000;
	rb->factors = (uint32_t *)xmalloc(rb->num_factors_alloc *
						sizeof(uint32_t));

    mpz_init(rb->n);
    mpz_init(rb->f1r);
    mpz_init(rb->f2r);
    mpz_init(rb->f1a);
    mpz_init(rb->f2a);
    mpz_init(rb->t0);
    mpz_init(rb->t1);
    mpz_init(rb->_small);
    mpz_init(rb->_large);
}

/*------------------------------------------------------------------*/
void relation_batch_free(relation_batch_t *rb) {

	mpz_clear(rb->prime_product);
    mpz_clear(rb->max_prime2);
    mpz_clear(rb->lp_cutoff_a2);
    mpz_clear(rb->lp_cutoff_r2);
	free(rb->relations);
	free(rb->factors);

    mpz_clear(rb->n);
    mpz_clear(rb->f1r);
    mpz_clear(rb->f2r);
    mpz_clear(rb->f1a);
    mpz_clear(rb->f2a);
    mpz_clear(rb->t0);
    mpz_clear(rb->t1);
    mpz_clear(rb->_small);
    mpz_clear(rb->_large);
}

/*------------------------------------------------------------------*/
void relation_batch_add(uint32_t a, uint32_t b, int32_t offset,
			uint32_t *factors_r, uint32_t num_factors_r, 
			mpz_t unfactored_r_in,
			uint32_t *factors_a, uint32_t num_factors_a, 
			mpz_t unfactored_a_in,
            mpz_t tmp_in,
			relation_batch_t *rb) {

	uint32_t i;
	uint32_t *f;
	cofactor_t *c;

	/* add one relation to the batch */
	if (rb->num_relations == rb->num_relations_alloc) {
		rb->num_relations_alloc *= 2;
		rb->relations = (cofactor_t *)xrealloc(rb->relations,
					rb->num_relations_alloc *
					sizeof(cofactor_t));
	}
	c = rb->relations + rb->num_relations++;
	c->a = a;
	c->b = b;
    c->signed_offset = offset;
	c->num_factors_r = num_factors_r;
	c->num_factors_a = num_factors_a;
	c->lp_r_num_words = 0;
	c->lp_a_num_words = 0;
	c->factor_list_word = rb->num_factors;
    c->lp_r[0] = c->lp_r[1] = c->lp_r[2] = 0;
    c->success = 0;

	/* add its small factors */

	if (rb->num_factors + num_factors_r + num_factors_a +
			2 * MAX_LARGE_PRIMES >= rb->num_factors_alloc) {
		rb->num_factors_alloc *= 2;
		rb->factors = (uint32_t *)xrealloc(rb->factors,
					rb->num_factors_alloc *
					sizeof(uint32_t));
	}
	f = rb->factors + rb->num_factors;

	for (i = 0; i < num_factors_r; i++)
		f[i] = factors_r[i];
	f += i;
	rb->num_factors += i;

	for (i = 0; i < num_factors_a; i++)
		f[i] = factors_a[i];
	f += i;
	rb->num_factors += i;

	/* add its large factors (if any) */

	if (mpz_cmp_ui(unfactored_r_in, 1) > 0) {

        mpz_set(rb->t0, unfactored_r_in);
        i = 0;
        while (mpz_cmp_ui(rb->t0, 0) > 0)
        {
            f[i++] = (uint32_t)mpz_get_ui(rb->t0);
            mpz_tdiv_q_2exp(rb->t0, rb->t0, 32);
        }

		f += i;
		rb->num_factors += i;
		c->lp_r_num_words = i;
	}

    if (mpz_cmp_ui(unfactored_a_in, 1) > 0) {
        mpz_set(rb->t0, unfactored_a_in);
        i = 0;
        while (mpz_cmp_ui(rb->t0, 0) > 0)
        {
            f[i++] = (uint32_t)mpz_get_ui(rb->t0);
            mpz_tdiv_q_2exp(rb->t0, rb->t0, 32);
        }
    
        f += i;
        rb->num_factors += i;
        c->lp_a_num_words = i;
    }
}
	
/*------------------------------------------------------------------*/



uint32_t relation_batch_run(relation_batch_t *rb, uint64_t *lcg_state) {
    // recursive batch GCD, with a tree storage enhancement
    // to avoid re-calculating many of the products, at a cost
    // of additional RAM.
    // with TreeCutoff = 16 and ~500k relations in a batch, the
    // tree occupies about 64MB per thread and yields a speedup
    // of about 2x for the whole batch processing.
	rb->num_success = 0;
	if (rb->num_relations > 0) {

        bintree_t tree;
        int i;

        tree.nodes = (bintree_element_t*)xmalloc(16 * sizeof(bintree_element_t));
        tree.alloc = 16;
        tree.size = 1;

        // root node holds the product of all relations and initially
        // is empty.
        mpz_init(tree.nodes[0].prod);
        tree.nodes[0].low = 0;
        tree.nodes[0].high = rb->num_relations - 1;
        tree.nodes[0].left_id = -1;
        tree.nodes[0].right_id = -1;
        tree.nodes[0].complete = 0;
        //printf("\n\n");
        //printf("created tree root at (%u:%u)\n", 0, rb->num_relations - 1);

        // this already traverses the tree... just build in
        // the capability to add and reuse nodes and we 
        // should be there.
		compute_remainder_tree(&tree, 0, rb->num_relations - 1,
					rb, lcg_state, rb->prime_product);

        uint32_t bytes = 0;
        for (i = 0; i < tree.size; i++)
        {
            bytes += mpz_sizeinbase(tree.nodes[i].prod, 2) / 8;
            mpz_clear(tree.nodes[i].prod);
        }

        //printf("cleaned up %u nodes occupying ~ %u bytes\n", 
        //    tree.size, bytes + tree.size * sizeof(bintree_element_t));
        
        free(tree.nodes);
	}

	/* wipe out batched relations */

	//rb->num_relations = 0;
	//rb->num_factors = 0;
	return rb->num_success;
}


