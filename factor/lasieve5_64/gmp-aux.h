/*1:*/
#line 12 "gmp-aux.w"

void adjust_mpz_bufsize(mpz_t**x,size_t*alloc_ptr,size_t size,size_t increment);
int string2mpz(mpz_t rop,char*x,int base);
#if __GNU_MP_VERSION < 3 || (__GNU_MP_VERSION == 3 && __GNU_MP_VERSION_MINOR == 0)




#endif
#line 21 "gmp-aux.w"

#define mpz_add_si(r,o,s) \
  ({ int _s; _s= (s); \
     _s>=0 ? mpz_add_ui(r,o,(mp_limb_t)(_s)) : mpz_sub_ui(r,o,(mp_limb_t)(-_s)); })

void mpz_set_ull(mpz_t targ,ullong src);
ullong mpz_get_ull(mpz_t src);
int mpz_cmp_ull(mpz_t op1,ullong op2);
#ifdef ULL_NO_UL
void mpz_ull_init();
void mpz_mul_ull(mpz_t rop,mpz_t op1,ullong op2);
int mpz_fits_sllong_p(mpz_t);
int mpz_fits_ullong_p(mpz_t);
long long int mpz_get_sll(mpz_t);
void mpz_set_sll(mpz_t,long long int);
unsigned long long mpz_get_ull(mpz_t);
void mpz_set_ull(mpz_t,unsigned long long);
void mpz_add_ull(mpz_t rop, mpz_t op1, ullong op2);
void mpz_tdiv_q_ull(mpz_t rop, mpz_t op1, ullong op2);
#else
#line 39 "gmp-aux.w"
#define mpz_ull_init()
#define mpz_mul_ull mpz_mul_ui
#define mpz_mul_sll mpz_mul_si
#define mpz_fits_sllong_p mpz_fits_slong_p
#define mpz_fits_ullong_p mpz_fits_ulong_p
#define mpz_get_sll mpz_get_si
#define mpz_set_sll mpz_set_si
#define mpz_get_ull mpz_get_ui
#define mpz_set_ull mpz_set_ui
#define mpz_add_ull mpz_add_ui
#define mpz_tdiv_q_ull mpz_tdiv_q_ui
#endif
#line 49 "gmp-aux.w"
/*:1*/
