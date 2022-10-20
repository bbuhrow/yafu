
void td_pprod_get(mpz_t rop, u64_t p0, u64_t p1);
void td_pprod_store(mpz_t op, u64_t p0, u64_t p1);
int td_pprod_load(mpz_t rop, u64_t p0, u64_t p1);
void td_init(mpz_t *op, u32_t nf);
void td_execute(mpz_t pprod);
void td_end(mpz_t *rop, mpz_t *op);

