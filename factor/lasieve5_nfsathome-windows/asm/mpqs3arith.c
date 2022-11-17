#ifdef MOD3
static u32_t mpqs3_mod3(u32_t *a, u32_t n)
{
  u64_t ha;

  ha=(u64_t)(a[2]); ha%=(u64_t)n; ha<<=32;
  ha+=(u64_t)(a[1]); ha%=(u64_t)n; ha<<=32;
  ha+=(u64_t)(a[0]); ha%=(u64_t)n;
  return (u32_t)ha;
}
#else
#define mpqs3_mod3 asm_mpqs3_mod3
#endif
