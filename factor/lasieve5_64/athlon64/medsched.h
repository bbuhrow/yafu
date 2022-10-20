/* SMJS Old style prototypes replaced
u32_t *medsched(u32_t*,u32_t*,u32_t*,u32_t**,u32_t,u32_t);
u32_t *medsched_1(u32_t*,u32_t*,u32_t*,u32_t,unsigned char *,
                  unsigned char);
*/
u32_t *
medsched(u32_t *ri, u32_t *ij_ptr, u32_t *ij_ptr_ub,
         u32_t **sched_ptr,u32_t fbi_offs, u32_t ot);
u32_t *
medsched_1(u32_t *ri, u32_t *ij_ptr, u32_t *ij_ptr_ub, u32_t ot,
           unsigned char *si, unsigned char lo);
