
static u16_t *
mmx_xmalloc(size_t n);
void
MMX_TdAllocate(int jps_arg,size_t s0,size_t s1);
u16_t*
MMX_TdInit(int side,u16_t *x,u16_t *x_ub,u32_t *pbound_ptr,
	   int initialize);
       void
MMX_TdUpdate(int side,int j_step);
u32_t *
MMX_Td(u32_t *pbuf,int side,u16_t strip_i);
