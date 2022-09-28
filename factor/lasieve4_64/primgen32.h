/*2:*/
#line 51 "primgen32.w"

typedef struct{
u32_t Pind;
u32_t first_in_sieve;
u32_t nPrim;
u32_t Prime;
unsigned char*PDiffs;
u32_t use_private;
u32_t PDiffs_allocated;
}pr32_struct;

void initprime32(pr32_struct*ps);
u32_t firstprime32(pr32_struct*ps);
u32_t nextprime32(pr32_struct*ps);
void clearprime32(pr32_struct*ps);
void p32_seek(pr32_struct*ps,u32_t lb);

/*:2*/
