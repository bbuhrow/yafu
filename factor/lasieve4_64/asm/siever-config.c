/*3:*/
#line 70 "siever-config.w"

#include <limits.h> 
#include "siever-config.h"
static void
siever_init(void)
{
}

/*:3*//*5:*/
#line 83 "siever-config.w"

#ifndef MPQS_ONLY
const unsigned long schedule_primebounds[N_PRIMEBOUNDS]= 
{0x100000,0x200000,0x400000,0x800000,0x1000000,0x2000000,0x4000000,
0x8000000,ULONG_MAX};

const unsigned long schedule_sizebits[N_PRIMEBOUNDS]= {20,21,22,23,24,25,26,27,32};
#endif/*:5*/
