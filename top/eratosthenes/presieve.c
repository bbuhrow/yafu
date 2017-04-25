#include "soe.h"
#include <immintrin.h> //<immintrin.h>

//masks for sieving multiple locations at once: small primes
//that hit a 64 bit interval more than once.  
int _64_MOD_P[9] = { 4, 1, 9, 12, 13, 7, 18, 6, 2 };

// 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89
int _256_MOD_P[22] = { 1, 4, 3, 9, 1, 9, 3, 24, 8, 34, 10, 41, 21, 44, 20, 12, 55, 43, 37, 19, 7, 78};

uint64 _5_MASKS_256[5 * 4] = {
    0xef7bdef7bdef7bdeULL, 0xdef7bdef7bdef7bdULL, 0xbdef7bdef7bdef7bULL, 0x7bdef7bdef7bdef7ULL,
    0xdef7bdef7bdef7bdULL, 0xbdef7bdef7bdef7bULL, 0x7bdef7bdef7bdef7ULL, 0xf7bdef7bdef7bdefULL,
    0xbdef7bdef7bdef7bULL, 0x7bdef7bdef7bdef7ULL, 0xf7bdef7bdef7bdefULL, 0xef7bdef7bdef7bdeULL,
    0x7bdef7bdef7bdef7ULL, 0xf7bdef7bdef7bdefULL, 0xef7bdef7bdef7bdeULL, 0xdef7bdef7bdef7bdULL,
    0xf7bdef7bdef7bdefULL, 0xef7bdef7bdef7bdeULL, 0xdef7bdef7bdef7bdULL, 0xbdef7bdef7bdef7bULL };

uint64 _7_MASKS_256[7 * 4] = {
    0x7efdfbf7efdfbf7eULL, 0xbf7efdfbf7efdfbfULL, 0xdfbf7efdfbf7efdfULL, 0xefdfbf7efdfbf7efULL,
    0xfdfbf7efdfbf7efdULL, 0x7efdfbf7efdfbf7eULL, 0xbf7efdfbf7efdfbfULL, 0xdfbf7efdfbf7efdfULL,
    0xfbf7efdfbf7efdfbULL, 0xfdfbf7efdfbf7efdULL, 0x7efdfbf7efdfbf7eULL, 0xbf7efdfbf7efdfbfULL,
    0xf7efdfbf7efdfbf7ULL, 0xfbf7efdfbf7efdfbULL, 0xfdfbf7efdfbf7efdULL, 0x7efdfbf7efdfbf7eULL,
    0xefdfbf7efdfbf7efULL, 0xf7efdfbf7efdfbf7ULL, 0xfbf7efdfbf7efdfbULL, 0xfdfbf7efdfbf7efdULL,
    0xdfbf7efdfbf7efdfULL, 0xefdfbf7efdfbf7efULL, 0xf7efdfbf7efdfbf7ULL, 0xfbf7efdfbf7efdfbULL,
    0xbf7efdfbf7efdfbfULL, 0xdfbf7efdfbf7efdfULL, 0xefdfbf7efdfbf7efULL, 0xf7efdfbf7efdfbf7ULL };

uint64 _11_MASKS_256[11 * 4] = {
    0xff7feffdffbff7feULL, 0xfdffbff7feffdffbULL, 0xf7feffdffbff7fefULL, 0xdffbff7feffdffbfULL,
    0xfeffdffbff7feffdULL, 0xfbff7feffdffbff7ULL, 0xeffdffbff7feffdfULL, 0xbff7feffdffbff7fULL,
    0xfdffbff7feffdffbULL, 0xf7feffdffbff7fefULL, 0xdffbff7feffdffbfULL, 0x7feffdffbff7feffULL,
    0xfbff7feffdffbff7ULL, 0xeffdffbff7feffdfULL, 0xbff7feffdffbff7fULL, 0xffdffbff7feffdffULL,
    0xf7feffdffbff7fefULL, 0xdffbff7feffdffbfULL, 0x7feffdffbff7feffULL, 0xffbff7feffdffbffULL,
    0xeffdffbff7feffdfULL, 0xbff7feffdffbff7fULL, 0xffdffbff7feffdffULL, 0xff7feffdffbff7feULL,
    0xdffbff7feffdffbfULL, 0x7feffdffbff7feffULL, 0xffbff7feffdffbffULL, 0xfeffdffbff7feffdULL,
    0xbff7feffdffbff7fULL, 0xffdffbff7feffdffULL, 0xff7feffdffbff7feULL, 0xfdffbff7feffdffbULL,
    0x7feffdffbff7feffULL, 0xffbff7feffdffbffULL, 0xfeffdffbff7feffdULL, 0xfbff7feffdffbff7ULL,
    0xffdffbff7feffdffULL, 0xff7feffdffbff7feULL, 0xfdffbff7feffdffbULL, 0xf7feffdffbff7fefULL,
    0xffbff7feffdffbffULL, 0xfeffdffbff7feffdULL, 0xfbff7feffdffbff7ULL, 0xeffdffbff7feffdfULL };

uint64 _13_MASKS_256[13 * 4] = {
    0xffefff7ffbffdffeULL, 0xffdffefff7ffbffdULL, 0xffbffdffefff7ffbULL, 0xff7ffbffdffefff7ULL,
    0xffdffefff7ffbffdULL, 0xffbffdffefff7ffbULL, 0xff7ffbffdffefff7ULL, 0xfefff7ffbffdffefULL,
    0xffbffdffefff7ffbULL, 0xff7ffbffdffefff7ULL, 0xfefff7ffbffdffefULL, 0xfdffefff7ffbffdfULL,
    0xff7ffbffdffefff7ULL, 0xfefff7ffbffdffefULL, 0xfdffefff7ffbffdfULL, 0xfbffdffefff7ffbfULL,
    0xfefff7ffbffdffefULL, 0xfdffefff7ffbffdfULL, 0xfbffdffefff7ffbfULL, 0xf7ffbffdffefff7fULL,
    0xfdffefff7ffbffdfULL, 0xfbffdffefff7ffbfULL, 0xf7ffbffdffefff7fULL, 0xefff7ffbffdffeffULL,
    0xfbffdffefff7ffbfULL, 0xf7ffbffdffefff7fULL, 0xefff7ffbffdffeffULL, 0xdffefff7ffbffdffULL,
    0xf7ffbffdffefff7fULL, 0xefff7ffbffdffeffULL, 0xdffefff7ffbffdffULL, 0xbffdffefff7ffbffULL,
    0xefff7ffbffdffeffULL, 0xdffefff7ffbffdffULL, 0xbffdffefff7ffbffULL, 0x7ffbffdffefff7ffULL,
    0xdffefff7ffbffdffULL, 0xbffdffefff7ffbffULL, 0x7ffbffdffefff7ffULL, 0xfff7ffbffdffefffULL,
    0xbffdffefff7ffbffULL, 0x7ffbffdffefff7ffULL, 0xfff7ffbffdffefffULL, 0xffefff7ffbffdffeULL,
    0x7ffbffdffefff7ffULL, 0xfff7ffbffdffefffULL, 0xffefff7ffbffdffeULL, 0xffdffefff7ffbffdULL,
    0xfff7ffbffdffefffULL, 0xffefff7ffbffdffeULL, 0xffdffefff7ffbffdULL, 0xffbffdffefff7ffbULL };

uint64 _17_MASKS_256[17 * 4] = {
    0xfff7fffbfffdfffeULL, 0xff7fffbfffdfffefULL, 0xf7fffbfffdfffeffULL, 0x7fffbfffdfffefffULL,
    0xffeffff7fffbfffdULL, 0xfeffff7fffbfffdfULL, 0xeffff7fffbfffdffULL, 0xffff7fffbfffdfffULL,
    0xffdfffeffff7fffbULL, 0xfdfffeffff7fffbfULL, 0xdfffeffff7fffbffULL, 0xfffeffff7fffbfffULL,
    0xffbfffdfffeffff7ULL, 0xfbfffdfffeffff7fULL, 0xbfffdfffeffff7ffULL, 0xfffdfffeffff7fffULL,
    0xff7fffbfffdfffefULL, 0xf7fffbfffdfffeffULL, 0x7fffbfffdfffefffULL, 0xfffbfffdfffeffffULL,
    0xfeffff7fffbfffdfULL, 0xeffff7fffbfffdffULL, 0xffff7fffbfffdfffULL, 0xfff7fffbfffdfffeULL,
    0xfdfffeffff7fffbfULL, 0xdfffeffff7fffbffULL, 0xfffeffff7fffbfffULL, 0xffeffff7fffbfffdULL,
    0xfbfffdfffeffff7fULL, 0xbfffdfffeffff7ffULL, 0xfffdfffeffff7fffULL, 0xffdfffeffff7fffbULL,
    0xf7fffbfffdfffeffULL, 0x7fffbfffdfffefffULL, 0xfffbfffdfffeffffULL, 0xffbfffdfffeffff7ULL,
    0xeffff7fffbfffdffULL, 0xffff7fffbfffdfffULL, 0xfff7fffbfffdfffeULL, 0xff7fffbfffdfffefULL,
    0xdfffeffff7fffbffULL, 0xfffeffff7fffbfffULL, 0xffeffff7fffbfffdULL, 0xfeffff7fffbfffdfULL,
    0xbfffdfffeffff7ffULL, 0xfffdfffeffff7fffULL, 0xffdfffeffff7fffbULL, 0xfdfffeffff7fffbfULL,
    0x7fffbfffdfffefffULL, 0xfffbfffdfffeffffULL, 0xffbfffdfffeffff7ULL, 0xfbfffdfffeffff7fULL,
    0xffff7fffbfffdfffULL, 0xfff7fffbfffdfffeULL, 0xff7fffbfffdfffefULL, 0xf7fffbfffdfffeffULL,
    0xfffeffff7fffbfffULL, 0xffeffff7fffbfffdULL, 0xfeffff7fffbfffdfULL, 0xeffff7fffbfffdffULL,
    0xfffdfffeffff7fffULL, 0xffdfffeffff7fffbULL, 0xfdfffeffff7fffbfULL, 0xdfffeffff7fffbffULL,
    0xfffbfffdfffeffffULL, 0xffbfffdfffeffff7ULL, 0xfbfffdfffeffff7fULL, 0xbfffdfffeffff7ffULL };

uint64 _19_MASKS_256[19 * 4] = {
    0xfdffffbffff7fffeULL, 0xfffbffff7fffefffULL, 0xbffff7fffeffffdfULL, 0xff7fffeffffdffffULL,
    0xfbffff7fffeffffdULL, 0xfff7fffeffffdfffULL, 0x7fffeffffdffffbfULL, 0xfeffffdffffbffffULL,
    0xf7fffeffffdffffbULL, 0xffeffffdffffbfffULL, 0xffffdffffbffff7fULL, 0xfdffffbffff7fffeULL,
    0xeffffdffffbffff7ULL, 0xffdffffbffff7fffULL, 0xffffbffff7fffeffULL, 0xfbffff7fffeffffdULL,
    0xdffffbffff7fffefULL, 0xffbffff7fffeffffULL, 0xffff7fffeffffdffULL, 0xf7fffeffffdffffbULL,
    0xbffff7fffeffffdfULL, 0xff7fffeffffdffffULL, 0xfffeffffdffffbffULL, 0xeffffdffffbffff7ULL,
    0x7fffeffffdffffbfULL, 0xfeffffdffffbffffULL, 0xfffdffffbffff7ffULL, 0xdffffbffff7fffefULL,
    0xffffdffffbffff7fULL, 0xfdffffbffff7fffeULL, 0xfffbffff7fffefffULL, 0xbffff7fffeffffdfULL,
    0xffffbffff7fffeffULL, 0xfbffff7fffeffffdULL, 0xfff7fffeffffdfffULL, 0x7fffeffffdffffbfULL,
    0xffff7fffeffffdffULL, 0xf7fffeffffdffffbULL, 0xffeffffdffffbfffULL, 0xffffdffffbffff7fULL,
    0xfffeffffdffffbffULL, 0xeffffdffffbffff7ULL, 0xffdffffbffff7fffULL, 0xffffbffff7fffeffULL,
    0xfffdffffbffff7ffULL, 0xdffffbffff7fffefULL, 0xffbffff7fffeffffULL, 0xffff7fffeffffdffULL,
    0xfffbffff7fffefffULL, 0xbffff7fffeffffdfULL, 0xff7fffeffffdffffULL, 0xfffeffffdffffbffULL,
    0xfff7fffeffffdfffULL, 0x7fffeffffdffffbfULL, 0xfeffffdffffbffffULL, 0xfffdffffbffff7ffULL,
    0xffeffffdffffbfffULL, 0xffffdffffbffff7fULL, 0xfdffffbffff7fffeULL, 0xfffbffff7fffefffULL,
    0xffdffffbffff7fffULL, 0xffffbffff7fffeffULL, 0xfbffff7fffeffffdULL, 0xfff7fffeffffdfffULL,
    0xffbffff7fffeffffULL, 0xffff7fffeffffdffULL, 0xf7fffeffffdffffbULL, 0xffeffffdffffbfffULL,
    0xff7fffeffffdffffULL, 0xfffeffffdffffbffULL, 0xeffffdffffbffff7ULL, 0xffdffffbffff7fffULL,
    0xfeffffdffffbffffULL, 0xfffdffffbffff7ffULL, 0xdffffbffff7fffefULL, 0xffbffff7fffeffffULL };

uint64 _23_MASKS_256[23 * 4] = {
    0xffffbfffff7ffffeULL, 0xfff7ffffefffffdfULL, 0xfefffffdfffffbffULL, 0xdfffffbfffff7fffULL,
    0xffff7ffffefffffdULL, 0xffefffffdfffffbfULL, 0xfdfffffbfffff7ffULL, 0xbfffff7ffffeffffULL,
    0xfffefffffdfffffbULL, 0xffdfffffbfffff7fULL, 0xfbfffff7ffffefffULL, 0x7ffffefffffdffffULL,
    0xfffdfffffbfffff7ULL, 0xffbfffff7ffffeffULL, 0xf7ffffefffffdfffULL, 0xfffffdfffffbffffULL,
    0xfffbfffff7ffffefULL, 0xff7ffffefffffdffULL, 0xefffffdfffffbfffULL, 0xfffffbfffff7ffffULL,
    0xfff7ffffefffffdfULL, 0xfefffffdfffffbffULL, 0xdfffffbfffff7fffULL, 0xfffff7ffffefffffULL,
    0xffefffffdfffffbfULL, 0xfdfffffbfffff7ffULL, 0xbfffff7ffffeffffULL, 0xffffefffffdfffffULL,
    0xffdfffffbfffff7fULL, 0xfbfffff7ffffefffULL, 0x7ffffefffffdffffULL, 0xffffdfffffbfffffULL,
    0xffbfffff7ffffeffULL, 0xf7ffffefffffdfffULL, 0xfffffdfffffbffffULL, 0xffffbfffff7ffffeULL,
    0xff7ffffefffffdffULL, 0xefffffdfffffbfffULL, 0xfffffbfffff7ffffULL, 0xffff7ffffefffffdULL,
    0xfefffffdfffffbffULL, 0xdfffffbfffff7fffULL, 0xfffff7ffffefffffULL, 0xfffefffffdfffffbULL,
    0xfdfffffbfffff7ffULL, 0xbfffff7ffffeffffULL, 0xffffefffffdfffffULL, 0xfffdfffffbfffff7ULL,
    0xfbfffff7ffffefffULL, 0x7ffffefffffdffffULL, 0xffffdfffffbfffffULL, 0xfffbfffff7ffffefULL,
    0xf7ffffefffffdfffULL, 0xfffffdfffffbffffULL, 0xffffbfffff7ffffeULL, 0xfff7ffffefffffdfULL,
    0xefffffdfffffbfffULL, 0xfffffbfffff7ffffULL, 0xffff7ffffefffffdULL, 0xffefffffdfffffbfULL,
    0xdfffffbfffff7fffULL, 0xfffff7ffffefffffULL, 0xfffefffffdfffffbULL, 0xffdfffffbfffff7fULL,
    0xbfffff7ffffeffffULL, 0xffffefffffdfffffULL, 0xfffdfffffbfffff7ULL, 0xffbfffff7ffffeffULL,
    0x7ffffefffffdffffULL, 0xffffdfffffbfffffULL, 0xfffbfffff7ffffefULL, 0xff7ffffefffffdffULL,
    0xfffffdfffffbffffULL, 0xffffbfffff7ffffeULL, 0xfff7ffffefffffdfULL, 0xfefffffdfffffbffULL,
    0xfffffbfffff7ffffULL, 0xffff7ffffefffffdULL, 0xffefffffdfffffbfULL, 0xfdfffffbfffff7ffULL,
    0xfffff7ffffefffffULL, 0xfffefffffdfffffbULL, 0xffdfffffbfffff7fULL, 0xfbfffff7ffffefffULL,
    0xffffefffffdfffffULL, 0xfffdfffffbfffff7ULL, 0xffbfffff7ffffeffULL, 0xf7ffffefffffdfffULL,
    0xffffdfffffbfffffULL, 0xfffbfffff7ffffefULL, 0xff7ffffefffffdffULL, 0xefffffdfffffbfffULL };

uint64 _29_MASKS_256[29 * 4] = {
    0xfbffffffdffffffeULL, 0xffefffffff7fffffULL, 0xffffbffffffdffffULL, 0xfffffefffffff7ffULL,
    0xf7ffffffbffffffdULL, 0xffdffffffeffffffULL, 0xffff7ffffffbffffULL, 0xfffffdffffffefffULL,
    0xefffffff7ffffffbULL, 0xffbffffffdffffffULL, 0xfffefffffff7ffffULL, 0xfffffbffffffdfffULL,
    0xdffffffefffffff7ULL, 0xff7ffffffbffffffULL, 0xfffdffffffefffffULL, 0xfffff7ffffffbfffULL,
    0xbffffffdffffffefULL, 0xfefffffff7ffffffULL, 0xfffbffffffdfffffULL, 0xffffefffffff7fffULL,
    0x7ffffffbffffffdfULL, 0xfdffffffefffffffULL, 0xfff7ffffffbfffffULL, 0xffffdffffffeffffULL,
    0xfffffff7ffffffbfULL, 0xfbffffffdffffffeULL, 0xffefffffff7fffffULL, 0xffffbffffffdffffULL,
    0xffffffefffffff7fULL, 0xf7ffffffbffffffdULL, 0xffdffffffeffffffULL, 0xffff7ffffffbffffULL,
    0xffffffdffffffeffULL, 0xefffffff7ffffffbULL, 0xffbffffffdffffffULL, 0xfffefffffff7ffffULL,
    0xffffffbffffffdffULL, 0xdffffffefffffff7ULL, 0xff7ffffffbffffffULL, 0xfffdffffffefffffULL,
    0xffffff7ffffffbffULL, 0xbffffffdffffffefULL, 0xfefffffff7ffffffULL, 0xfffbffffffdfffffULL,
    0xfffffefffffff7ffULL, 0x7ffffffbffffffdfULL, 0xfdffffffefffffffULL, 0xfff7ffffffbfffffULL,
    0xfffffdffffffefffULL, 0xfffffff7ffffffbfULL, 0xfbffffffdffffffeULL, 0xffefffffff7fffffULL,
    0xfffffbffffffdfffULL, 0xffffffefffffff7fULL, 0xf7ffffffbffffffdULL, 0xffdffffffeffffffULL,
    0xfffff7ffffffbfffULL, 0xffffffdffffffeffULL, 0xefffffff7ffffffbULL, 0xffbffffffdffffffULL,
    0xffffefffffff7fffULL, 0xffffffbffffffdffULL, 0xdffffffefffffff7ULL, 0xff7ffffffbffffffULL,
    0xffffdffffffeffffULL, 0xffffff7ffffffbffULL, 0xbffffffdffffffefULL, 0xfefffffff7ffffffULL,
    0xffffbffffffdffffULL, 0xfffffefffffff7ffULL, 0x7ffffffbffffffdfULL, 0xfdffffffefffffffULL,
    0xffff7ffffffbffffULL, 0xfffffdffffffefffULL, 0xfffffff7ffffffbfULL, 0xfbffffffdffffffeULL,
    0xfffefffffff7ffffULL, 0xfffffbffffffdfffULL, 0xffffffefffffff7fULL, 0xf7ffffffbffffffdULL,
    0xfffdffffffefffffULL, 0xfffff7ffffffbfffULL, 0xffffffdffffffeffULL, 0xefffffff7ffffffbULL,
    0xfffbffffffdfffffULL, 0xffffefffffff7fffULL, 0xffffffbffffffdffULL, 0xdffffffefffffff7ULL,
    0xfff7ffffffbfffffULL, 0xffffdffffffeffffULL, 0xffffff7ffffffbffULL, 0xbffffffdffffffefULL,
    0xffefffffff7fffffULL, 0xffffbffffffdffffULL, 0xfffffefffffff7ffULL, 0x7ffffffbffffffdfULL,
    0xffdffffffeffffffULL, 0xffff7ffffffbffffULL, 0xfffffdffffffefffULL, 0xfffffff7ffffffbfULL,
    0xffbffffffdffffffULL, 0xfffefffffff7ffffULL, 0xfffffbffffffdfffULL, 0xffffffefffffff7fULL,
    0xff7ffffffbffffffULL, 0xfffdffffffefffffULL, 0xfffff7ffffffbfffULL, 0xffffffdffffffeffULL,
    0xfefffffff7ffffffULL, 0xfffbffffffdfffffULL, 0xffffefffffff7fffULL, 0xffffffbffffffdffULL,
    0xfdffffffefffffffULL, 0xfff7ffffffbfffffULL, 0xffffdffffffeffffULL, 0xffffff7ffffffbffULL };

uint64 _31_MASKS_256[31 * 4] = {
    0xbfffffff7ffffffeULL, 0xefffffffdfffffffULL, 0xfbfffffff7ffffffULL, 0xfefffffffdffffffULL,
    0x7ffffffefffffffdULL, 0xdfffffffbfffffffULL, 0xf7ffffffefffffffULL, 0xfdfffffffbffffffULL,
    0xfffffffdfffffffbULL, 0xbfffffff7ffffffeULL, 0xefffffffdfffffffULL, 0xfbfffffff7ffffffULL,
    0xfffffffbfffffff7ULL, 0x7ffffffefffffffdULL, 0xdfffffffbfffffffULL, 0xf7ffffffefffffffULL,
    0xfffffff7ffffffefULL, 0xfffffffdfffffffbULL, 0xbfffffff7ffffffeULL, 0xefffffffdfffffffULL,
    0xffffffefffffffdfULL, 0xfffffffbfffffff7ULL, 0x7ffffffefffffffdULL, 0xdfffffffbfffffffULL,
    0xffffffdfffffffbfULL, 0xfffffff7ffffffefULL, 0xfffffffdfffffffbULL, 0xbfffffff7ffffffeULL,
    0xffffffbfffffff7fULL, 0xffffffefffffffdfULL, 0xfffffffbfffffff7ULL, 0x7ffffffefffffffdULL,
    0xffffff7ffffffeffULL, 0xffffffdfffffffbfULL, 0xfffffff7ffffffefULL, 0xfffffffdfffffffbULL,
    0xfffffefffffffdffULL, 0xffffffbfffffff7fULL, 0xffffffefffffffdfULL, 0xfffffffbfffffff7ULL,
    0xfffffdfffffffbffULL, 0xffffff7ffffffeffULL, 0xffffffdfffffffbfULL, 0xfffffff7ffffffefULL,
    0xfffffbfffffff7ffULL, 0xfffffefffffffdffULL, 0xffffffbfffffff7fULL, 0xffffffefffffffdfULL,
    0xfffff7ffffffefffULL, 0xfffffdfffffffbffULL, 0xffffff7ffffffeffULL, 0xffffffdfffffffbfULL,
    0xffffefffffffdfffULL, 0xfffffbfffffff7ffULL, 0xfffffefffffffdffULL, 0xffffffbfffffff7fULL,
    0xffffdfffffffbfffULL, 0xfffff7ffffffefffULL, 0xfffffdfffffffbffULL, 0xffffff7ffffffeffULL,
    0xffffbfffffff7fffULL, 0xffffefffffffdfffULL, 0xfffffbfffffff7ffULL, 0xfffffefffffffdffULL,
    0xffff7ffffffeffffULL, 0xffffdfffffffbfffULL, 0xfffff7ffffffefffULL, 0xfffffdfffffffbffULL,
    0xfffefffffffdffffULL, 0xffffbfffffff7fffULL, 0xffffefffffffdfffULL, 0xfffffbfffffff7ffULL,
    0xfffdfffffffbffffULL, 0xffff7ffffffeffffULL, 0xffffdfffffffbfffULL, 0xfffff7ffffffefffULL,
    0xfffbfffffff7ffffULL, 0xfffefffffffdffffULL, 0xffffbfffffff7fffULL, 0xffffefffffffdfffULL,
    0xfff7ffffffefffffULL, 0xfffdfffffffbffffULL, 0xffff7ffffffeffffULL, 0xffffdfffffffbfffULL,
    0xffefffffffdfffffULL, 0xfffbfffffff7ffffULL, 0xfffefffffffdffffULL, 0xffffbfffffff7fffULL,
    0xffdfffffffbfffffULL, 0xfff7ffffffefffffULL, 0xfffdfffffffbffffULL, 0xffff7ffffffeffffULL,
    0xffbfffffff7fffffULL, 0xffefffffffdfffffULL, 0xfffbfffffff7ffffULL, 0xfffefffffffdffffULL,
    0xff7ffffffeffffffULL, 0xffdfffffffbfffffULL, 0xfff7ffffffefffffULL, 0xfffdfffffffbffffULL,
    0xfefffffffdffffffULL, 0xffbfffffff7fffffULL, 0xffefffffffdfffffULL, 0xfffbfffffff7ffffULL,
    0xfdfffffffbffffffULL, 0xff7ffffffeffffffULL, 0xffdfffffffbfffffULL, 0xfff7ffffffefffffULL,
    0xfbfffffff7ffffffULL, 0xfefffffffdffffffULL, 0xffbfffffff7fffffULL, 0xffefffffffdfffffULL,
    0xf7ffffffefffffffULL, 0xfdfffffffbffffffULL, 0xff7ffffffeffffffULL, 0xffdfffffffbfffffULL,
    0xefffffffdfffffffULL, 0xfbfffffff7ffffffULL, 0xfefffffffdffffffULL, 0xffbfffffff7fffffULL,
    0xdfffffffbfffffffULL, 0xf7ffffffefffffffULL, 0xfdfffffffbffffffULL, 0xff7ffffffeffffffULL };

uint64 _37_MASKS_256[37 * 4] = {
    0xffffffdffffffffeULL, 0xffff7ffffffffbffULL, 0xfdffffffffefffffULL, 0xffffffffbfffffffULL,
    0xffffffbffffffffdULL, 0xfffefffffffff7ffULL, 0xfbffffffffdfffffULL, 0xffffffff7fffffffULL,
    0xffffff7ffffffffbULL, 0xfffdffffffffefffULL, 0xf7ffffffffbfffffULL, 0xfffffffeffffffffULL,
    0xfffffefffffffff7ULL, 0xfffbffffffffdfffULL, 0xefffffffff7fffffULL, 0xfffffffdffffffffULL,
    0xfffffdffffffffefULL, 0xfff7ffffffffbfffULL, 0xdffffffffeffffffULL, 0xfffffffbffffffffULL,
    0xfffffbffffffffdfULL, 0xffefffffffff7fffULL, 0xbffffffffdffffffULL, 0xfffffff7ffffffffULL,
    0xfffff7ffffffffbfULL, 0xffdffffffffeffffULL, 0x7ffffffffbffffffULL, 0xffffffefffffffffULL,
    0xffffefffffffff7fULL, 0xffbffffffffdffffULL, 0xfffffffff7ffffffULL, 0xffffffdffffffffeULL,
    0xffffdffffffffeffULL, 0xff7ffffffffbffffULL, 0xffffffffefffffffULL, 0xffffffbffffffffdULL,
    0xffffbffffffffdffULL, 0xfefffffffff7ffffULL, 0xffffffffdfffffffULL, 0xffffff7ffffffffbULL,
    0xffff7ffffffffbffULL, 0xfdffffffffefffffULL, 0xffffffffbfffffffULL, 0xfffffefffffffff7ULL,
    0xfffefffffffff7ffULL, 0xfbffffffffdfffffULL, 0xffffffff7fffffffULL, 0xfffffdffffffffefULL,
    0xfffdffffffffefffULL, 0xf7ffffffffbfffffULL, 0xfffffffeffffffffULL, 0xfffffbffffffffdfULL,
    0xfffbffffffffdfffULL, 0xefffffffff7fffffULL, 0xfffffffdffffffffULL, 0xfffff7ffffffffbfULL,
    0xfff7ffffffffbfffULL, 0xdffffffffeffffffULL, 0xfffffffbffffffffULL, 0xffffefffffffff7fULL,
    0xffefffffffff7fffULL, 0xbffffffffdffffffULL, 0xfffffff7ffffffffULL, 0xffffdffffffffeffULL,
    0xffdffffffffeffffULL, 0x7ffffffffbffffffULL, 0xffffffefffffffffULL, 0xffffbffffffffdffULL,
    0xffbffffffffdffffULL, 0xfffffffff7ffffffULL, 0xffffffdffffffffeULL, 0xffff7ffffffffbffULL,
    0xff7ffffffffbffffULL, 0xffffffffefffffffULL, 0xffffffbffffffffdULL, 0xfffefffffffff7ffULL,
    0xfefffffffff7ffffULL, 0xffffffffdfffffffULL, 0xffffff7ffffffffbULL, 0xfffdffffffffefffULL,
    0xfdffffffffefffffULL, 0xffffffffbfffffffULL, 0xfffffefffffffff7ULL, 0xfffbffffffffdfffULL,
    0xfbffffffffdfffffULL, 0xffffffff7fffffffULL, 0xfffffdffffffffefULL, 0xfff7ffffffffbfffULL,
    0xf7ffffffffbfffffULL, 0xfffffffeffffffffULL, 0xfffffbffffffffdfULL, 0xffefffffffff7fffULL,
    0xefffffffff7fffffULL, 0xfffffffdffffffffULL, 0xfffff7ffffffffbfULL, 0xffdffffffffeffffULL,
    0xdffffffffeffffffULL, 0xfffffffbffffffffULL, 0xffffefffffffff7fULL, 0xffbffffffffdffffULL,
    0xbffffffffdffffffULL, 0xfffffff7ffffffffULL, 0xffffdffffffffeffULL, 0xff7ffffffffbffffULL,
    0x7ffffffffbffffffULL, 0xffffffefffffffffULL, 0xffffbffffffffdffULL, 0xfefffffffff7ffffULL,
    0xfffffffff7ffffffULL, 0xffffffdffffffffeULL, 0xffff7ffffffffbffULL, 0xfdffffffffefffffULL,
    0xffffffffefffffffULL, 0xffffffbffffffffdULL, 0xfffefffffffff7ffULL, 0xfbffffffffdfffffULL,
    0xffffffffdfffffffULL, 0xffffff7ffffffffbULL, 0xfffdffffffffefffULL, 0xf7ffffffffbfffffULL,
    0xffffffffbfffffffULL, 0xfffffefffffffff7ULL, 0xfffbffffffffdfffULL, 0xefffffffff7fffffULL,
    0xffffffff7fffffffULL, 0xfffffdffffffffefULL, 0xfff7ffffffffbfffULL, 0xdffffffffeffffffULL,
    0xfffffffeffffffffULL, 0xfffffbffffffffdfULL, 0xffefffffffff7fffULL, 0xbffffffffdffffffULL,
    0xfffffffdffffffffULL, 0xfffff7ffffffffbfULL, 0xffdffffffffeffffULL, 0x7ffffffffbffffffULL,
    0xfffffffbffffffffULL, 0xffffefffffffff7fULL, 0xffbffffffffdffffULL, 0xfffffffff7ffffffULL,
    0xfffffff7ffffffffULL, 0xffffdffffffffeffULL, 0xff7ffffffffbffffULL, 0xffffffffefffffffULL,
    0xffffffefffffffffULL, 0xffffbffffffffdffULL, 0xfefffffffff7ffffULL, 0xffffffffdfffffffULL };

uint64 _41_MASKS_256[41 * 4] = {
    0xfffffdfffffffffeULL, 0xf7fffffffffbffffULL, 0xffffffefffffffffULL, 0xffbfffffffffdfffULL,
    0xfffffbfffffffffdULL, 0xeffffffffff7ffffULL, 0xffffffdfffffffffULL, 0xff7fffffffffbfffULL,
    0xfffff7fffffffffbULL, 0xdfffffffffefffffULL, 0xffffffbfffffffffULL, 0xfeffffffffff7fffULL,
    0xffffeffffffffff7ULL, 0xbfffffffffdfffffULL, 0xffffff7fffffffffULL, 0xfdfffffffffeffffULL,
    0xffffdfffffffffefULL, 0x7fffffffffbfffffULL, 0xfffffeffffffffffULL, 0xfbfffffffffdffffULL,
    0xffffbfffffffffdfULL, 0xffffffffff7fffffULL, 0xfffffdfffffffffeULL, 0xf7fffffffffbffffULL,
    0xffff7fffffffffbfULL, 0xfffffffffeffffffULL, 0xfffffbfffffffffdULL, 0xeffffffffff7ffffULL,
    0xfffeffffffffff7fULL, 0xfffffffffdffffffULL, 0xfffff7fffffffffbULL, 0xdfffffffffefffffULL,
    0xfffdfffffffffeffULL, 0xfffffffffbffffffULL, 0xffffeffffffffff7ULL, 0xbfffffffffdfffffULL,
    0xfffbfffffffffdffULL, 0xfffffffff7ffffffULL, 0xffffdfffffffffefULL, 0x7fffffffffbfffffULL,
    0xfff7fffffffffbffULL, 0xffffffffefffffffULL, 0xffffbfffffffffdfULL, 0xffffffffff7fffffULL,
    0xffeffffffffff7ffULL, 0xffffffffdfffffffULL, 0xffff7fffffffffbfULL, 0xfffffffffeffffffULL,
    0xffdfffffffffefffULL, 0xffffffffbfffffffULL, 0xfffeffffffffff7fULL, 0xfffffffffdffffffULL,
    0xffbfffffffffdfffULL, 0xffffffff7fffffffULL, 0xfffdfffffffffeffULL, 0xfffffffffbffffffULL,
    0xff7fffffffffbfffULL, 0xfffffffeffffffffULL, 0xfffbfffffffffdffULL, 0xfffffffff7ffffffULL,
    0xfeffffffffff7fffULL, 0xfffffffdffffffffULL, 0xfff7fffffffffbffULL, 0xffffffffefffffffULL,
    0xfdfffffffffeffffULL, 0xfffffffbffffffffULL, 0xffeffffffffff7ffULL, 0xffffffffdfffffffULL,
    0xfbfffffffffdffffULL, 0xfffffff7ffffffffULL, 0xffdfffffffffefffULL, 0xffffffffbfffffffULL,
    0xf7fffffffffbffffULL, 0xffffffefffffffffULL, 0xffbfffffffffdfffULL, 0xffffffff7fffffffULL,
    0xeffffffffff7ffffULL, 0xffffffdfffffffffULL, 0xff7fffffffffbfffULL, 0xfffffffeffffffffULL,
    0xdfffffffffefffffULL, 0xffffffbfffffffffULL, 0xfeffffffffff7fffULL, 0xfffffffdffffffffULL,
    0xbfffffffffdfffffULL, 0xffffff7fffffffffULL, 0xfdfffffffffeffffULL, 0xfffffffbffffffffULL,
    0x7fffffffffbfffffULL, 0xfffffeffffffffffULL, 0xfbfffffffffdffffULL, 0xfffffff7ffffffffULL,
    0xffffffffff7fffffULL, 0xfffffdfffffffffeULL, 0xf7fffffffffbffffULL, 0xffffffefffffffffULL,
    0xfffffffffeffffffULL, 0xfffffbfffffffffdULL, 0xeffffffffff7ffffULL, 0xffffffdfffffffffULL,
    0xfffffffffdffffffULL, 0xfffff7fffffffffbULL, 0xdfffffffffefffffULL, 0xffffffbfffffffffULL,
    0xfffffffffbffffffULL, 0xffffeffffffffff7ULL, 0xbfffffffffdfffffULL, 0xffffff7fffffffffULL,
    0xfffffffff7ffffffULL, 0xffffdfffffffffefULL, 0x7fffffffffbfffffULL, 0xfffffeffffffffffULL,
    0xffffffffefffffffULL, 0xffffbfffffffffdfULL, 0xffffffffff7fffffULL, 0xfffffdfffffffffeULL,
    0xffffffffdfffffffULL, 0xffff7fffffffffbfULL, 0xfffffffffeffffffULL, 0xfffffbfffffffffdULL,
    0xffffffffbfffffffULL, 0xfffeffffffffff7fULL, 0xfffffffffdffffffULL, 0xfffff7fffffffffbULL,
    0xffffffff7fffffffULL, 0xfffdfffffffffeffULL, 0xfffffffffbffffffULL, 0xffffeffffffffff7ULL,
    0xfffffffeffffffffULL, 0xfffbfffffffffdffULL, 0xfffffffff7ffffffULL, 0xffffdfffffffffefULL,
    0xfffffffdffffffffULL, 0xfff7fffffffffbffULL, 0xffffffffefffffffULL, 0xffffbfffffffffdfULL,
    0xfffffffbffffffffULL, 0xffeffffffffff7ffULL, 0xffffffffdfffffffULL, 0xffff7fffffffffbfULL,
    0xfffffff7ffffffffULL, 0xffdfffffffffefffULL, 0xffffffffbfffffffULL, 0xfffeffffffffff7fULL,
    0xffffffefffffffffULL, 0xffbfffffffffdfffULL, 0xffffffff7fffffffULL, 0xfffdfffffffffeffULL,
    0xffffffdfffffffffULL, 0xff7fffffffffbfffULL, 0xfffffffeffffffffULL, 0xfffbfffffffffdffULL,
    0xffffffbfffffffffULL, 0xfeffffffffff7fffULL, 0xfffffffdffffffffULL, 0xfff7fffffffffbffULL,
    0xffffff7fffffffffULL, 0xfdfffffffffeffffULL, 0xfffffffbffffffffULL, 0xffeffffffffff7ffULL,
    0xfffffeffffffffffULL, 0xfbfffffffffdffffULL, 0xfffffff7ffffffffULL, 0xffdfffffffffefffULL };

uint64 _43_MASKS_256[43 * 4] = {
    0xfffff7fffffffffeULL, 0xffffffffffbfffffULL, 0xffffeffffffffffdULL, 0xffffffffff7fffffULL,
    0xffffeffffffffffdULL, 0xffffffffff7fffffULL, 0xffffdffffffffffbULL, 0xfffffffffeffffffULL,
    0xffffdffffffffffbULL, 0xfffffffffeffffffULL, 0xffffbffffffffff7ULL, 0xfffffffffdffffffULL,
    0xffffbffffffffff7ULL, 0xfffffffffdffffffULL, 0xffff7fffffffffefULL, 0xfffffffffbffffffULL,
    0xffff7fffffffffefULL, 0xfffffffffbffffffULL, 0xfffeffffffffffdfULL, 0xfffffffff7ffffffULL,
    0xfffeffffffffffdfULL, 0xfffffffff7ffffffULL, 0xfffdffffffffffbfULL, 0xffffffffefffffffULL,
    0xfffdffffffffffbfULL, 0xffffffffefffffffULL, 0xfffbffffffffff7fULL, 0xffffffffdfffffffULL,
    0xfffbffffffffff7fULL, 0xffffffffdfffffffULL, 0xfff7fffffffffeffULL, 0xffffffffbfffffffULL,
    0xfff7fffffffffeffULL, 0xffffffffbfffffffULL, 0xffeffffffffffdffULL, 0xffffffff7fffffffULL,
    0xffeffffffffffdffULL, 0xffffffff7fffffffULL, 0xffdffffffffffbffULL, 0xfffffffeffffffffULL,
    0xffdffffffffffbffULL, 0xfffffffeffffffffULL, 0xffbffffffffff7ffULL, 0xfffffffdffffffffULL,
    0xffbffffffffff7ffULL, 0xfffffffdffffffffULL, 0xff7fffffffffefffULL, 0xfffffffbffffffffULL,
    0xff7fffffffffefffULL, 0xfffffffbffffffffULL, 0xfeffffffffffdfffULL, 0xfffffff7ffffffffULL,
    0xfeffffffffffdfffULL, 0xfffffff7ffffffffULL, 0xfdffffffffffbfffULL, 0xffffffefffffffffULL,
    0xfdffffffffffbfffULL, 0xffffffefffffffffULL, 0xfbffffffffff7fffULL, 0xffffffdfffffffffULL,
    0xfbffffffffff7fffULL, 0xffffffdfffffffffULL, 0xf7fffffffffeffffULL, 0xffffffbfffffffffULL,
    0xf7fffffffffeffffULL, 0xffffffbfffffffffULL, 0xeffffffffffdffffULL, 0xffffff7fffffffffULL,
    0xeffffffffffdffffULL, 0xffffff7fffffffffULL, 0xdffffffffffbffffULL, 0xfffffeffffffffffULL,
    0xdffffffffffbffffULL, 0xfffffeffffffffffULL, 0xbffffffffff7ffffULL, 0xfffffdffffffffffULL,
    0xbffffffffff7ffffULL, 0xfffffdffffffffffULL, 0x7fffffffffefffffULL, 0xfffffbffffffffffULL,
    0x7fffffffffefffffULL, 0xfffffbffffffffffULL, 0xffffffffffdfffffULL, 0xfffff7fffffffffeULL,
    0xffffffffffdfffffULL, 0xfffff7fffffffffeULL, 0xffffffffffbfffffULL, 0xffffeffffffffffdULL,
    0xffffffffffbfffffULL, 0xffffeffffffffffdULL, 0xffffffffff7fffffULL, 0xffffdffffffffffbULL,
    0xffffffffff7fffffULL, 0xffffdffffffffffbULL, 0xfffffffffeffffffULL, 0xffffbffffffffff7ULL,
    0xfffffffffeffffffULL, 0xffffbffffffffff7ULL, 0xfffffffffdffffffULL, 0xffff7fffffffffefULL,
    0xfffffffffdffffffULL, 0xffff7fffffffffefULL, 0xfffffffffbffffffULL, 0xfffeffffffffffdfULL,
    0xfffffffffbffffffULL, 0xfffeffffffffffdfULL, 0xfffffffff7ffffffULL, 0xfffdffffffffffbfULL,
    0xfffffffff7ffffffULL, 0xfffdffffffffffbfULL, 0xffffffffefffffffULL, 0xfffbffffffffff7fULL,
    0xffffffffefffffffULL, 0xfffbffffffffff7fULL, 0xffffffffdfffffffULL, 0xfff7fffffffffeffULL,
    0xffffffffdfffffffULL, 0xfff7fffffffffeffULL, 0xffffffffbfffffffULL, 0xffeffffffffffdffULL,
    0xffffffffbfffffffULL, 0xffeffffffffffdffULL, 0xffffffff7fffffffULL, 0xffdffffffffffbffULL,
    0xffffffff7fffffffULL, 0xffdffffffffffbffULL, 0xfffffffeffffffffULL, 0xffbffffffffff7ffULL,
    0xfffffffeffffffffULL, 0xffbffffffffff7ffULL, 0xfffffffdffffffffULL, 0xff7fffffffffefffULL,
    0xfffffffdffffffffULL, 0xff7fffffffffefffULL, 0xfffffffbffffffffULL, 0xfeffffffffffdfffULL,
    0xfffffffbffffffffULL, 0xfeffffffffffdfffULL, 0xfffffff7ffffffffULL, 0xfdffffffffffbfffULL,
    0xfffffff7ffffffffULL, 0xfdffffffffffbfffULL, 0xffffffefffffffffULL, 0xfbffffffffff7fffULL,
    0xffffffefffffffffULL, 0xfbffffffffff7fffULL, 0xffffffdfffffffffULL, 0xf7fffffffffeffffULL,
    0xffffffdfffffffffULL, 0xf7fffffffffeffffULL, 0xffffffbfffffffffULL, 0xeffffffffffdffffULL,
    0xffffffbfffffffffULL, 0xeffffffffffdffffULL, 0xffffff7fffffffffULL, 0xdffffffffffbffffULL,
    0xffffff7fffffffffULL, 0xdffffffffffbffffULL, 0xfffffeffffffffffULL, 0xbffffffffff7ffffULL,
    0xfffffeffffffffffULL, 0xbffffffffff7ffffULL, 0xfffffdffffffffffULL, 0x7fffffffffefffffULL,
    0xfffffdffffffffffULL, 0x7fffffffffefffffULL, 0xfffffbffffffffffULL, 0xffffffffffdfffffULL,
    0xfffffbffffffffffULL, 0xffffffffffdfffffULL, 0xfffff7fffffffffeULL, 0xffffffffffbfffffULL };

uint64 _47_MASKS_256[47 * 4] = {
    0xffff7ffffffffffeULL, 0xffffffffbfffffffULL, 0xefffffffffffdfffULL, 0xfffff7ffffffffffULL,
    0xfffefffffffffffdULL, 0xffffffff7fffffffULL, 0xdfffffffffffbfffULL, 0xffffefffffffffffULL,
    0xfffdfffffffffffbULL, 0xfffffffeffffffffULL, 0xbfffffffffff7fffULL, 0xffffdfffffffffffULL,
    0xfffbfffffffffff7ULL, 0xfffffffdffffffffULL, 0x7ffffffffffeffffULL, 0xffffbfffffffffffULL,
    0xfff7ffffffffffefULL, 0xfffffffbffffffffULL, 0xfffffffffffdffffULL, 0xffff7ffffffffffeULL,
    0xffefffffffffffdfULL, 0xfffffff7ffffffffULL, 0xfffffffffffbffffULL, 0xfffefffffffffffdULL,
    0xffdfffffffffffbfULL, 0xffffffefffffffffULL, 0xfffffffffff7ffffULL, 0xfffdfffffffffffbULL,
    0xffbfffffffffff7fULL, 0xffffffdfffffffffULL, 0xffffffffffefffffULL, 0xfffbfffffffffff7ULL,
    0xff7ffffffffffeffULL, 0xffffffbfffffffffULL, 0xffffffffffdfffffULL, 0xfff7ffffffffffefULL,
    0xfefffffffffffdffULL, 0xffffff7fffffffffULL, 0xffffffffffbfffffULL, 0xffefffffffffffdfULL,
    0xfdfffffffffffbffULL, 0xfffffeffffffffffULL, 0xffffffffff7fffffULL, 0xffdfffffffffffbfULL,
    0xfbfffffffffff7ffULL, 0xfffffdffffffffffULL, 0xfffffffffeffffffULL, 0xffbfffffffffff7fULL,
    0xf7ffffffffffefffULL, 0xfffffbffffffffffULL, 0xfffffffffdffffffULL, 0xff7ffffffffffeffULL,
    0xefffffffffffdfffULL, 0xfffff7ffffffffffULL, 0xfffffffffbffffffULL, 0xfefffffffffffdffULL,
    0xdfffffffffffbfffULL, 0xffffefffffffffffULL, 0xfffffffff7ffffffULL, 0xfdfffffffffffbffULL,
    0xbfffffffffff7fffULL, 0xffffdfffffffffffULL, 0xffffffffefffffffULL, 0xfbfffffffffff7ffULL,
    0x7ffffffffffeffffULL, 0xffffbfffffffffffULL, 0xffffffffdfffffffULL, 0xf7ffffffffffefffULL,
    0xfffffffffffdffffULL, 0xffff7ffffffffffeULL, 0xffffffffbfffffffULL, 0xefffffffffffdfffULL,
    0xfffffffffffbffffULL, 0xfffefffffffffffdULL, 0xffffffff7fffffffULL, 0xdfffffffffffbfffULL,
    0xfffffffffff7ffffULL, 0xfffdfffffffffffbULL, 0xfffffffeffffffffULL, 0xbfffffffffff7fffULL,
    0xffffffffffefffffULL, 0xfffbfffffffffff7ULL, 0xfffffffdffffffffULL, 0x7ffffffffffeffffULL,
    0xffffffffffdfffffULL, 0xfff7ffffffffffefULL, 0xfffffffbffffffffULL, 0xfffffffffffdffffULL,
    0xffffffffffbfffffULL, 0xffefffffffffffdfULL, 0xfffffff7ffffffffULL, 0xfffffffffffbffffULL,
    0xffffffffff7fffffULL, 0xffdfffffffffffbfULL, 0xffffffefffffffffULL, 0xfffffffffff7ffffULL,
    0xfffffffffeffffffULL, 0xffbfffffffffff7fULL, 0xffffffdfffffffffULL, 0xffffffffffefffffULL,
    0xfffffffffdffffffULL, 0xff7ffffffffffeffULL, 0xffffffbfffffffffULL, 0xffffffffffdfffffULL,
    0xfffffffffbffffffULL, 0xfefffffffffffdffULL, 0xffffff7fffffffffULL, 0xffffffffffbfffffULL,
    0xfffffffff7ffffffULL, 0xfdfffffffffffbffULL, 0xfffffeffffffffffULL, 0xffffffffff7fffffULL,
    0xffffffffefffffffULL, 0xfbfffffffffff7ffULL, 0xfffffdffffffffffULL, 0xfffffffffeffffffULL,
    0xffffffffdfffffffULL, 0xf7ffffffffffefffULL, 0xfffffbffffffffffULL, 0xfffffffffdffffffULL,
    0xffffffffbfffffffULL, 0xefffffffffffdfffULL, 0xfffff7ffffffffffULL, 0xfffffffffbffffffULL,
    0xffffffff7fffffffULL, 0xdfffffffffffbfffULL, 0xffffefffffffffffULL, 0xfffffffff7ffffffULL,
    0xfffffffeffffffffULL, 0xbfffffffffff7fffULL, 0xffffdfffffffffffULL, 0xffffffffefffffffULL,
    0xfffffffdffffffffULL, 0x7ffffffffffeffffULL, 0xffffbfffffffffffULL, 0xffffffffdfffffffULL,
    0xfffffffbffffffffULL, 0xfffffffffffdffffULL, 0xffff7ffffffffffeULL, 0xffffffffbfffffffULL,
    0xfffffff7ffffffffULL, 0xfffffffffffbffffULL, 0xfffefffffffffffdULL, 0xffffffff7fffffffULL,
    0xffffffefffffffffULL, 0xfffffffffff7ffffULL, 0xfffdfffffffffffbULL, 0xfffffffeffffffffULL,
    0xffffffdfffffffffULL, 0xffffffffffefffffULL, 0xfffbfffffffffff7ULL, 0xfffffffdffffffffULL,
    0xffffffbfffffffffULL, 0xffffffffffdfffffULL, 0xfff7ffffffffffefULL, 0xfffffffbffffffffULL,
    0xffffff7fffffffffULL, 0xffffffffffbfffffULL, 0xffefffffffffffdfULL, 0xfffffff7ffffffffULL,
    0xfffffeffffffffffULL, 0xffffffffff7fffffULL, 0xffdfffffffffffbfULL, 0xffffffefffffffffULL,
    0xfffffdffffffffffULL, 0xfffffffffeffffffULL, 0xffbfffffffffff7fULL, 0xffffffdfffffffffULL,
    0xfffffbffffffffffULL, 0xfffffffffdffffffULL, 0xff7ffffffffffeffULL, 0xffffffbfffffffffULL,
    0xfffff7ffffffffffULL, 0xfffffffffbffffffULL, 0xfefffffffffffdffULL, 0xffffff7fffffffffULL,
    0xffffefffffffffffULL, 0xfffffffff7ffffffULL, 0xfdfffffffffffbffULL, 0xfffffeffffffffffULL,
    0xffffdfffffffffffULL, 0xffffffffefffffffULL, 0xfbfffffffffff7ffULL, 0xfffffdffffffffffULL,
    0xffffbfffffffffffULL, 0xffffffffdfffffffULL, 0xf7ffffffffffefffULL, 0xfffffbffffffffffULL };

uint64 _53_MASKS_256[53 * 4] = {
    0xffdffffffffffffeULL, 0xfffffbffffffffffULL, 0xffffffff7fffffffULL, 0xffffffffffefffffULL,
    0xffbffffffffffffdULL, 0xfffff7ffffffffffULL, 0xfffffffeffffffffULL, 0xffffffffffdfffffULL,
    0xff7ffffffffffffbULL, 0xffffefffffffffffULL, 0xfffffffdffffffffULL, 0xffffffffffbfffffULL,
    0xfefffffffffffff7ULL, 0xffffdfffffffffffULL, 0xfffffffbffffffffULL, 0xffffffffff7fffffULL,
    0xfdffffffffffffefULL, 0xffffbfffffffffffULL, 0xfffffff7ffffffffULL, 0xfffffffffeffffffULL,
    0xfbffffffffffffdfULL, 0xffff7fffffffffffULL, 0xffffffefffffffffULL, 0xfffffffffdffffffULL,
    0xf7ffffffffffffbfULL, 0xfffeffffffffffffULL, 0xffffffdfffffffffULL, 0xfffffffffbffffffULL,
    0xefffffffffffff7fULL, 0xfffdffffffffffffULL, 0xffffffbfffffffffULL, 0xfffffffff7ffffffULL,
    0xdffffffffffffeffULL, 0xfffbffffffffffffULL, 0xffffff7fffffffffULL, 0xffffffffefffffffULL,
    0xbffffffffffffdffULL, 0xfff7ffffffffffffULL, 0xfffffeffffffffffULL, 0xffffffffdfffffffULL,
    0x7ffffffffffffbffULL, 0xffefffffffffffffULL, 0xfffffdffffffffffULL, 0xffffffffbfffffffULL,
    0xfffffffffffff7ffULL, 0xffdffffffffffffeULL, 0xfffffbffffffffffULL, 0xffffffff7fffffffULL,
    0xffffffffffffefffULL, 0xffbffffffffffffdULL, 0xfffff7ffffffffffULL, 0xfffffffeffffffffULL,
    0xffffffffffffdfffULL, 0xff7ffffffffffffbULL, 0xffffefffffffffffULL, 0xfffffffdffffffffULL,
    0xffffffffffffbfffULL, 0xfefffffffffffff7ULL, 0xffffdfffffffffffULL, 0xfffffffbffffffffULL,
    0xffffffffffff7fffULL, 0xfdffffffffffffefULL, 0xffffbfffffffffffULL, 0xfffffff7ffffffffULL,
    0xfffffffffffeffffULL, 0xfbffffffffffffdfULL, 0xffff7fffffffffffULL, 0xffffffefffffffffULL,
    0xfffffffffffdffffULL, 0xf7ffffffffffffbfULL, 0xfffeffffffffffffULL, 0xffffffdfffffffffULL,
    0xfffffffffffbffffULL, 0xefffffffffffff7fULL, 0xfffdffffffffffffULL, 0xffffffbfffffffffULL,
    0xfffffffffff7ffffULL, 0xdffffffffffffeffULL, 0xfffbffffffffffffULL, 0xffffff7fffffffffULL,
    0xffffffffffefffffULL, 0xbffffffffffffdffULL, 0xfff7ffffffffffffULL, 0xfffffeffffffffffULL,
    0xffffffffffdfffffULL, 0x7ffffffffffffbffULL, 0xffefffffffffffffULL, 0xfffffdffffffffffULL,
    0xffffffffffbfffffULL, 0xfffffffffffff7ffULL, 0xffdffffffffffffeULL, 0xfffffbffffffffffULL,
    0xffffffffff7fffffULL, 0xffffffffffffefffULL, 0xffbffffffffffffdULL, 0xfffff7ffffffffffULL,
    0xfffffffffeffffffULL, 0xffffffffffffdfffULL, 0xff7ffffffffffffbULL, 0xffffefffffffffffULL,
    0xfffffffffdffffffULL, 0xffffffffffffbfffULL, 0xfefffffffffffff7ULL, 0xffffdfffffffffffULL,
    0xfffffffffbffffffULL, 0xffffffffffff7fffULL, 0xfdffffffffffffefULL, 0xffffbfffffffffffULL,
    0xfffffffff7ffffffULL, 0xfffffffffffeffffULL, 0xfbffffffffffffdfULL, 0xffff7fffffffffffULL,
    0xffffffffefffffffULL, 0xfffffffffffdffffULL, 0xf7ffffffffffffbfULL, 0xfffeffffffffffffULL,
    0xffffffffdfffffffULL, 0xfffffffffffbffffULL, 0xefffffffffffff7fULL, 0xfffdffffffffffffULL,
    0xffffffffbfffffffULL, 0xfffffffffff7ffffULL, 0xdffffffffffffeffULL, 0xfffbffffffffffffULL,
    0xffffffff7fffffffULL, 0xffffffffffefffffULL, 0xbffffffffffffdffULL, 0xfff7ffffffffffffULL,
    0xfffffffeffffffffULL, 0xffffffffffdfffffULL, 0x7ffffffffffffbffULL, 0xffefffffffffffffULL,
    0xfffffffdffffffffULL, 0xffffffffffbfffffULL, 0xfffffffffffff7ffULL, 0xffdffffffffffffeULL,
    0xfffffffbffffffffULL, 0xffffffffff7fffffULL, 0xffffffffffffefffULL, 0xffbffffffffffffdULL,
    0xfffffff7ffffffffULL, 0xfffffffffeffffffULL, 0xffffffffffffdfffULL, 0xff7ffffffffffffbULL,
    0xffffffefffffffffULL, 0xfffffffffdffffffULL, 0xffffffffffffbfffULL, 0xfefffffffffffff7ULL,
    0xffffffdfffffffffULL, 0xfffffffffbffffffULL, 0xffffffffffff7fffULL, 0xfdffffffffffffefULL,
    0xffffffbfffffffffULL, 0xfffffffff7ffffffULL, 0xfffffffffffeffffULL, 0xfbffffffffffffdfULL,
    0xffffff7fffffffffULL, 0xffffffffefffffffULL, 0xfffffffffffdffffULL, 0xf7ffffffffffffbfULL,
    0xfffffeffffffffffULL, 0xffffffffdfffffffULL, 0xfffffffffffbffffULL, 0xefffffffffffff7fULL,
    0xfffffdffffffffffULL, 0xffffffffbfffffffULL, 0xfffffffffff7ffffULL, 0xdffffffffffffeffULL,
    0xfffffbffffffffffULL, 0xffffffff7fffffffULL, 0xffffffffffefffffULL, 0xbffffffffffffdffULL,
    0xfffff7ffffffffffULL, 0xfffffffeffffffffULL, 0xffffffffffdfffffULL, 0x7ffffffffffffbffULL,
    0xffffefffffffffffULL, 0xfffffffdffffffffULL, 0xffffffffffbfffffULL, 0xfffffffffffff7ffULL,
    0xffffdfffffffffffULL, 0xfffffffbffffffffULL, 0xffffffffff7fffffULL, 0xffffffffffffefffULL,
    0xffffbfffffffffffULL, 0xfffffff7ffffffffULL, 0xfffffffffeffffffULL, 0xffffffffffffdfffULL,
    0xffff7fffffffffffULL, 0xffffffefffffffffULL, 0xfffffffffdffffffULL, 0xffffffffffffbfffULL,
    0xfffeffffffffffffULL, 0xffffffdfffffffffULL, 0xfffffffffbffffffULL, 0xffffffffffff7fffULL,
    0xfffdffffffffffffULL, 0xffffffbfffffffffULL, 0xfffffffff7ffffffULL, 0xfffffffffffeffffULL,
    0xfffbffffffffffffULL, 0xffffff7fffffffffULL, 0xffffffffefffffffULL, 0xfffffffffffdffffULL,
    0xfff7ffffffffffffULL, 0xfffffeffffffffffULL, 0xffffffffdfffffffULL, 0xfffffffffffbffffULL,
    0xffefffffffffffffULL, 0xfffffdffffffffffULL, 0xffffffffbfffffffULL, 0xfffffffffff7ffffULL };

uint64 _59_MASKS_256[59 * 4] = {
    0xf7fffffffffffffeULL, 0xffbfffffffffffffULL, 0xfffdffffffffffffULL, 0xffffefffffffffffULL,
    0xeffffffffffffffdULL, 0xff7fffffffffffffULL, 0xfffbffffffffffffULL, 0xffffdfffffffffffULL,
    0xdffffffffffffffbULL, 0xfeffffffffffffffULL, 0xfff7ffffffffffffULL, 0xffffbfffffffffffULL,
    0xbffffffffffffff7ULL, 0xfdffffffffffffffULL, 0xffefffffffffffffULL, 0xffff7fffffffffffULL,
    0x7fffffffffffffefULL, 0xfbffffffffffffffULL, 0xffdfffffffffffffULL, 0xfffeffffffffffffULL,
    0xffffffffffffffdfULL, 0xf7fffffffffffffeULL, 0xffbfffffffffffffULL, 0xfffdffffffffffffULL,
    0xffffffffffffffbfULL, 0xeffffffffffffffdULL, 0xff7fffffffffffffULL, 0xfffbffffffffffffULL,
    0xffffffffffffff7fULL, 0xdffffffffffffffbULL, 0xfeffffffffffffffULL, 0xfff7ffffffffffffULL,
    0xfffffffffffffeffULL, 0xbffffffffffffff7ULL, 0xfdffffffffffffffULL, 0xffefffffffffffffULL,
    0xfffffffffffffdffULL, 0x7fffffffffffffefULL, 0xfbffffffffffffffULL, 0xffdfffffffffffffULL,
    0xfffffffffffffbffULL, 0xffffffffffffffdfULL, 0xf7fffffffffffffeULL, 0xffbfffffffffffffULL,
    0xfffffffffffff7ffULL, 0xffffffffffffffbfULL, 0xeffffffffffffffdULL, 0xff7fffffffffffffULL,
    0xffffffffffffefffULL, 0xffffffffffffff7fULL, 0xdffffffffffffffbULL, 0xfeffffffffffffffULL,
    0xffffffffffffdfffULL, 0xfffffffffffffeffULL, 0xbffffffffffffff7ULL, 0xfdffffffffffffffULL,
    0xffffffffffffbfffULL, 0xfffffffffffffdffULL, 0x7fffffffffffffefULL, 0xfbffffffffffffffULL,
    0xffffffffffff7fffULL, 0xfffffffffffffbffULL, 0xffffffffffffffdfULL, 0xf7fffffffffffffeULL,
    0xfffffffffffeffffULL, 0xfffffffffffff7ffULL, 0xffffffffffffffbfULL, 0xeffffffffffffffdULL,
    0xfffffffffffdffffULL, 0xffffffffffffefffULL, 0xffffffffffffff7fULL, 0xdffffffffffffffbULL,
    0xfffffffffffbffffULL, 0xffffffffffffdfffULL, 0xfffffffffffffeffULL, 0xbffffffffffffff7ULL,
    0xfffffffffff7ffffULL, 0xffffffffffffbfffULL, 0xfffffffffffffdffULL, 0x7fffffffffffffefULL,
    0xffffffffffefffffULL, 0xffffffffffff7fffULL, 0xfffffffffffffbffULL, 0xffffffffffffffdfULL,
    0xffffffffffdfffffULL, 0xfffffffffffeffffULL, 0xfffffffffffff7ffULL, 0xffffffffffffffbfULL,
    0xffffffffffbfffffULL, 0xfffffffffffdffffULL, 0xffffffffffffefffULL, 0xffffffffffffff7fULL,
    0xffffffffff7fffffULL, 0xfffffffffffbffffULL, 0xffffffffffffdfffULL, 0xfffffffffffffeffULL,
    0xfffffffffeffffffULL, 0xfffffffffff7ffffULL, 0xffffffffffffbfffULL, 0xfffffffffffffdffULL,
    0xfffffffffdffffffULL, 0xffffffffffefffffULL, 0xffffffffffff7fffULL, 0xfffffffffffffbffULL,
    0xfffffffffbffffffULL, 0xffffffffffdfffffULL, 0xfffffffffffeffffULL, 0xfffffffffffff7ffULL,
    0xfffffffff7ffffffULL, 0xffffffffffbfffffULL, 0xfffffffffffdffffULL, 0xffffffffffffefffULL,
    0xffffffffefffffffULL, 0xffffffffff7fffffULL, 0xfffffffffffbffffULL, 0xffffffffffffdfffULL,
    0xffffffffdfffffffULL, 0xfffffffffeffffffULL, 0xfffffffffff7ffffULL, 0xffffffffffffbfffULL,
    0xffffffffbfffffffULL, 0xfffffffffdffffffULL, 0xffffffffffefffffULL, 0xffffffffffff7fffULL,
    0xffffffff7fffffffULL, 0xfffffffffbffffffULL, 0xffffffffffdfffffULL, 0xfffffffffffeffffULL,
    0xfffffffeffffffffULL, 0xfffffffff7ffffffULL, 0xffffffffffbfffffULL, 0xfffffffffffdffffULL,
    0xfffffffdffffffffULL, 0xffffffffefffffffULL, 0xffffffffff7fffffULL, 0xfffffffffffbffffULL,
    0xfffffffbffffffffULL, 0xffffffffdfffffffULL, 0xfffffffffeffffffULL, 0xfffffffffff7ffffULL,
    0xfffffff7ffffffffULL, 0xffffffffbfffffffULL, 0xfffffffffdffffffULL, 0xffffffffffefffffULL,
    0xffffffefffffffffULL, 0xffffffff7fffffffULL, 0xfffffffffbffffffULL, 0xffffffffffdfffffULL,
    0xffffffdfffffffffULL, 0xfffffffeffffffffULL, 0xfffffffff7ffffffULL, 0xffffffffffbfffffULL,
    0xffffffbfffffffffULL, 0xfffffffdffffffffULL, 0xffffffffefffffffULL, 0xffffffffff7fffffULL,
    0xffffff7fffffffffULL, 0xfffffffbffffffffULL, 0xffffffffdfffffffULL, 0xfffffffffeffffffULL,
    0xfffffeffffffffffULL, 0xfffffff7ffffffffULL, 0xffffffffbfffffffULL, 0xfffffffffdffffffULL,
    0xfffffdffffffffffULL, 0xffffffefffffffffULL, 0xffffffff7fffffffULL, 0xfffffffffbffffffULL,
    0xfffffbffffffffffULL, 0xffffffdfffffffffULL, 0xfffffffeffffffffULL, 0xfffffffff7ffffffULL,
    0xfffff7ffffffffffULL, 0xffffffbfffffffffULL, 0xfffffffdffffffffULL, 0xffffffffefffffffULL,
    0xffffefffffffffffULL, 0xffffff7fffffffffULL, 0xfffffffbffffffffULL, 0xffffffffdfffffffULL,
    0xffffdfffffffffffULL, 0xfffffeffffffffffULL, 0xfffffff7ffffffffULL, 0xffffffffbfffffffULL,
    0xffffbfffffffffffULL, 0xfffffdffffffffffULL, 0xffffffefffffffffULL, 0xffffffff7fffffffULL,
    0xffff7fffffffffffULL, 0xfffffbffffffffffULL, 0xffffffdfffffffffULL, 0xfffffffeffffffffULL,
    0xfffeffffffffffffULL, 0xfffff7ffffffffffULL, 0xffffffbfffffffffULL, 0xfffffffdffffffffULL,
    0xfffdffffffffffffULL, 0xffffefffffffffffULL, 0xffffff7fffffffffULL, 0xfffffffbffffffffULL,
    0xfffbffffffffffffULL, 0xffffdfffffffffffULL, 0xfffffeffffffffffULL, 0xfffffff7ffffffffULL,
    0xfff7ffffffffffffULL, 0xffffbfffffffffffULL, 0xfffffdffffffffffULL, 0xffffffefffffffffULL,
    0xffefffffffffffffULL, 0xffff7fffffffffffULL, 0xfffffbffffffffffULL, 0xffffffdfffffffffULL,
    0xffdfffffffffffffULL, 0xfffeffffffffffffULL, 0xfffff7ffffffffffULL, 0xffffffbfffffffffULL,
    0xffbfffffffffffffULL, 0xfffdffffffffffffULL, 0xffffefffffffffffULL, 0xffffff7fffffffffULL,
    0xff7fffffffffffffULL, 0xfffbffffffffffffULL, 0xffffdfffffffffffULL, 0xfffffeffffffffffULL,
    0xfeffffffffffffffULL, 0xfff7ffffffffffffULL, 0xffffbfffffffffffULL, 0xfffffdffffffffffULL,
    0xfdffffffffffffffULL, 0xffefffffffffffffULL, 0xffff7fffffffffffULL, 0xfffffbffffffffffULL,
    0xfbffffffffffffffULL, 0xffdfffffffffffffULL, 0xfffeffffffffffffULL, 0xfffff7ffffffffffULL };

uint64 _61_MASKS_256[61 * 4] = {
    0xdffffffffffffffeULL, 0xfbffffffffffffffULL, 0xff7fffffffffffffULL, 0xffefffffffffffffULL,
    0xbffffffffffffffdULL, 0xf7ffffffffffffffULL, 0xfeffffffffffffffULL, 0xffdfffffffffffffULL,
    0x7ffffffffffffffbULL, 0xefffffffffffffffULL, 0xfdffffffffffffffULL, 0xffbfffffffffffffULL,
    0xfffffffffffffff7ULL, 0xdffffffffffffffeULL, 0xfbffffffffffffffULL, 0xff7fffffffffffffULL,
    0xffffffffffffffefULL, 0xbffffffffffffffdULL, 0xf7ffffffffffffffULL, 0xfeffffffffffffffULL,
    0xffffffffffffffdfULL, 0x7ffffffffffffffbULL, 0xefffffffffffffffULL, 0xfdffffffffffffffULL,
    0xffffffffffffffbfULL, 0xfffffffffffffff7ULL, 0xdffffffffffffffeULL, 0xfbffffffffffffffULL,
    0xffffffffffffff7fULL, 0xffffffffffffffefULL, 0xbffffffffffffffdULL, 0xf7ffffffffffffffULL,
    0xfffffffffffffeffULL, 0xffffffffffffffdfULL, 0x7ffffffffffffffbULL, 0xefffffffffffffffULL,
    0xfffffffffffffdffULL, 0xffffffffffffffbfULL, 0xfffffffffffffff7ULL, 0xdffffffffffffffeULL,
    0xfffffffffffffbffULL, 0xffffffffffffff7fULL, 0xffffffffffffffefULL, 0xbffffffffffffffdULL,
    0xfffffffffffff7ffULL, 0xfffffffffffffeffULL, 0xffffffffffffffdfULL, 0x7ffffffffffffffbULL,
    0xffffffffffffefffULL, 0xfffffffffffffdffULL, 0xffffffffffffffbfULL, 0xfffffffffffffff7ULL,
    0xffffffffffffdfffULL, 0xfffffffffffffbffULL, 0xffffffffffffff7fULL, 0xffffffffffffffefULL,
    0xffffffffffffbfffULL, 0xfffffffffffff7ffULL, 0xfffffffffffffeffULL, 0xffffffffffffffdfULL,
    0xffffffffffff7fffULL, 0xffffffffffffefffULL, 0xfffffffffffffdffULL, 0xffffffffffffffbfULL,
    0xfffffffffffeffffULL, 0xffffffffffffdfffULL, 0xfffffffffffffbffULL, 0xffffffffffffff7fULL,
    0xfffffffffffdffffULL, 0xffffffffffffbfffULL, 0xfffffffffffff7ffULL, 0xfffffffffffffeffULL,
    0xfffffffffffbffffULL, 0xffffffffffff7fffULL, 0xffffffffffffefffULL, 0xfffffffffffffdffULL,
    0xfffffffffff7ffffULL, 0xfffffffffffeffffULL, 0xffffffffffffdfffULL, 0xfffffffffffffbffULL,
    0xffffffffffefffffULL, 0xfffffffffffdffffULL, 0xffffffffffffbfffULL, 0xfffffffffffff7ffULL,
    0xffffffffffdfffffULL, 0xfffffffffffbffffULL, 0xffffffffffff7fffULL, 0xffffffffffffefffULL,
    0xffffffffffbfffffULL, 0xfffffffffff7ffffULL, 0xfffffffffffeffffULL, 0xffffffffffffdfffULL,
    0xffffffffff7fffffULL, 0xffffffffffefffffULL, 0xfffffffffffdffffULL, 0xffffffffffffbfffULL,
    0xfffffffffeffffffULL, 0xffffffffffdfffffULL, 0xfffffffffffbffffULL, 0xffffffffffff7fffULL,
    0xfffffffffdffffffULL, 0xffffffffffbfffffULL, 0xfffffffffff7ffffULL, 0xfffffffffffeffffULL,
    0xfffffffffbffffffULL, 0xffffffffff7fffffULL, 0xffffffffffefffffULL, 0xfffffffffffdffffULL,
    0xfffffffff7ffffffULL, 0xfffffffffeffffffULL, 0xffffffffffdfffffULL, 0xfffffffffffbffffULL,
    0xffffffffefffffffULL, 0xfffffffffdffffffULL, 0xffffffffffbfffffULL, 0xfffffffffff7ffffULL,
    0xffffffffdfffffffULL, 0xfffffffffbffffffULL, 0xffffffffff7fffffULL, 0xffffffffffefffffULL,
    0xffffffffbfffffffULL, 0xfffffffff7ffffffULL, 0xfffffffffeffffffULL, 0xffffffffffdfffffULL,
    0xffffffff7fffffffULL, 0xffffffffefffffffULL, 0xfffffffffdffffffULL, 0xffffffffffbfffffULL,
    0xfffffffeffffffffULL, 0xffffffffdfffffffULL, 0xfffffffffbffffffULL, 0xffffffffff7fffffULL,
    0xfffffffdffffffffULL, 0xffffffffbfffffffULL, 0xfffffffff7ffffffULL, 0xfffffffffeffffffULL,
    0xfffffffbffffffffULL, 0xffffffff7fffffffULL, 0xffffffffefffffffULL, 0xfffffffffdffffffULL,
    0xfffffff7ffffffffULL, 0xfffffffeffffffffULL, 0xffffffffdfffffffULL, 0xfffffffffbffffffULL,
    0xffffffefffffffffULL, 0xfffffffdffffffffULL, 0xffffffffbfffffffULL, 0xfffffffff7ffffffULL,
    0xffffffdfffffffffULL, 0xfffffffbffffffffULL, 0xffffffff7fffffffULL, 0xffffffffefffffffULL,
    0xffffffbfffffffffULL, 0xfffffff7ffffffffULL, 0xfffffffeffffffffULL, 0xffffffffdfffffffULL,
    0xffffff7fffffffffULL, 0xffffffefffffffffULL, 0xfffffffdffffffffULL, 0xffffffffbfffffffULL,
    0xfffffeffffffffffULL, 0xffffffdfffffffffULL, 0xfffffffbffffffffULL, 0xffffffff7fffffffULL,
    0xfffffdffffffffffULL, 0xffffffbfffffffffULL, 0xfffffff7ffffffffULL, 0xfffffffeffffffffULL,
    0xfffffbffffffffffULL, 0xffffff7fffffffffULL, 0xffffffefffffffffULL, 0xfffffffdffffffffULL,
    0xfffff7ffffffffffULL, 0xfffffeffffffffffULL, 0xffffffdfffffffffULL, 0xfffffffbffffffffULL,
    0xffffefffffffffffULL, 0xfffffdffffffffffULL, 0xffffffbfffffffffULL, 0xfffffff7ffffffffULL,
    0xffffdfffffffffffULL, 0xfffffbffffffffffULL, 0xffffff7fffffffffULL, 0xffffffefffffffffULL,
    0xffffbfffffffffffULL, 0xfffff7ffffffffffULL, 0xfffffeffffffffffULL, 0xffffffdfffffffffULL,
    0xffff7fffffffffffULL, 0xffffefffffffffffULL, 0xfffffdffffffffffULL, 0xffffffbfffffffffULL,
    0xfffeffffffffffffULL, 0xffffdfffffffffffULL, 0xfffffbffffffffffULL, 0xffffff7fffffffffULL,
    0xfffdffffffffffffULL, 0xffffbfffffffffffULL, 0xfffff7ffffffffffULL, 0xfffffeffffffffffULL,
    0xfffbffffffffffffULL, 0xffff7fffffffffffULL, 0xffffefffffffffffULL, 0xfffffdffffffffffULL,
    0xfff7ffffffffffffULL, 0xfffeffffffffffffULL, 0xffffdfffffffffffULL, 0xfffffbffffffffffULL,
    0xffefffffffffffffULL, 0xfffdffffffffffffULL, 0xffffbfffffffffffULL, 0xfffff7ffffffffffULL,
    0xffdfffffffffffffULL, 0xfffbffffffffffffULL, 0xffff7fffffffffffULL, 0xffffefffffffffffULL,
    0xffbfffffffffffffULL, 0xfff7ffffffffffffULL, 0xfffeffffffffffffULL, 0xffffdfffffffffffULL,
    0xff7fffffffffffffULL, 0xffefffffffffffffULL, 0xfffdffffffffffffULL, 0xffffbfffffffffffULL,
    0xfeffffffffffffffULL, 0xffdfffffffffffffULL, 0xfffbffffffffffffULL, 0xffff7fffffffffffULL,
    0xfdffffffffffffffULL, 0xffbfffffffffffffULL, 0xfff7ffffffffffffULL, 0xfffeffffffffffffULL,
    0xfbffffffffffffffULL, 0xff7fffffffffffffULL, 0xffefffffffffffffULL, 0xfffdffffffffffffULL,
    0xf7ffffffffffffffULL, 0xfeffffffffffffffULL, 0xffdfffffffffffffULL, 0xfffbffffffffffffULL,
    0xefffffffffffffffULL, 0xfdffffffffffffffULL, 0xffbfffffffffffffULL, 0xfff7ffffffffffffULL };

uint64 _67_MASKS_256[67 * 4] = {
    0xfffffffffffffffeULL, 0xfffffffffffffff7ULL, 0xffffffffffffffbfULL, 0xfffffffffffffdffULL,
    0xfffffffffffffffdULL, 0xffffffffffffffefULL, 0xffffffffffffff7fULL, 0xfffffffffffffbffULL,
    0xfffffffffffffffbULL, 0xffffffffffffffdfULL, 0xfffffffffffffeffULL, 0xfffffffffffff7ffULL,
    0xfffffffffffffff7ULL, 0xffffffffffffffbfULL, 0xfffffffffffffdffULL, 0xffffffffffffefffULL,
    0xffffffffffffffefULL, 0xffffffffffffff7fULL, 0xfffffffffffffbffULL, 0xffffffffffffdfffULL,
    0xffffffffffffffdfULL, 0xfffffffffffffeffULL, 0xfffffffffffff7ffULL, 0xffffffffffffbfffULL,
    0xffffffffffffffbfULL, 0xfffffffffffffdffULL, 0xffffffffffffefffULL, 0xffffffffffff7fffULL,
    0xffffffffffffff7fULL, 0xfffffffffffffbffULL, 0xffffffffffffdfffULL, 0xfffffffffffeffffULL,
    0xfffffffffffffeffULL, 0xfffffffffffff7ffULL, 0xffffffffffffbfffULL, 0xfffffffffffdffffULL,
    0xfffffffffffffdffULL, 0xffffffffffffefffULL, 0xffffffffffff7fffULL, 0xfffffffffffbffffULL,
    0xfffffffffffffbffULL, 0xffffffffffffdfffULL, 0xfffffffffffeffffULL, 0xfffffffffff7ffffULL,
    0xfffffffffffff7ffULL, 0xffffffffffffbfffULL, 0xfffffffffffdffffULL, 0xffffffffffefffffULL,
    0xffffffffffffefffULL, 0xffffffffffff7fffULL, 0xfffffffffffbffffULL, 0xffffffffffdfffffULL,
    0xffffffffffffdfffULL, 0xfffffffffffeffffULL, 0xfffffffffff7ffffULL, 0xffffffffffbfffffULL,
    0xffffffffffffbfffULL, 0xfffffffffffdffffULL, 0xffffffffffefffffULL, 0xffffffffff7fffffULL,
    0xffffffffffff7fffULL, 0xfffffffffffbffffULL, 0xffffffffffdfffffULL, 0xfffffffffeffffffULL,
    0xfffffffffffeffffULL, 0xfffffffffff7ffffULL, 0xffffffffffbfffffULL, 0xfffffffffdffffffULL,
    0xfffffffffffdffffULL, 0xffffffffffefffffULL, 0xffffffffff7fffffULL, 0xfffffffffbffffffULL,
    0xfffffffffffbffffULL, 0xffffffffffdfffffULL, 0xfffffffffeffffffULL, 0xfffffffff7ffffffULL,
    0xfffffffffff7ffffULL, 0xffffffffffbfffffULL, 0xfffffffffdffffffULL, 0xffffffffefffffffULL,
    0xffffffffffefffffULL, 0xffffffffff7fffffULL, 0xfffffffffbffffffULL, 0xffffffffdfffffffULL,
    0xffffffffffdfffffULL, 0xfffffffffeffffffULL, 0xfffffffff7ffffffULL, 0xffffffffbfffffffULL,
    0xffffffffffbfffffULL, 0xfffffffffdffffffULL, 0xffffffffefffffffULL, 0xffffffff7fffffffULL,
    0xffffffffff7fffffULL, 0xfffffffffbffffffULL, 0xffffffffdfffffffULL, 0xfffffffeffffffffULL,
    0xfffffffffeffffffULL, 0xfffffffff7ffffffULL, 0xffffffffbfffffffULL, 0xfffffffdffffffffULL,
    0xfffffffffdffffffULL, 0xffffffffefffffffULL, 0xffffffff7fffffffULL, 0xfffffffbffffffffULL,
    0xfffffffffbffffffULL, 0xffffffffdfffffffULL, 0xfffffffeffffffffULL, 0xfffffff7ffffffffULL,
    0xfffffffff7ffffffULL, 0xffffffffbfffffffULL, 0xfffffffdffffffffULL, 0xffffffefffffffffULL,
    0xffffffffefffffffULL, 0xffffffff7fffffffULL, 0xfffffffbffffffffULL, 0xffffffdfffffffffULL,
    0xffffffffdfffffffULL, 0xfffffffeffffffffULL, 0xfffffff7ffffffffULL, 0xffffffbfffffffffULL,
    0xffffffffbfffffffULL, 0xfffffffdffffffffULL, 0xffffffefffffffffULL, 0xffffff7fffffffffULL,
    0xffffffff7fffffffULL, 0xfffffffbffffffffULL, 0xffffffdfffffffffULL, 0xfffffeffffffffffULL,
    0xfffffffeffffffffULL, 0xfffffff7ffffffffULL, 0xffffffbfffffffffULL, 0xfffffdffffffffffULL,
    0xfffffffdffffffffULL, 0xffffffefffffffffULL, 0xffffff7fffffffffULL, 0xfffffbffffffffffULL,
    0xfffffffbffffffffULL, 0xffffffdfffffffffULL, 0xfffffeffffffffffULL, 0xfffff7ffffffffffULL,
    0xfffffff7ffffffffULL, 0xffffffbfffffffffULL, 0xfffffdffffffffffULL, 0xffffefffffffffffULL,
    0xffffffefffffffffULL, 0xffffff7fffffffffULL, 0xfffffbffffffffffULL, 0xffffdfffffffffffULL,
    0xffffffdfffffffffULL, 0xfffffeffffffffffULL, 0xfffff7ffffffffffULL, 0xffffbfffffffffffULL,
    0xffffffbfffffffffULL, 0xfffffdffffffffffULL, 0xffffefffffffffffULL, 0xffff7fffffffffffULL,
    0xffffff7fffffffffULL, 0xfffffbffffffffffULL, 0xffffdfffffffffffULL, 0xfffeffffffffffffULL,
    0xfffffeffffffffffULL, 0xfffff7ffffffffffULL, 0xffffbfffffffffffULL, 0xfffdffffffffffffULL,
    0xfffffdffffffffffULL, 0xffffefffffffffffULL, 0xffff7fffffffffffULL, 0xfffbffffffffffffULL,
    0xfffffbffffffffffULL, 0xffffdfffffffffffULL, 0xfffeffffffffffffULL, 0xfff7ffffffffffffULL,
    0xfffff7ffffffffffULL, 0xffffbfffffffffffULL, 0xfffdffffffffffffULL, 0xffefffffffffffffULL,
    0xffffefffffffffffULL, 0xffff7fffffffffffULL, 0xfffbffffffffffffULL, 0xffdfffffffffffffULL,
    0xffffdfffffffffffULL, 0xfffeffffffffffffULL, 0xfff7ffffffffffffULL, 0xffbfffffffffffffULL,
    0xffffbfffffffffffULL, 0xfffdffffffffffffULL, 0xffefffffffffffffULL, 0xff7fffffffffffffULL,
    0xffff7fffffffffffULL, 0xfffbffffffffffffULL, 0xffdfffffffffffffULL, 0xfeffffffffffffffULL,
    0xfffeffffffffffffULL, 0xfff7ffffffffffffULL, 0xffbfffffffffffffULL, 0xfdffffffffffffffULL,
    0xfffdffffffffffffULL, 0xffefffffffffffffULL, 0xff7fffffffffffffULL, 0xfbffffffffffffffULL,
    0xfffbffffffffffffULL, 0xffdfffffffffffffULL, 0xfeffffffffffffffULL, 0xf7ffffffffffffffULL,
    0xfff7ffffffffffffULL, 0xffbfffffffffffffULL, 0xfdffffffffffffffULL, 0xefffffffffffffffULL,
    0xffefffffffffffffULL, 0xff7fffffffffffffULL, 0xfbffffffffffffffULL, 0xdfffffffffffffffULL,
    0xffdfffffffffffffULL, 0xfeffffffffffffffULL, 0xf7ffffffffffffffULL, 0xbfffffffffffffffULL,
    0xffbfffffffffffffULL, 0xfdffffffffffffffULL, 0xefffffffffffffffULL, 0x7fffffffffffffffULL,
    0xff7fffffffffffffULL, 0xfbffffffffffffffULL, 0xdfffffffffffffffULL, 0xffffffffffffffffULL,
    0xfeffffffffffffffULL, 0xf7ffffffffffffffULL, 0xbfffffffffffffffULL, 0xffffffffffffffffULL,
    0xfdffffffffffffffULL, 0xefffffffffffffffULL, 0x7fffffffffffffffULL, 0xffffffffffffffffULL,
    0xfbffffffffffffffULL, 0xdfffffffffffffffULL, 0xffffffffffffffffULL, 0xfffffffffffffffeULL,
    0xf7ffffffffffffffULL, 0xbfffffffffffffffULL, 0xffffffffffffffffULL, 0xfffffffffffffffdULL,
    0xefffffffffffffffULL, 0x7fffffffffffffffULL, 0xffffffffffffffffULL, 0xfffffffffffffffbULL,
    0xdfffffffffffffffULL, 0xffffffffffffffffULL, 0xfffffffffffffffeULL, 0xfffffffffffffff7ULL,
    0xbfffffffffffffffULL, 0xffffffffffffffffULL, 0xfffffffffffffffdULL, 0xffffffffffffffefULL,
    0x7fffffffffffffffULL, 0xffffffffffffffffULL, 0xfffffffffffffffbULL, 0xffffffffffffffdfULL,
    0xffffffffffffffffULL, 0xfffffffffffffffeULL, 0xfffffffffffffff7ULL, 0xffffffffffffffbfULL,
    0xffffffffffffffffULL, 0xfffffffffffffffdULL, 0xffffffffffffffefULL, 0xffffffffffffff7fULL,
    0xffffffffffffffffULL, 0xfffffffffffffffbULL, 0xffffffffffffffdfULL, 0xfffffffffffffeffULL };

uint64 _71_MASKS_256[71 * 4] = {
    0xfffffffffffffffeULL, 0xffffffffffffff7fULL, 0xffffffffffffbfffULL, 0xffffffffffdfffffULL,
    0xfffffffffffffffdULL, 0xfffffffffffffeffULL, 0xffffffffffff7fffULL, 0xffffffffffbfffffULL,
    0xfffffffffffffffbULL, 0xfffffffffffffdffULL, 0xfffffffffffeffffULL, 0xffffffffff7fffffULL,
    0xfffffffffffffff7ULL, 0xfffffffffffffbffULL, 0xfffffffffffdffffULL, 0xfffffffffeffffffULL,
    0xffffffffffffffefULL, 0xfffffffffffff7ffULL, 0xfffffffffffbffffULL, 0xfffffffffdffffffULL,
    0xffffffffffffffdfULL, 0xffffffffffffefffULL, 0xfffffffffff7ffffULL, 0xfffffffffbffffffULL,
    0xffffffffffffffbfULL, 0xffffffffffffdfffULL, 0xffffffffffefffffULL, 0xfffffffff7ffffffULL,
    0xffffffffffffff7fULL, 0xffffffffffffbfffULL, 0xffffffffffdfffffULL, 0xffffffffefffffffULL,
    0xfffffffffffffeffULL, 0xffffffffffff7fffULL, 0xffffffffffbfffffULL, 0xffffffffdfffffffULL,
    0xfffffffffffffdffULL, 0xfffffffffffeffffULL, 0xffffffffff7fffffULL, 0xffffffffbfffffffULL,
    0xfffffffffffffbffULL, 0xfffffffffffdffffULL, 0xfffffffffeffffffULL, 0xffffffff7fffffffULL,
    0xfffffffffffff7ffULL, 0xfffffffffffbffffULL, 0xfffffffffdffffffULL, 0xfffffffeffffffffULL,
    0xffffffffffffefffULL, 0xfffffffffff7ffffULL, 0xfffffffffbffffffULL, 0xfffffffdffffffffULL,
    0xffffffffffffdfffULL, 0xffffffffffefffffULL, 0xfffffffff7ffffffULL, 0xfffffffbffffffffULL,
    0xffffffffffffbfffULL, 0xffffffffffdfffffULL, 0xffffffffefffffffULL, 0xfffffff7ffffffffULL,
    0xffffffffffff7fffULL, 0xffffffffffbfffffULL, 0xffffffffdfffffffULL, 0xffffffefffffffffULL,
    0xfffffffffffeffffULL, 0xffffffffff7fffffULL, 0xffffffffbfffffffULL, 0xffffffdfffffffffULL,
    0xfffffffffffdffffULL, 0xfffffffffeffffffULL, 0xffffffff7fffffffULL, 0xffffffbfffffffffULL,
    0xfffffffffffbffffULL, 0xfffffffffdffffffULL, 0xfffffffeffffffffULL, 0xffffff7fffffffffULL,
    0xfffffffffff7ffffULL, 0xfffffffffbffffffULL, 0xfffffffdffffffffULL, 0xfffffeffffffffffULL,
    0xffffffffffefffffULL, 0xfffffffff7ffffffULL, 0xfffffffbffffffffULL, 0xfffffdffffffffffULL,
    0xffffffffffdfffffULL, 0xffffffffefffffffULL, 0xfffffff7ffffffffULL, 0xfffffbffffffffffULL,
    0xffffffffffbfffffULL, 0xffffffffdfffffffULL, 0xffffffefffffffffULL, 0xfffff7ffffffffffULL,
    0xffffffffff7fffffULL, 0xffffffffbfffffffULL, 0xffffffdfffffffffULL, 0xffffefffffffffffULL,
    0xfffffffffeffffffULL, 0xffffffff7fffffffULL, 0xffffffbfffffffffULL, 0xffffdfffffffffffULL,
    0xfffffffffdffffffULL, 0xfffffffeffffffffULL, 0xffffff7fffffffffULL, 0xffffbfffffffffffULL,
    0xfffffffffbffffffULL, 0xfffffffdffffffffULL, 0xfffffeffffffffffULL, 0xffff7fffffffffffULL,
    0xfffffffff7ffffffULL, 0xfffffffbffffffffULL, 0xfffffdffffffffffULL, 0xfffeffffffffffffULL,
    0xffffffffefffffffULL, 0xfffffff7ffffffffULL, 0xfffffbffffffffffULL, 0xfffdffffffffffffULL,
    0xffffffffdfffffffULL, 0xffffffefffffffffULL, 0xfffff7ffffffffffULL, 0xfffbffffffffffffULL,
    0xffffffffbfffffffULL, 0xffffffdfffffffffULL, 0xffffefffffffffffULL, 0xfff7ffffffffffffULL,
    0xffffffff7fffffffULL, 0xffffffbfffffffffULL, 0xffffdfffffffffffULL, 0xffefffffffffffffULL,
    0xfffffffeffffffffULL, 0xffffff7fffffffffULL, 0xffffbfffffffffffULL, 0xffdfffffffffffffULL,
    0xfffffffdffffffffULL, 0xfffffeffffffffffULL, 0xffff7fffffffffffULL, 0xffbfffffffffffffULL,
    0xfffffffbffffffffULL, 0xfffffdffffffffffULL, 0xfffeffffffffffffULL, 0xff7fffffffffffffULL,
    0xfffffff7ffffffffULL, 0xfffffbffffffffffULL, 0xfffdffffffffffffULL, 0xfeffffffffffffffULL,
    0xffffffefffffffffULL, 0xfffff7ffffffffffULL, 0xfffbffffffffffffULL, 0xfdffffffffffffffULL,
    0xffffffdfffffffffULL, 0xffffefffffffffffULL, 0xfff7ffffffffffffULL, 0xfbffffffffffffffULL,
    0xffffffbfffffffffULL, 0xffffdfffffffffffULL, 0xffefffffffffffffULL, 0xf7ffffffffffffffULL,
    0xffffff7fffffffffULL, 0xffffbfffffffffffULL, 0xffdfffffffffffffULL, 0xefffffffffffffffULL,
    0xfffffeffffffffffULL, 0xffff7fffffffffffULL, 0xffbfffffffffffffULL, 0xdfffffffffffffffULL,
    0xfffffdffffffffffULL, 0xfffeffffffffffffULL, 0xff7fffffffffffffULL, 0xbfffffffffffffffULL,
    0xfffffbffffffffffULL, 0xfffdffffffffffffULL, 0xfeffffffffffffffULL, 0x7fffffffffffffffULL,
    0xfffff7ffffffffffULL, 0xfffbffffffffffffULL, 0xfdffffffffffffffULL, 0xffffffffffffffffULL,
    0xffffefffffffffffULL, 0xfff7ffffffffffffULL, 0xfbffffffffffffffULL, 0xffffffffffffffffULL,
    0xffffdfffffffffffULL, 0xffefffffffffffffULL, 0xf7ffffffffffffffULL, 0xffffffffffffffffULL,
    0xffffbfffffffffffULL, 0xffdfffffffffffffULL, 0xefffffffffffffffULL, 0xffffffffffffffffULL,
    0xffff7fffffffffffULL, 0xffbfffffffffffffULL, 0xdfffffffffffffffULL, 0xffffffffffffffffULL,
    0xfffeffffffffffffULL, 0xff7fffffffffffffULL, 0xbfffffffffffffffULL, 0xffffffffffffffffULL,
    0xfffdffffffffffffULL, 0xfeffffffffffffffULL, 0x7fffffffffffffffULL, 0xffffffffffffffffULL,
    0xfffbffffffffffffULL, 0xfdffffffffffffffULL, 0xffffffffffffffffULL, 0xfffffffffffffffeULL,
    0xfff7ffffffffffffULL, 0xfbffffffffffffffULL, 0xffffffffffffffffULL, 0xfffffffffffffffdULL,
    0xffefffffffffffffULL, 0xf7ffffffffffffffULL, 0xffffffffffffffffULL, 0xfffffffffffffffbULL,
    0xffdfffffffffffffULL, 0xefffffffffffffffULL, 0xffffffffffffffffULL, 0xfffffffffffffff7ULL,
    0xffbfffffffffffffULL, 0xdfffffffffffffffULL, 0xffffffffffffffffULL, 0xffffffffffffffefULL,
    0xff7fffffffffffffULL, 0xbfffffffffffffffULL, 0xffffffffffffffffULL, 0xffffffffffffffdfULL,
    0xfeffffffffffffffULL, 0x7fffffffffffffffULL, 0xffffffffffffffffULL, 0xffffffffffffffbfULL,
    0xfdffffffffffffffULL, 0xffffffffffffffffULL, 0xfffffffffffffffeULL, 0xffffffffffffff7fULL,
    0xfbffffffffffffffULL, 0xffffffffffffffffULL, 0xfffffffffffffffdULL, 0xfffffffffffffeffULL,
    0xf7ffffffffffffffULL, 0xffffffffffffffffULL, 0xfffffffffffffffbULL, 0xfffffffffffffdffULL,
    0xefffffffffffffffULL, 0xffffffffffffffffULL, 0xfffffffffffffff7ULL, 0xfffffffffffffbffULL,
    0xdfffffffffffffffULL, 0xffffffffffffffffULL, 0xffffffffffffffefULL, 0xfffffffffffff7ffULL,
    0xbfffffffffffffffULL, 0xffffffffffffffffULL, 0xffffffffffffffdfULL, 0xffffffffffffefffULL,
    0x7fffffffffffffffULL, 0xffffffffffffffffULL, 0xffffffffffffffbfULL, 0xffffffffffffdfffULL,
    0xffffffffffffffffULL, 0xfffffffffffffffeULL, 0xffffffffffffff7fULL, 0xffffffffffffbfffULL,
    0xffffffffffffffffULL, 0xfffffffffffffffdULL, 0xfffffffffffffeffULL, 0xffffffffffff7fffULL,
    0xffffffffffffffffULL, 0xfffffffffffffffbULL, 0xfffffffffffffdffULL, 0xfffffffffffeffffULL,
    0xffffffffffffffffULL, 0xfffffffffffffff7ULL, 0xfffffffffffffbffULL, 0xfffffffffffdffffULL,
    0xffffffffffffffffULL, 0xffffffffffffffefULL, 0xfffffffffffff7ffULL, 0xfffffffffffbffffULL,
    0xffffffffffffffffULL, 0xffffffffffffffdfULL, 0xffffffffffffefffULL, 0xfffffffffff7ffffULL,
    0xffffffffffffffffULL, 0xffffffffffffffbfULL, 0xffffffffffffdfffULL, 0xffffffffffefffffULL };

uint64 _73_MASKS_256[73 * 4] = {
    0xfffffffffffffffeULL, 0xfffffffffffffdffULL, 0xfffffffffffbffffULL, 0xfffffffff7ffffffULL,
    0xfffffffffffffffdULL, 0xfffffffffffffbffULL, 0xfffffffffff7ffffULL, 0xffffffffefffffffULL,
    0xfffffffffffffffbULL, 0xfffffffffffff7ffULL, 0xffffffffffefffffULL, 0xffffffffdfffffffULL,
    0xfffffffffffffff7ULL, 0xffffffffffffefffULL, 0xffffffffffdfffffULL, 0xffffffffbfffffffULL,
    0xffffffffffffffefULL, 0xffffffffffffdfffULL, 0xffffffffffbfffffULL, 0xffffffff7fffffffULL,
    0xffffffffffffffdfULL, 0xffffffffffffbfffULL, 0xffffffffff7fffffULL, 0xfffffffeffffffffULL,
    0xffffffffffffffbfULL, 0xffffffffffff7fffULL, 0xfffffffffeffffffULL, 0xfffffffdffffffffULL,
    0xffffffffffffff7fULL, 0xfffffffffffeffffULL, 0xfffffffffdffffffULL, 0xfffffffbffffffffULL,
    0xfffffffffffffeffULL, 0xfffffffffffdffffULL, 0xfffffffffbffffffULL, 0xfffffff7ffffffffULL,
    0xfffffffffffffdffULL, 0xfffffffffffbffffULL, 0xfffffffff7ffffffULL, 0xffffffefffffffffULL,
    0xfffffffffffffbffULL, 0xfffffffffff7ffffULL, 0xffffffffefffffffULL, 0xffffffdfffffffffULL,
    0xfffffffffffff7ffULL, 0xffffffffffefffffULL, 0xffffffffdfffffffULL, 0xffffffbfffffffffULL,
    0xffffffffffffefffULL, 0xffffffffffdfffffULL, 0xffffffffbfffffffULL, 0xffffff7fffffffffULL,
    0xffffffffffffdfffULL, 0xffffffffffbfffffULL, 0xffffffff7fffffffULL, 0xfffffeffffffffffULL,
    0xffffffffffffbfffULL, 0xffffffffff7fffffULL, 0xfffffffeffffffffULL, 0xfffffdffffffffffULL,
    0xffffffffffff7fffULL, 0xfffffffffeffffffULL, 0xfffffffdffffffffULL, 0xfffffbffffffffffULL,
    0xfffffffffffeffffULL, 0xfffffffffdffffffULL, 0xfffffffbffffffffULL, 0xfffff7ffffffffffULL,
    0xfffffffffffdffffULL, 0xfffffffffbffffffULL, 0xfffffff7ffffffffULL, 0xffffefffffffffffULL,
    0xfffffffffffbffffULL, 0xfffffffff7ffffffULL, 0xffffffefffffffffULL, 0xffffdfffffffffffULL,
    0xfffffffffff7ffffULL, 0xffffffffefffffffULL, 0xffffffdfffffffffULL, 0xffffbfffffffffffULL,
    0xffffffffffefffffULL, 0xffffffffdfffffffULL, 0xffffffbfffffffffULL, 0xffff7fffffffffffULL,
    0xffffffffffdfffffULL, 0xffffffffbfffffffULL, 0xffffff7fffffffffULL, 0xfffeffffffffffffULL,
    0xffffffffffbfffffULL, 0xffffffff7fffffffULL, 0xfffffeffffffffffULL, 0xfffdffffffffffffULL,
    0xffffffffff7fffffULL, 0xfffffffeffffffffULL, 0xfffffdffffffffffULL, 0xfffbffffffffffffULL,
    0xfffffffffeffffffULL, 0xfffffffdffffffffULL, 0xfffffbffffffffffULL, 0xfff7ffffffffffffULL,
    0xfffffffffdffffffULL, 0xfffffffbffffffffULL, 0xfffff7ffffffffffULL, 0xffefffffffffffffULL,
    0xfffffffffbffffffULL, 0xfffffff7ffffffffULL, 0xffffefffffffffffULL, 0xffdfffffffffffffULL,
    0xfffffffff7ffffffULL, 0xffffffefffffffffULL, 0xffffdfffffffffffULL, 0xffbfffffffffffffULL,
    0xffffffffefffffffULL, 0xffffffdfffffffffULL, 0xffffbfffffffffffULL, 0xff7fffffffffffffULL,
    0xffffffffdfffffffULL, 0xffffffbfffffffffULL, 0xffff7fffffffffffULL, 0xfeffffffffffffffULL,
    0xffffffffbfffffffULL, 0xffffff7fffffffffULL, 0xfffeffffffffffffULL, 0xfdffffffffffffffULL,
    0xffffffff7fffffffULL, 0xfffffeffffffffffULL, 0xfffdffffffffffffULL, 0xfbffffffffffffffULL,
    0xfffffffeffffffffULL, 0xfffffdffffffffffULL, 0xfffbffffffffffffULL, 0xf7ffffffffffffffULL,
    0xfffffffdffffffffULL, 0xfffffbffffffffffULL, 0xfff7ffffffffffffULL, 0xefffffffffffffffULL,
    0xfffffffbffffffffULL, 0xfffff7ffffffffffULL, 0xffefffffffffffffULL, 0xdfffffffffffffffULL,
    0xfffffff7ffffffffULL, 0xffffefffffffffffULL, 0xffdfffffffffffffULL, 0xbfffffffffffffffULL,
    0xffffffefffffffffULL, 0xffffdfffffffffffULL, 0xffbfffffffffffffULL, 0x7fffffffffffffffULL,
    0xffffffdfffffffffULL, 0xffffbfffffffffffULL, 0xff7fffffffffffffULL, 0xffffffffffffffffULL,
    0xffffffbfffffffffULL, 0xffff7fffffffffffULL, 0xfeffffffffffffffULL, 0xffffffffffffffffULL,
    0xffffff7fffffffffULL, 0xfffeffffffffffffULL, 0xfdffffffffffffffULL, 0xffffffffffffffffULL,
    0xfffffeffffffffffULL, 0xfffdffffffffffffULL, 0xfbffffffffffffffULL, 0xffffffffffffffffULL,
    0xfffffdffffffffffULL, 0xfffbffffffffffffULL, 0xf7ffffffffffffffULL, 0xffffffffffffffffULL,
    0xfffffbffffffffffULL, 0xfff7ffffffffffffULL, 0xefffffffffffffffULL, 0xffffffffffffffffULL,
    0xfffff7ffffffffffULL, 0xffefffffffffffffULL, 0xdfffffffffffffffULL, 0xffffffffffffffffULL,
    0xffffefffffffffffULL, 0xffdfffffffffffffULL, 0xbfffffffffffffffULL, 0xffffffffffffffffULL,
    0xffffdfffffffffffULL, 0xffbfffffffffffffULL, 0x7fffffffffffffffULL, 0xffffffffffffffffULL,
    0xffffbfffffffffffULL, 0xff7fffffffffffffULL, 0xffffffffffffffffULL, 0xfffffffffffffffeULL,
    0xffff7fffffffffffULL, 0xfeffffffffffffffULL, 0xffffffffffffffffULL, 0xfffffffffffffffdULL,
    0xfffeffffffffffffULL, 0xfdffffffffffffffULL, 0xffffffffffffffffULL, 0xfffffffffffffffbULL,
    0xfffdffffffffffffULL, 0xfbffffffffffffffULL, 0xffffffffffffffffULL, 0xfffffffffffffff7ULL,
    0xfffbffffffffffffULL, 0xf7ffffffffffffffULL, 0xffffffffffffffffULL, 0xffffffffffffffefULL,
    0xfff7ffffffffffffULL, 0xefffffffffffffffULL, 0xffffffffffffffffULL, 0xffffffffffffffdfULL,
    0xffefffffffffffffULL, 0xdfffffffffffffffULL, 0xffffffffffffffffULL, 0xffffffffffffffbfULL,
    0xffdfffffffffffffULL, 0xbfffffffffffffffULL, 0xffffffffffffffffULL, 0xffffffffffffff7fULL,
    0xffbfffffffffffffULL, 0x7fffffffffffffffULL, 0xffffffffffffffffULL, 0xfffffffffffffeffULL,
    0xff7fffffffffffffULL, 0xffffffffffffffffULL, 0xfffffffffffffffeULL, 0xfffffffffffffdffULL,
    0xfeffffffffffffffULL, 0xffffffffffffffffULL, 0xfffffffffffffffdULL, 0xfffffffffffffbffULL,
    0xfdffffffffffffffULL, 0xffffffffffffffffULL, 0xfffffffffffffffbULL, 0xfffffffffffff7ffULL,
    0xfbffffffffffffffULL, 0xffffffffffffffffULL, 0xfffffffffffffff7ULL, 0xffffffffffffefffULL,
    0xf7ffffffffffffffULL, 0xffffffffffffffffULL, 0xffffffffffffffefULL, 0xffffffffffffdfffULL,
    0xefffffffffffffffULL, 0xffffffffffffffffULL, 0xffffffffffffffdfULL, 0xffffffffffffbfffULL,
    0xdfffffffffffffffULL, 0xffffffffffffffffULL, 0xffffffffffffffbfULL, 0xffffffffffff7fffULL,
    0xbfffffffffffffffULL, 0xffffffffffffffffULL, 0xffffffffffffff7fULL, 0xfffffffffffeffffULL,
    0x7fffffffffffffffULL, 0xffffffffffffffffULL, 0xfffffffffffffeffULL, 0xfffffffffffdffffULL,
    0xffffffffffffffffULL, 0xfffffffffffffffeULL, 0xfffffffffffffdffULL, 0xfffffffffffbffffULL,
    0xffffffffffffffffULL, 0xfffffffffffffffdULL, 0xfffffffffffffbffULL, 0xfffffffffff7ffffULL,
    0xffffffffffffffffULL, 0xfffffffffffffffbULL, 0xfffffffffffff7ffULL, 0xffffffffffefffffULL,
    0xffffffffffffffffULL, 0xfffffffffffffff7ULL, 0xffffffffffffefffULL, 0xffffffffffdfffffULL,
    0xffffffffffffffffULL, 0xffffffffffffffefULL, 0xffffffffffffdfffULL, 0xffffffffffbfffffULL,
    0xffffffffffffffffULL, 0xffffffffffffffdfULL, 0xffffffffffffbfffULL, 0xffffffffff7fffffULL,
    0xffffffffffffffffULL, 0xffffffffffffffbfULL, 0xffffffffffff7fffULL, 0xfffffffffeffffffULL,
    0xffffffffffffffffULL, 0xffffffffffffff7fULL, 0xfffffffffffeffffULL, 0xfffffffffdffffffULL,
    0xffffffffffffffffULL, 0xfffffffffffffeffULL, 0xfffffffffffdffffULL, 0xfffffffffbffffffULL };

uint64 _79_MASKS_256[79 * 4] = {
    0xfffffffffffffffeULL, 0xffffffffffff7fffULL, 0xffffffffbfffffffULL, 0xffffdfffffffffffULL,
    0xfffffffffffffffdULL, 0xfffffffffffeffffULL, 0xffffffff7fffffffULL, 0xffffbfffffffffffULL,
    0xfffffffffffffffbULL, 0xfffffffffffdffffULL, 0xfffffffeffffffffULL, 0xffff7fffffffffffULL,
    0xfffffffffffffff7ULL, 0xfffffffffffbffffULL, 0xfffffffdffffffffULL, 0xfffeffffffffffffULL,
    0xffffffffffffffefULL, 0xfffffffffff7ffffULL, 0xfffffffbffffffffULL, 0xfffdffffffffffffULL,
    0xffffffffffffffdfULL, 0xffffffffffefffffULL, 0xfffffff7ffffffffULL, 0xfffbffffffffffffULL,
    0xffffffffffffffbfULL, 0xffffffffffdfffffULL, 0xffffffefffffffffULL, 0xfff7ffffffffffffULL,
    0xffffffffffffff7fULL, 0xffffffffffbfffffULL, 0xffffffdfffffffffULL, 0xffefffffffffffffULL,
    0xfffffffffffffeffULL, 0xffffffffff7fffffULL, 0xffffffbfffffffffULL, 0xffdfffffffffffffULL,
    0xfffffffffffffdffULL, 0xfffffffffeffffffULL, 0xffffff7fffffffffULL, 0xffbfffffffffffffULL,
    0xfffffffffffffbffULL, 0xfffffffffdffffffULL, 0xfffffeffffffffffULL, 0xff7fffffffffffffULL,
    0xfffffffffffff7ffULL, 0xfffffffffbffffffULL, 0xfffffdffffffffffULL, 0xfeffffffffffffffULL,
    0xffffffffffffefffULL, 0xfffffffff7ffffffULL, 0xfffffbffffffffffULL, 0xfdffffffffffffffULL,
    0xffffffffffffdfffULL, 0xffffffffefffffffULL, 0xfffff7ffffffffffULL, 0xfbffffffffffffffULL,
    0xffffffffffffbfffULL, 0xffffffffdfffffffULL, 0xffffefffffffffffULL, 0xf7ffffffffffffffULL,
    0xffffffffffff7fffULL, 0xffffffffbfffffffULL, 0xffffdfffffffffffULL, 0xefffffffffffffffULL,
    0xfffffffffffeffffULL, 0xffffffff7fffffffULL, 0xffffbfffffffffffULL, 0xdfffffffffffffffULL,
    0xfffffffffffdffffULL, 0xfffffffeffffffffULL, 0xffff7fffffffffffULL, 0xbfffffffffffffffULL,
    0xfffffffffffbffffULL, 0xfffffffdffffffffULL, 0xfffeffffffffffffULL, 0x7fffffffffffffffULL,
    0xfffffffffff7ffffULL, 0xfffffffbffffffffULL, 0xfffdffffffffffffULL, 0xffffffffffffffffULL,
    0xffffffffffefffffULL, 0xfffffff7ffffffffULL, 0xfffbffffffffffffULL, 0xffffffffffffffffULL,
    0xffffffffffdfffffULL, 0xffffffefffffffffULL, 0xfff7ffffffffffffULL, 0xffffffffffffffffULL,
    0xffffffffffbfffffULL, 0xffffffdfffffffffULL, 0xffefffffffffffffULL, 0xffffffffffffffffULL,
    0xffffffffff7fffffULL, 0xffffffbfffffffffULL, 0xffdfffffffffffffULL, 0xffffffffffffffffULL,
    0xfffffffffeffffffULL, 0xffffff7fffffffffULL, 0xffbfffffffffffffULL, 0xffffffffffffffffULL,
    0xfffffffffdffffffULL, 0xfffffeffffffffffULL, 0xff7fffffffffffffULL, 0xffffffffffffffffULL,
    0xfffffffffbffffffULL, 0xfffffdffffffffffULL, 0xfeffffffffffffffULL, 0xffffffffffffffffULL,
    0xfffffffff7ffffffULL, 0xfffffbffffffffffULL, 0xfdffffffffffffffULL, 0xffffffffffffffffULL,
    0xffffffffefffffffULL, 0xfffff7ffffffffffULL, 0xfbffffffffffffffULL, 0xffffffffffffffffULL,
    0xffffffffdfffffffULL, 0xffffefffffffffffULL, 0xf7ffffffffffffffULL, 0xffffffffffffffffULL,
    0xffffffffbfffffffULL, 0xffffdfffffffffffULL, 0xefffffffffffffffULL, 0xffffffffffffffffULL,
    0xffffffff7fffffffULL, 0xffffbfffffffffffULL, 0xdfffffffffffffffULL, 0xffffffffffffffffULL,
    0xfffffffeffffffffULL, 0xffff7fffffffffffULL, 0xbfffffffffffffffULL, 0xffffffffffffffffULL,
    0xfffffffdffffffffULL, 0xfffeffffffffffffULL, 0x7fffffffffffffffULL, 0xffffffffffffffffULL,
    0xfffffffbffffffffULL, 0xfffdffffffffffffULL, 0xffffffffffffffffULL, 0xfffffffffffffffeULL,
    0xfffffff7ffffffffULL, 0xfffbffffffffffffULL, 0xffffffffffffffffULL, 0xfffffffffffffffdULL,
    0xffffffefffffffffULL, 0xfff7ffffffffffffULL, 0xffffffffffffffffULL, 0xfffffffffffffffbULL,
    0xffffffdfffffffffULL, 0xffefffffffffffffULL, 0xffffffffffffffffULL, 0xfffffffffffffff7ULL,
    0xffffffbfffffffffULL, 0xffdfffffffffffffULL, 0xffffffffffffffffULL, 0xffffffffffffffefULL,
    0xffffff7fffffffffULL, 0xffbfffffffffffffULL, 0xffffffffffffffffULL, 0xffffffffffffffdfULL,
    0xfffffeffffffffffULL, 0xff7fffffffffffffULL, 0xffffffffffffffffULL, 0xffffffffffffffbfULL,
    0xfffffdffffffffffULL, 0xfeffffffffffffffULL, 0xffffffffffffffffULL, 0xffffffffffffff7fULL,
    0xfffffbffffffffffULL, 0xfdffffffffffffffULL, 0xffffffffffffffffULL, 0xfffffffffffffeffULL,
    0xfffff7ffffffffffULL, 0xfbffffffffffffffULL, 0xffffffffffffffffULL, 0xfffffffffffffdffULL,
    0xffffefffffffffffULL, 0xf7ffffffffffffffULL, 0xffffffffffffffffULL, 0xfffffffffffffbffULL,
    0xffffdfffffffffffULL, 0xefffffffffffffffULL, 0xffffffffffffffffULL, 0xfffffffffffff7ffULL,
    0xffffbfffffffffffULL, 0xdfffffffffffffffULL, 0xffffffffffffffffULL, 0xffffffffffffefffULL,
    0xffff7fffffffffffULL, 0xbfffffffffffffffULL, 0xffffffffffffffffULL, 0xffffffffffffdfffULL,
    0xfffeffffffffffffULL, 0x7fffffffffffffffULL, 0xffffffffffffffffULL, 0xffffffffffffbfffULL,
    0xfffdffffffffffffULL, 0xffffffffffffffffULL, 0xfffffffffffffffeULL, 0xffffffffffff7fffULL,
    0xfffbffffffffffffULL, 0xffffffffffffffffULL, 0xfffffffffffffffdULL, 0xfffffffffffeffffULL,
    0xfff7ffffffffffffULL, 0xffffffffffffffffULL, 0xfffffffffffffffbULL, 0xfffffffffffdffffULL,
    0xffefffffffffffffULL, 0xffffffffffffffffULL, 0xfffffffffffffff7ULL, 0xfffffffffffbffffULL,
    0xffdfffffffffffffULL, 0xffffffffffffffffULL, 0xffffffffffffffefULL, 0xfffffffffff7ffffULL,
    0xffbfffffffffffffULL, 0xffffffffffffffffULL, 0xffffffffffffffdfULL, 0xffffffffffefffffULL,
    0xff7fffffffffffffULL, 0xffffffffffffffffULL, 0xffffffffffffffbfULL, 0xffffffffffdfffffULL,
    0xfeffffffffffffffULL, 0xffffffffffffffffULL, 0xffffffffffffff7fULL, 0xffffffffffbfffffULL,
    0xfdffffffffffffffULL, 0xffffffffffffffffULL, 0xfffffffffffffeffULL, 0xffffffffff7fffffULL,
    0xfbffffffffffffffULL, 0xffffffffffffffffULL, 0xfffffffffffffdffULL, 0xfffffffffeffffffULL,
    0xf7ffffffffffffffULL, 0xffffffffffffffffULL, 0xfffffffffffffbffULL, 0xfffffffffdffffffULL,
    0xefffffffffffffffULL, 0xffffffffffffffffULL, 0xfffffffffffff7ffULL, 0xfffffffffbffffffULL,
    0xdfffffffffffffffULL, 0xffffffffffffffffULL, 0xffffffffffffefffULL, 0xfffffffff7ffffffULL,
    0xbfffffffffffffffULL, 0xffffffffffffffffULL, 0xffffffffffffdfffULL, 0xffffffffefffffffULL,
    0x7fffffffffffffffULL, 0xffffffffffffffffULL, 0xffffffffffffbfffULL, 0xffffffffdfffffffULL,
    0xffffffffffffffffULL, 0xfffffffffffffffeULL, 0xffffffffffff7fffULL, 0xffffffffbfffffffULL,
    0xffffffffffffffffULL, 0xfffffffffffffffdULL, 0xfffffffffffeffffULL, 0xffffffff7fffffffULL,
    0xffffffffffffffffULL, 0xfffffffffffffffbULL, 0xfffffffffffdffffULL, 0xfffffffeffffffffULL,
    0xffffffffffffffffULL, 0xfffffffffffffff7ULL, 0xfffffffffffbffffULL, 0xfffffffdffffffffULL,
    0xffffffffffffffffULL, 0xffffffffffffffefULL, 0xfffffffffff7ffffULL, 0xfffffffbffffffffULL,
    0xffffffffffffffffULL, 0xffffffffffffffdfULL, 0xffffffffffefffffULL, 0xfffffff7ffffffffULL,
    0xffffffffffffffffULL, 0xffffffffffffffbfULL, 0xffffffffffdfffffULL, 0xffffffefffffffffULL,
    0xffffffffffffffffULL, 0xffffffffffffff7fULL, 0xffffffffffbfffffULL, 0xffffffdfffffffffULL,
    0xffffffffffffffffULL, 0xfffffffffffffeffULL, 0xffffffffff7fffffULL, 0xffffffbfffffffffULL,
    0xffffffffffffffffULL, 0xfffffffffffffdffULL, 0xfffffffffeffffffULL, 0xffffff7fffffffffULL,
    0xffffffffffffffffULL, 0xfffffffffffffbffULL, 0xfffffffffdffffffULL, 0xfffffeffffffffffULL,
    0xffffffffffffffffULL, 0xfffffffffffff7ffULL, 0xfffffffffbffffffULL, 0xfffffdffffffffffULL,
    0xffffffffffffffffULL, 0xffffffffffffefffULL, 0xfffffffff7ffffffULL, 0xfffffbffffffffffULL,
    0xffffffffffffffffULL, 0xffffffffffffdfffULL, 0xffffffffefffffffULL, 0xfffff7ffffffffffULL,
    0xffffffffffffffffULL, 0xffffffffffffbfffULL, 0xffffffffdfffffffULL, 0xffffefffffffffffULL };

uint64 _83_MASKS_256[83 * 4] = {
    0xfffffffffffffffeULL, 0xfffffffffff7ffffULL, 0xffffffbfffffffffULL, 0xfdffffffffffffffULL,
    0xfffffffffffffffdULL, 0xffffffffffefffffULL, 0xffffff7fffffffffULL, 0xfbffffffffffffffULL,
    0xfffffffffffffffbULL, 0xffffffffffdfffffULL, 0xfffffeffffffffffULL, 0xf7ffffffffffffffULL,
    0xfffffffffffffff7ULL, 0xffffffffffbfffffULL, 0xfffffdffffffffffULL, 0xefffffffffffffffULL,
    0xffffffffffffffefULL, 0xffffffffff7fffffULL, 0xfffffbffffffffffULL, 0xdfffffffffffffffULL,
    0xffffffffffffffdfULL, 0xfffffffffeffffffULL, 0xfffff7ffffffffffULL, 0xbfffffffffffffffULL,
    0xffffffffffffffbfULL, 0xfffffffffdffffffULL, 0xffffefffffffffffULL, 0x7fffffffffffffffULL,
    0xffffffffffffff7fULL, 0xfffffffffbffffffULL, 0xffffdfffffffffffULL, 0xffffffffffffffffULL,
    0xfffffffffffffeffULL, 0xfffffffff7ffffffULL, 0xffffbfffffffffffULL, 0xffffffffffffffffULL,
    0xfffffffffffffdffULL, 0xffffffffefffffffULL, 0xffff7fffffffffffULL, 0xffffffffffffffffULL,
    0xfffffffffffffbffULL, 0xffffffffdfffffffULL, 0xfffeffffffffffffULL, 0xffffffffffffffffULL,
    0xfffffffffffff7ffULL, 0xffffffffbfffffffULL, 0xfffdffffffffffffULL, 0xffffffffffffffffULL,
    0xffffffffffffefffULL, 0xffffffff7fffffffULL, 0xfffbffffffffffffULL, 0xffffffffffffffffULL,
    0xffffffffffffdfffULL, 0xfffffffeffffffffULL, 0xfff7ffffffffffffULL, 0xffffffffffffffffULL,
    0xffffffffffffbfffULL, 0xfffffffdffffffffULL, 0xffefffffffffffffULL, 0xffffffffffffffffULL,
    0xffffffffffff7fffULL, 0xfffffffbffffffffULL, 0xffdfffffffffffffULL, 0xffffffffffffffffULL,
    0xfffffffffffeffffULL, 0xfffffff7ffffffffULL, 0xffbfffffffffffffULL, 0xffffffffffffffffULL,
    0xfffffffffffdffffULL, 0xffffffefffffffffULL, 0xff7fffffffffffffULL, 0xffffffffffffffffULL,
    0xfffffffffffbffffULL, 0xffffffdfffffffffULL, 0xfeffffffffffffffULL, 0xffffffffffffffffULL,
    0xfffffffffff7ffffULL, 0xffffffbfffffffffULL, 0xfdffffffffffffffULL, 0xffffffffffffffffULL,
    0xffffffffffefffffULL, 0xffffff7fffffffffULL, 0xfbffffffffffffffULL, 0xffffffffffffffffULL,
    0xffffffffffdfffffULL, 0xfffffeffffffffffULL, 0xf7ffffffffffffffULL, 0xffffffffffffffffULL,
    0xffffffffffbfffffULL, 0xfffffdffffffffffULL, 0xefffffffffffffffULL, 0xffffffffffffffffULL,
    0xffffffffff7fffffULL, 0xfffffbffffffffffULL, 0xdfffffffffffffffULL, 0xffffffffffffffffULL,
    0xfffffffffeffffffULL, 0xfffff7ffffffffffULL, 0xbfffffffffffffffULL, 0xffffffffffffffffULL,
    0xfffffffffdffffffULL, 0xffffefffffffffffULL, 0x7fffffffffffffffULL, 0xffffffffffffffffULL,
    0xfffffffffbffffffULL, 0xffffdfffffffffffULL, 0xffffffffffffffffULL, 0xfffffffffffffffeULL,
    0xfffffffff7ffffffULL, 0xffffbfffffffffffULL, 0xffffffffffffffffULL, 0xfffffffffffffffdULL,
    0xffffffffefffffffULL, 0xffff7fffffffffffULL, 0xffffffffffffffffULL, 0xfffffffffffffffbULL,
    0xffffffffdfffffffULL, 0xfffeffffffffffffULL, 0xffffffffffffffffULL, 0xfffffffffffffff7ULL,
    0xffffffffbfffffffULL, 0xfffdffffffffffffULL, 0xffffffffffffffffULL, 0xffffffffffffffefULL,
    0xffffffff7fffffffULL, 0xfffbffffffffffffULL, 0xffffffffffffffffULL, 0xffffffffffffffdfULL,
    0xfffffffeffffffffULL, 0xfff7ffffffffffffULL, 0xffffffffffffffffULL, 0xffffffffffffffbfULL,
    0xfffffffdffffffffULL, 0xffefffffffffffffULL, 0xffffffffffffffffULL, 0xffffffffffffff7fULL,
    0xfffffffbffffffffULL, 0xffdfffffffffffffULL, 0xffffffffffffffffULL, 0xfffffffffffffeffULL,
    0xfffffff7ffffffffULL, 0xffbfffffffffffffULL, 0xffffffffffffffffULL, 0xfffffffffffffdffULL,
    0xffffffefffffffffULL, 0xff7fffffffffffffULL, 0xffffffffffffffffULL, 0xfffffffffffffbffULL,
    0xffffffdfffffffffULL, 0xfeffffffffffffffULL, 0xffffffffffffffffULL, 0xfffffffffffff7ffULL,
    0xffffffbfffffffffULL, 0xfdffffffffffffffULL, 0xffffffffffffffffULL, 0xffffffffffffefffULL,
    0xffffff7fffffffffULL, 0xfbffffffffffffffULL, 0xffffffffffffffffULL, 0xffffffffffffdfffULL,
    0xfffffeffffffffffULL, 0xf7ffffffffffffffULL, 0xffffffffffffffffULL, 0xffffffffffffbfffULL,
    0xfffffdffffffffffULL, 0xefffffffffffffffULL, 0xffffffffffffffffULL, 0xffffffffffff7fffULL,
    0xfffffbffffffffffULL, 0xdfffffffffffffffULL, 0xffffffffffffffffULL, 0xfffffffffffeffffULL,
    0xfffff7ffffffffffULL, 0xbfffffffffffffffULL, 0xffffffffffffffffULL, 0xfffffffffffdffffULL,
    0xffffefffffffffffULL, 0x7fffffffffffffffULL, 0xffffffffffffffffULL, 0xfffffffffffbffffULL,
    0xffffdfffffffffffULL, 0xffffffffffffffffULL, 0xfffffffffffffffeULL, 0xfffffffffff7ffffULL,
    0xffffbfffffffffffULL, 0xffffffffffffffffULL, 0xfffffffffffffffdULL, 0xffffffffffefffffULL,
    0xffff7fffffffffffULL, 0xffffffffffffffffULL, 0xfffffffffffffffbULL, 0xffffffffffdfffffULL,
    0xfffeffffffffffffULL, 0xffffffffffffffffULL, 0xfffffffffffffff7ULL, 0xffffffffffbfffffULL,
    0xfffdffffffffffffULL, 0xffffffffffffffffULL, 0xffffffffffffffefULL, 0xffffffffff7fffffULL,
    0xfffbffffffffffffULL, 0xffffffffffffffffULL, 0xffffffffffffffdfULL, 0xfffffffffeffffffULL,
    0xfff7ffffffffffffULL, 0xffffffffffffffffULL, 0xffffffffffffffbfULL, 0xfffffffffdffffffULL,
    0xffefffffffffffffULL, 0xffffffffffffffffULL, 0xffffffffffffff7fULL, 0xfffffffffbffffffULL,
    0xffdfffffffffffffULL, 0xffffffffffffffffULL, 0xfffffffffffffeffULL, 0xfffffffff7ffffffULL,
    0xffbfffffffffffffULL, 0xffffffffffffffffULL, 0xfffffffffffffdffULL, 0xffffffffefffffffULL,
    0xff7fffffffffffffULL, 0xffffffffffffffffULL, 0xfffffffffffffbffULL, 0xffffffffdfffffffULL,
    0xfeffffffffffffffULL, 0xffffffffffffffffULL, 0xfffffffffffff7ffULL, 0xffffffffbfffffffULL,
    0xfdffffffffffffffULL, 0xffffffffffffffffULL, 0xffffffffffffefffULL, 0xffffffff7fffffffULL,
    0xfbffffffffffffffULL, 0xffffffffffffffffULL, 0xffffffffffffdfffULL, 0xfffffffeffffffffULL,
    0xf7ffffffffffffffULL, 0xffffffffffffffffULL, 0xffffffffffffbfffULL, 0xfffffffdffffffffULL,
    0xefffffffffffffffULL, 0xffffffffffffffffULL, 0xffffffffffff7fffULL, 0xfffffffbffffffffULL,
    0xdfffffffffffffffULL, 0xffffffffffffffffULL, 0xfffffffffffeffffULL, 0xfffffff7ffffffffULL,
    0xbfffffffffffffffULL, 0xffffffffffffffffULL, 0xfffffffffffdffffULL, 0xffffffefffffffffULL,
    0x7fffffffffffffffULL, 0xffffffffffffffffULL, 0xfffffffffffbffffULL, 0xffffffdfffffffffULL,
    0xffffffffffffffffULL, 0xfffffffffffffffeULL, 0xfffffffffff7ffffULL, 0xffffffbfffffffffULL,
    0xffffffffffffffffULL, 0xfffffffffffffffdULL, 0xffffffffffefffffULL, 0xffffff7fffffffffULL,
    0xffffffffffffffffULL, 0xfffffffffffffffbULL, 0xffffffffffdfffffULL, 0xfffffeffffffffffULL,
    0xffffffffffffffffULL, 0xfffffffffffffff7ULL, 0xffffffffffbfffffULL, 0xfffffdffffffffffULL,
    0xffffffffffffffffULL, 0xffffffffffffffefULL, 0xffffffffff7fffffULL, 0xfffffbffffffffffULL,
    0xffffffffffffffffULL, 0xffffffffffffffdfULL, 0xfffffffffeffffffULL, 0xfffff7ffffffffffULL,
    0xffffffffffffffffULL, 0xffffffffffffffbfULL, 0xfffffffffdffffffULL, 0xffffefffffffffffULL,
    0xffffffffffffffffULL, 0xffffffffffffff7fULL, 0xfffffffffbffffffULL, 0xffffdfffffffffffULL,
    0xffffffffffffffffULL, 0xfffffffffffffeffULL, 0xfffffffff7ffffffULL, 0xffffbfffffffffffULL,
    0xffffffffffffffffULL, 0xfffffffffffffdffULL, 0xffffffffefffffffULL, 0xffff7fffffffffffULL,
    0xffffffffffffffffULL, 0xfffffffffffffbffULL, 0xffffffffdfffffffULL, 0xfffeffffffffffffULL,
    0xffffffffffffffffULL, 0xfffffffffffff7ffULL, 0xffffffffbfffffffULL, 0xfffdffffffffffffULL,
    0xffffffffffffffffULL, 0xffffffffffffefffULL, 0xffffffff7fffffffULL, 0xfffbffffffffffffULL,
    0xffffffffffffffffULL, 0xffffffffffffdfffULL, 0xfffffffeffffffffULL, 0xfff7ffffffffffffULL,
    0xffffffffffffffffULL, 0xffffffffffffbfffULL, 0xfffffffdffffffffULL, 0xffefffffffffffffULL,
    0xffffffffffffffffULL, 0xffffffffffff7fffULL, 0xfffffffbffffffffULL, 0xffdfffffffffffffULL,
    0xffffffffffffffffULL, 0xfffffffffffeffffULL, 0xfffffff7ffffffffULL, 0xffbfffffffffffffULL,
    0xffffffffffffffffULL, 0xfffffffffffdffffULL, 0xffffffefffffffffULL, 0xff7fffffffffffffULL,
    0xffffffffffffffffULL, 0xfffffffffffbffffULL, 0xffffffdfffffffffULL, 0xfeffffffffffffffULL };

uint64 _89_MASKS_256[89 * 4] = {
    0xfffffffffffffffeULL, 0xfffffffffdffffffULL, 0xfffbffffffffffffULL, 0xffffffffffffffffULL,
    0xfffffffffffffffdULL, 0xfffffffffbffffffULL, 0xfff7ffffffffffffULL, 0xffffffffffffffffULL,
    0xfffffffffffffffbULL, 0xfffffffff7ffffffULL, 0xffefffffffffffffULL, 0xffffffffffffffffULL,
    0xfffffffffffffff7ULL, 0xffffffffefffffffULL, 0xffdfffffffffffffULL, 0xffffffffffffffffULL,
    0xffffffffffffffefULL, 0xffffffffdfffffffULL, 0xffbfffffffffffffULL, 0xffffffffffffffffULL,
    0xffffffffffffffdfULL, 0xffffffffbfffffffULL, 0xff7fffffffffffffULL, 0xffffffffffffffffULL,
    0xffffffffffffffbfULL, 0xffffffff7fffffffULL, 0xfeffffffffffffffULL, 0xffffffffffffffffULL,
    0xffffffffffffff7fULL, 0xfffffffeffffffffULL, 0xfdffffffffffffffULL, 0xffffffffffffffffULL,
    0xfffffffffffffeffULL, 0xfffffffdffffffffULL, 0xfbffffffffffffffULL, 0xffffffffffffffffULL,
    0xfffffffffffffdffULL, 0xfffffffbffffffffULL, 0xf7ffffffffffffffULL, 0xffffffffffffffffULL,
    0xfffffffffffffbffULL, 0xfffffff7ffffffffULL, 0xefffffffffffffffULL, 0xffffffffffffffffULL,
    0xfffffffffffff7ffULL, 0xffffffefffffffffULL, 0xdfffffffffffffffULL, 0xffffffffffffffffULL,
    0xffffffffffffefffULL, 0xffffffdfffffffffULL, 0xbfffffffffffffffULL, 0xffffffffffffffffULL,
    0xffffffffffffdfffULL, 0xffffffbfffffffffULL, 0x7fffffffffffffffULL, 0xffffffffffffffffULL,
    0xffffffffffffbfffULL, 0xffffff7fffffffffULL, 0xffffffffffffffffULL, 0xfffffffffffffffeULL,
    0xffffffffffff7fffULL, 0xfffffeffffffffffULL, 0xffffffffffffffffULL, 0xfffffffffffffffdULL,
    0xfffffffffffeffffULL, 0xfffffdffffffffffULL, 0xffffffffffffffffULL, 0xfffffffffffffffbULL,
    0xfffffffffffdffffULL, 0xfffffbffffffffffULL, 0xffffffffffffffffULL, 0xfffffffffffffff7ULL,
    0xfffffffffffbffffULL, 0xfffff7ffffffffffULL, 0xffffffffffffffffULL, 0xffffffffffffffefULL,
    0xfffffffffff7ffffULL, 0xffffefffffffffffULL, 0xffffffffffffffffULL, 0xffffffffffffffdfULL,
    0xffffffffffefffffULL, 0xffffdfffffffffffULL, 0xffffffffffffffffULL, 0xffffffffffffffbfULL,
    0xffffffffffdfffffULL, 0xffffbfffffffffffULL, 0xffffffffffffffffULL, 0xffffffffffffff7fULL,
    0xffffffffffbfffffULL, 0xffff7fffffffffffULL, 0xffffffffffffffffULL, 0xfffffffffffffeffULL,
    0xffffffffff7fffffULL, 0xfffeffffffffffffULL, 0xffffffffffffffffULL, 0xfffffffffffffdffULL,
    0xfffffffffeffffffULL, 0xfffdffffffffffffULL, 0xffffffffffffffffULL, 0xfffffffffffffbffULL,
    0xfffffffffdffffffULL, 0xfffbffffffffffffULL, 0xffffffffffffffffULL, 0xfffffffffffff7ffULL,
    0xfffffffffbffffffULL, 0xfff7ffffffffffffULL, 0xffffffffffffffffULL, 0xffffffffffffefffULL,
    0xfffffffff7ffffffULL, 0xffefffffffffffffULL, 0xffffffffffffffffULL, 0xffffffffffffdfffULL,
    0xffffffffefffffffULL, 0xffdfffffffffffffULL, 0xffffffffffffffffULL, 0xffffffffffffbfffULL,
    0xffffffffdfffffffULL, 0xffbfffffffffffffULL, 0xffffffffffffffffULL, 0xffffffffffff7fffULL,
    0xffffffffbfffffffULL, 0xff7fffffffffffffULL, 0xffffffffffffffffULL, 0xfffffffffffeffffULL,
    0xffffffff7fffffffULL, 0xfeffffffffffffffULL, 0xffffffffffffffffULL, 0xfffffffffffdffffULL,
    0xfffffffeffffffffULL, 0xfdffffffffffffffULL, 0xffffffffffffffffULL, 0xfffffffffffbffffULL,
    0xfffffffdffffffffULL, 0xfbffffffffffffffULL, 0xffffffffffffffffULL, 0xfffffffffff7ffffULL,
    0xfffffffbffffffffULL, 0xf7ffffffffffffffULL, 0xffffffffffffffffULL, 0xffffffffffefffffULL,
    0xfffffff7ffffffffULL, 0xefffffffffffffffULL, 0xffffffffffffffffULL, 0xffffffffffdfffffULL,
    0xffffffefffffffffULL, 0xdfffffffffffffffULL, 0xffffffffffffffffULL, 0xffffffffffbfffffULL,
    0xffffffdfffffffffULL, 0xbfffffffffffffffULL, 0xffffffffffffffffULL, 0xffffffffff7fffffULL,
    0xffffffbfffffffffULL, 0x7fffffffffffffffULL, 0xffffffffffffffffULL, 0xfffffffffeffffffULL,
    0xffffff7fffffffffULL, 0xffffffffffffffffULL, 0xfffffffffffffffeULL, 0xfffffffffdffffffULL,
    0xfffffeffffffffffULL, 0xffffffffffffffffULL, 0xfffffffffffffffdULL, 0xfffffffffbffffffULL,
    0xfffffdffffffffffULL, 0xffffffffffffffffULL, 0xfffffffffffffffbULL, 0xfffffffff7ffffffULL,
    0xfffffbffffffffffULL, 0xffffffffffffffffULL, 0xfffffffffffffff7ULL, 0xffffffffefffffffULL,
    0xfffff7ffffffffffULL, 0xffffffffffffffffULL, 0xffffffffffffffefULL, 0xffffffffdfffffffULL,
    0xffffefffffffffffULL, 0xffffffffffffffffULL, 0xffffffffffffffdfULL, 0xffffffffbfffffffULL,
    0xffffdfffffffffffULL, 0xffffffffffffffffULL, 0xffffffffffffffbfULL, 0xffffffff7fffffffULL,
    0xffffbfffffffffffULL, 0xffffffffffffffffULL, 0xffffffffffffff7fULL, 0xfffffffeffffffffULL,
    0xffff7fffffffffffULL, 0xffffffffffffffffULL, 0xfffffffffffffeffULL, 0xfffffffdffffffffULL,
    0xfffeffffffffffffULL, 0xffffffffffffffffULL, 0xfffffffffffffdffULL, 0xfffffffbffffffffULL,
    0xfffdffffffffffffULL, 0xffffffffffffffffULL, 0xfffffffffffffbffULL, 0xfffffff7ffffffffULL,
    0xfffbffffffffffffULL, 0xffffffffffffffffULL, 0xfffffffffffff7ffULL, 0xffffffefffffffffULL,
    0xfff7ffffffffffffULL, 0xffffffffffffffffULL, 0xffffffffffffefffULL, 0xffffffdfffffffffULL,
    0xffefffffffffffffULL, 0xffffffffffffffffULL, 0xffffffffffffdfffULL, 0xffffffbfffffffffULL,
    0xffdfffffffffffffULL, 0xffffffffffffffffULL, 0xffffffffffffbfffULL, 0xffffff7fffffffffULL,
    0xffbfffffffffffffULL, 0xffffffffffffffffULL, 0xffffffffffff7fffULL, 0xfffffeffffffffffULL,
    0xff7fffffffffffffULL, 0xffffffffffffffffULL, 0xfffffffffffeffffULL, 0xfffffdffffffffffULL,
    0xfeffffffffffffffULL, 0xffffffffffffffffULL, 0xfffffffffffdffffULL, 0xfffffbffffffffffULL,
    0xfdffffffffffffffULL, 0xffffffffffffffffULL, 0xfffffffffffbffffULL, 0xfffff7ffffffffffULL,
    0xfbffffffffffffffULL, 0xffffffffffffffffULL, 0xfffffffffff7ffffULL, 0xffffefffffffffffULL,
    0xf7ffffffffffffffULL, 0xffffffffffffffffULL, 0xffffffffffefffffULL, 0xffffdfffffffffffULL,
    0xefffffffffffffffULL, 0xffffffffffffffffULL, 0xffffffffffdfffffULL, 0xffffbfffffffffffULL,
    0xdfffffffffffffffULL, 0xffffffffffffffffULL, 0xffffffffffbfffffULL, 0xffff7fffffffffffULL,
    0xbfffffffffffffffULL, 0xffffffffffffffffULL, 0xffffffffff7fffffULL, 0xfffeffffffffffffULL,
    0x7fffffffffffffffULL, 0xffffffffffffffffULL, 0xfffffffffeffffffULL, 0xfffdffffffffffffULL,
    0xffffffffffffffffULL, 0xfffffffffffffffeULL, 0xfffffffffdffffffULL, 0xfffbffffffffffffULL,
    0xffffffffffffffffULL, 0xfffffffffffffffdULL, 0xfffffffffbffffffULL, 0xfff7ffffffffffffULL,
    0xffffffffffffffffULL, 0xfffffffffffffffbULL, 0xfffffffff7ffffffULL, 0xffefffffffffffffULL,
    0xffffffffffffffffULL, 0xfffffffffffffff7ULL, 0xffffffffefffffffULL, 0xffdfffffffffffffULL,
    0xffffffffffffffffULL, 0xffffffffffffffefULL, 0xffffffffdfffffffULL, 0xffbfffffffffffffULL,
    0xffffffffffffffffULL, 0xffffffffffffffdfULL, 0xffffffffbfffffffULL, 0xff7fffffffffffffULL,
    0xffffffffffffffffULL, 0xffffffffffffffbfULL, 0xffffffff7fffffffULL, 0xfeffffffffffffffULL,
    0xffffffffffffffffULL, 0xffffffffffffff7fULL, 0xfffffffeffffffffULL, 0xfdffffffffffffffULL,
    0xffffffffffffffffULL, 0xfffffffffffffeffULL, 0xfffffffdffffffffULL, 0xfbffffffffffffffULL,
    0xffffffffffffffffULL, 0xfffffffffffffdffULL, 0xfffffffbffffffffULL, 0xf7ffffffffffffffULL,
    0xffffffffffffffffULL, 0xfffffffffffffbffULL, 0xfffffff7ffffffffULL, 0xefffffffffffffffULL,
    0xffffffffffffffffULL, 0xfffffffffffff7ffULL, 0xffffffefffffffffULL, 0xdfffffffffffffffULL,
    0xffffffffffffffffULL, 0xffffffffffffefffULL, 0xffffffdfffffffffULL, 0xbfffffffffffffffULL,
    0xffffffffffffffffULL, 0xffffffffffffdfffULL, 0xffffffbfffffffffULL, 0x7fffffffffffffffULL,
    0xffffffffffffffffULL, 0xffffffffffffbfffULL, 0xffffff7fffffffffULL, 0xffffffffffffffffULL,
    0xffffffffffffffffULL, 0xffffffffffff7fffULL, 0xfffffeffffffffffULL, 0xffffffffffffffffULL,
    0xffffffffffffffffULL, 0xfffffffffffeffffULL, 0xfffffdffffffffffULL, 0xffffffffffffffffULL,
    0xffffffffffffffffULL, 0xfffffffffffdffffULL, 0xfffffbffffffffffULL, 0xffffffffffffffffULL,
    0xffffffffffffffffULL, 0xfffffffffffbffffULL, 0xfffff7ffffffffffULL, 0xffffffffffffffffULL,
    0xffffffffffffffffULL, 0xfffffffffff7ffffULL, 0xffffefffffffffffULL, 0xffffffffffffffffULL,
    0xffffffffffffffffULL, 0xffffffffffefffffULL, 0xffffdfffffffffffULL, 0xffffffffffffffffULL,
    0xffffffffffffffffULL, 0xffffffffffdfffffULL, 0xffffbfffffffffffULL, 0xffffffffffffffffULL,
    0xffffffffffffffffULL, 0xffffffffffbfffffULL, 0xffff7fffffffffffULL, 0xffffffffffffffffULL,
    0xffffffffffffffffULL, 0xffffffffff7fffffULL, 0xfffeffffffffffffULL, 0xffffffffffffffffULL,
    0xffffffffffffffffULL, 0xfffffffffeffffffULL, 0xfffdffffffffffffULL, 0xffffffffffffffffULL };


uint64 _5_MASKS[5] = {
    0xef7bdef7bdef7bdeULL,
    0xdef7bdef7bdef7bdULL,
    0xbdef7bdef7bdef7bULL,
    0x7bdef7bdef7bdef7ULL,
    0xf7bdef7bdef7bdefULL };

uint64 _7_MASKS[7] = {
    0x7efdfbf7efdfbf7eULL,
    0xfdfbf7efdfbf7efdULL,
    0xfbf7efdfbf7efdfbULL,
    0xf7efdfbf7efdfbf7ULL,
    0xefdfbf7efdfbf7efULL,
    0xdfbf7efdfbf7efdfULL,
    0xbf7efdfbf7efdfbfULL };

uint64 _11_MASKS[11] = {
    0xff7feffdffbff7feULL,
    0xfeffdffbff7feffdULL,
    0xfdffbff7feffdffbULL,
    0xfbff7feffdffbff7ULL,
    0xf7feffdffbff7fefULL,
    0xeffdffbff7feffdfULL,
    0xdffbff7feffdffbfULL,
    0xbff7feffdffbff7fULL,
    0x7feffdffbff7feffULL,
    0xffdffbff7feffdffULL,
    0xffbff7feffdffbffULL };

uint64 _13_MASKS[13] = {
    0xffefff7ffbffdffeULL,
    0xffdffefff7ffbffdULL,
    0xffbffdffefff7ffbULL,
    0xff7ffbffdffefff7ULL,
    0xfefff7ffbffdffefULL,
    0xfdffefff7ffbffdfULL,
    0xfbffdffefff7ffbfULL,
    0xf7ffbffdffefff7fULL,
    0xefff7ffbffdffeffULL,
    0xdffefff7ffbffdffULL,
    0xbffdffefff7ffbffULL,
    0x7ffbffdffefff7ffULL,
    0xfff7ffbffdffefffULL };

uint64 _17_MASKS[17] = {
    0xfff7fffbfffdfffeULL,
    0xffeffff7fffbfffdULL,
    0xffdfffeffff7fffbULL,
    0xffbfffdfffeffff7ULL,
    0xff7fffbfffdfffefULL,
    0xfeffff7fffbfffdfULL,
    0xfdfffeffff7fffbfULL,
    0xfbfffdfffeffff7fULL,
    0xf7fffbfffdfffeffULL,
    0xeffff7fffbfffdffULL,
    0xdfffeffff7fffbffULL,
    0xbfffdfffeffff7ffULL,
    0x7fffbfffdfffefffULL,
    0xffff7fffbfffdfffULL,
    0xfffeffff7fffbfffULL,
    0xfffdfffeffff7fffULL,
    0xfffbfffdfffeffffULL };

uint64 _19_MASKS[19] = {
    0xfdffffbffff7fffeULL,
    0xfbffff7fffeffffdULL,
    0xf7fffeffffdffffbULL,
    0xeffffdffffbffff7ULL,
    0xdffffbffff7fffefULL,
    0xbffff7fffeffffdfULL,
    0x7fffeffffdffffbfULL,
    0xffffdffffbffff7fULL,
    0xffffbffff7fffeffULL,
    0xffff7fffeffffdffULL,
    0xfffeffffdffffbffULL,
    0xfffdffffbffff7ffULL,
    0xfffbffff7fffefffULL,
    0xfff7fffeffffdfffULL,
    0xffeffffdffffbfffULL,
    0xffdffffbffff7fffULL,
    0xffbffff7fffeffffULL,
    0xff7fffeffffdffffULL,
    0xfeffffdffffbffffULL };

uint64 _23_MASKS[23] = {
    0xffffbfffff7ffffeULL,
    0xffff7ffffefffffdULL,
    0xfffefffffdfffffbULL,
    0xfffdfffffbfffff7ULL,
    0xfffbfffff7ffffefULL,
    0xfff7ffffefffffdfULL,
    0xffefffffdfffffbfULL,
    0xffdfffffbfffff7fULL,
    0xffbfffff7ffffeffULL,
    0xff7ffffefffffdffULL,
    0xfefffffdfffffbffULL,
    0xfdfffffbfffff7ffULL,
    0xfbfffff7ffffefffULL,
    0xf7ffffefffffdfffULL,
    0xefffffdfffffbfffULL,
    0xdfffffbfffff7fffULL,
    0xbfffff7ffffeffffULL,
    0x7ffffefffffdffffULL,
    0xfffffdfffffbffffULL,
    0xfffffbfffff7ffffULL,
    0xfffff7ffffefffffULL,
    0xffffefffffdfffffULL,
    0xffffdfffffbfffffULL };

uint64 _29_MASKS[29] = {
    0xfbffffffdffffffeULL,
    0xf7ffffffbffffffdULL,
    0xefffffff7ffffffbULL,
    0xdffffffefffffff7ULL,
    0xbffffffdffffffefULL,
    0x7ffffffbffffffdfULL,
    0xfffffff7ffffffbfULL,
    0xffffffefffffff7fULL,
    0xffffffdffffffeffULL,
    0xffffffbffffffdffULL,
    0xffffff7ffffffbffULL,
    0xfffffefffffff7ffULL,
    0xfffffdffffffefffULL,
    0xfffffbffffffdfffULL,
    0xfffff7ffffffbfffULL,
    0xffffefffffff7fffULL,
    0xffffdffffffeffffULL,
    0xffffbffffffdffffULL,
    0xffff7ffffffbffffULL,
    0xfffefffffff7ffffULL,
    0xfffdffffffefffffULL,
    0xfffbffffffdfffffULL,
    0xfff7ffffffbfffffULL,
    0xffefffffff7fffffULL,
    0xffdffffffeffffffULL,
    0xffbffffffdffffffULL,
    0xff7ffffffbffffffULL,
    0xfefffffff7ffffffULL,
    0xfdffffffefffffffULL };

#ifdef USE_AVX2
void pre_sieve_avx2(soe_dynamicdata_t *ddata, soe_staticdata_t *sdata, uint8 *flagblock)
{
    uint64 startprime = sdata->startprime;
    uint64 k;
    int i, j;
    int mask_step, mask_step2, mask_step3, mask_step4;
    int mask_num, mask_num2, mask_num3, mask_num4;
    uint64 *flagblock64;

    // flagblock is always a multiple of 8 bytes
    flagblock64 = (uint64 *)flagblock;

    // do the smallest primes in predetermined batches.
    // first handle the different startprime cases, going up through prime 13.
    if (startprime == 2)
    {
        for (k = 0, mask_step = _256_MOD_P[0], mask_num = ddata->offsets[2],
            mask_step2 = _256_MOD_P[1], mask_num2 = ddata->offsets[3];
            k < FLAGSIZE >> 8; k++)
        {
            __m256i vsieve = _mm256_load_si256((__m256i *)(&flagblock64[k * 4]));
            __m256i vmask1 = _mm256_load_si256((__m256i *)(&_5_MASKS_256[mask_num * 4]));
            __m256i vmask2 = _mm256_load_si256((__m256i *)(&_7_MASKS_256[mask_num2 * 4]));
            vsieve = _mm256_and_si256(vmask1, vsieve);
            vsieve = _mm256_and_si256(vmask2, vsieve);
            _mm256_store_si256((__m256i *)(&flagblock64[k * 4]), vsieve);

            mask_num -= mask_step;
            if (mask_num < 0) mask_num = 5 + mask_num;
            mask_num2 -= mask_step2;
            if (mask_num2 < 0) mask_num2 = 7 + mask_num2;
        }
        ddata->offsets[2] = (uint32)mask_num;
        ddata->offsets[3] = (uint32)mask_num2;

        for (k = 0, mask_step = _256_MOD_P[2], mask_num = ddata->offsets[4],
            mask_step2 = _256_MOD_P[3], mask_num2 = ddata->offsets[5];
            k < FLAGSIZE >> 8; k++)
        {
            __m256i vsieve = _mm256_load_si256((__m256i *)(&flagblock64[k * 4]));
            __m256i vmask1 = _mm256_load_si256((__m256i *)(&_11_MASKS_256[mask_num * 4]));
            __m256i vmask2 = _mm256_load_si256((__m256i *)(&_13_MASKS_256[mask_num2 * 4]));
            vsieve = _mm256_and_si256(vmask1, vsieve);
            vsieve = _mm256_and_si256(vmask2, vsieve);
            _mm256_store_si256((__m256i *)(&flagblock64[k * 4]), vsieve);

            mask_num -= mask_step;
            if (mask_num < 0) mask_num = 11 + mask_num;
            mask_num2 -= mask_step2;
            if (mask_num2 < 0) mask_num2 = 13 + mask_num2;
        }
        ddata->offsets[4] = (uint32)mask_num;
        ddata->offsets[5] = (uint32)mask_num2;



    }
    else if (startprime == 3)
    {
        for (k = 0, mask_step = _256_MOD_P[1], mask_num = ddata->offsets[3];
            k < FLAGSIZE >> 8; k++)
        {
            __m256i vsieve = _mm256_load_si256((__m256i *)(&flagblock64[k * 4]));
            __m256i vmask = _mm256_load_si256((__m256i *)(&_7_MASKS_256[mask_num * 4]));
            vsieve = _mm256_and_si256(vmask, vsieve);
            _mm256_store_si256((__m256i *)(&flagblock64[k * 4]), vsieve);
            mask_num -= mask_step;
            if (mask_num < 0) mask_num = 7 + mask_num;
        }
        ddata->offsets[3] = (uint32)mask_num;

        for (k = 0, mask_step = _256_MOD_P[2], mask_num = ddata->offsets[4],
            mask_step2 = _256_MOD_P[3], mask_num2 = ddata->offsets[5];
            k < FLAGSIZE >> 8; k++)
        {
            __m256i vsieve = _mm256_load_si256((__m256i *)(&flagblock64[k * 4]));
            __m256i vmask1 = _mm256_load_si256((__m256i *)(&_11_MASKS_256[mask_num * 4]));
            __m256i vmask2 = _mm256_load_si256((__m256i *)(&_13_MASKS_256[mask_num2 * 4]));
            vsieve = _mm256_and_si256(vmask1, vsieve);
            vsieve = _mm256_and_si256(vmask2, vsieve);
            _mm256_store_si256((__m256i *)(&flagblock64[k * 4]), vsieve);

            mask_num -= mask_step;
            if (mask_num < 0) mask_num = 11 + mask_num;
            mask_num2 -= mask_step2;
            if (mask_num2 < 0) mask_num2 = 13 + mask_num2;
        }
        ddata->offsets[4] = (uint32)mask_num;
        ddata->offsets[5] = (uint32)mask_num2;

    }
    else if (startprime == 4)
    {
        for (k = 0, mask_step = _256_MOD_P[2], mask_num = ddata->offsets[4],
            mask_step2 = _256_MOD_P[3], mask_num2 = ddata->offsets[5];
            k < FLAGSIZE >> 8; k++)
        {
            __m256i vsieve = _mm256_load_si256((__m256i *)(&flagblock64[k * 4]));
            __m256i vmask1 = _mm256_load_si256((__m256i *)(&_11_MASKS_256[mask_num * 4]));
            __m256i vmask2 = _mm256_load_si256((__m256i *)(&_13_MASKS_256[mask_num2 * 4]));
            vsieve = _mm256_and_si256(vmask1, vsieve);
            vsieve = _mm256_and_si256(vmask2, vsieve);
            _mm256_store_si256((__m256i *)(&flagblock64[k * 4]), vsieve);

            mask_num -= mask_step;
            if (mask_num < 0) mask_num = 11 + mask_num;
            mask_num2 -= mask_step2;
            if (mask_num2 < 0) mask_num2 = 13 + mask_num2;
        }

        ddata->offsets[4] = (uint32)mask_num;
        ddata->offsets[5] = (uint32)mask_num2;

    }
    else if (startprime == 5)
    {
        for (k = 0, mask_step = _256_MOD_P[3], mask_num = ddata->offsets[5];
            k < FLAGSIZE >> 8; k++)
        {
            __m256i vsieve = _mm256_load_si256((__m256i *)(&flagblock64[k * 4]));
            __m256i vmask = _mm256_load_si256((__m256i *)(&_13_MASKS_256[mask_num * 4]));
            vsieve = _mm256_and_si256(vmask, vsieve);
            _mm256_store_si256((__m256i *)(&flagblock64[k * 4]), vsieve);
            mask_num -= mask_step;
            if (mask_num < 0) mask_num = 13 + mask_num;
        }
        ddata->offsets[5] = (uint32)mask_num;
    }

    // then do some more presieving, starting at prime 17.
    for (k = 0, mask_step = _256_MOD_P[4], mask_num = ddata->offsets[6],
        mask_step2 = _256_MOD_P[5], mask_num2 = ddata->offsets[7];
        k < FLAGSIZE >> 8; k++)
    {
        __m256i vsieve = _mm256_load_si256((__m256i *)(&flagblock64[k * 4]));
        __m256i vmask1 = _mm256_load_si256((__m256i *)(&_17_MASKS_256[mask_num * 4]));
        __m256i vmask2 = _mm256_load_si256((__m256i *)(&_19_MASKS_256[mask_num2 * 4]));
        vsieve = _mm256_and_si256(vmask1, vsieve);
        vsieve = _mm256_and_si256(vmask2, vsieve);
        _mm256_store_si256((__m256i *)(&flagblock64[k * 4]), vsieve);

        mask_num -= mask_step;
        if (mask_num < 0) mask_num = 17 + mask_num;
        mask_num2 -= mask_step2;
        if (mask_num2 < 0) mask_num2 = 19 + mask_num2;
    }
    ddata->offsets[6] = (uint32)mask_num;
    ddata->offsets[7] = (uint32)mask_num2;

    {
        __m256i vmasknums = _mm256_load_si256((__m256i *)(ddata->offsets + 8));
        __m256i vmaskp = _mm256_setr_epi32(23, 29, 31, 37, 41, 43, 47, 53);
        __m256i vmaskp1 = _mm256_setr_epi32(22, 28, 30, 36, 40, 42, 46, 52);
        __m256i vzero = _mm256_set1_epi32(0);
        __m256i vmaskstep = _mm256_setr_epi32(3, 24, 8, 34, 10, 41, 21, 44);
        __m256i vmask, vmask2;
#ifdef __INTEL_COMPILER
        __declspec(align(64)) uint32 t[8];
#else
        uint32 *t = ddata->presieve_scratch;
#endif

        for (k = 0; k < FLAGSIZE >> 8; k++)
        {
            __m256i vsieve = _mm256_load_si256((__m256i *)(&flagblock64[k * 4]));
            _mm256_store_si256((__m256i *)t, vmasknums);

            vmask = _mm256_load_si256((__m256i *)(&_23_MASKS_256[t[0] * 4]));
            vsieve = _mm256_and_si256(vmask, vsieve);

            vmask = _mm256_load_si256((__m256i *)(&_29_MASKS_256[t[1] * 4]));
            vsieve = _mm256_and_si256(vmask, vsieve);

            vmask = _mm256_load_si256((__m256i *)(&_31_MASKS_256[t[2] * 4]));
            vsieve = _mm256_and_si256(vmask, vsieve);

            vmask = _mm256_load_si256((__m256i *)(&_37_MASKS_256[t[3] * 4]));
            vsieve = _mm256_and_si256(vmask, vsieve);

            vmask = _mm256_load_si256((__m256i *)(&_41_MASKS_256[t[4] * 4]));
            vsieve = _mm256_and_si256(vmask, vsieve);

            vmask = _mm256_load_si256((__m256i *)(&_43_MASKS_256[t[5] * 4]));
            vsieve = _mm256_and_si256(vmask, vsieve);

            vmask = _mm256_load_si256((__m256i *)(&_47_MASKS_256[t[6] * 4]));
            vsieve = _mm256_and_si256(vmask, vsieve);

            vmask = _mm256_load_si256((__m256i *)(&_53_MASKS_256[t[7] * 4]));
            vsieve = _mm256_and_si256(vmask, vsieve);

            _mm256_store_si256((__m256i *)(&flagblock64[k * 4]), vsieve);

            vmasknums = _mm256_sub_epi32(vmasknums, vmaskstep);
            vmasknums = _mm256_add_epi32(vmasknums, vmaskp);
            vmask = _mm256_cmpgt_epi32(vmasknums, vmaskp1);
            vmask = _mm256_and_si256(vmask, vmaskp);
            vmasknums = _mm256_sub_epi32(vmasknums, vmask);
        }

        _mm256_store_si256((__m256i *)(ddata->offsets + 8), vmasknums);

        vmaskp = _mm256_setr_epi32(59, 61, 67, 71, 73, 79, 83, 89);
        vmaskp1 = _mm256_setr_epi32(58, 60, 66, 70, 72, 78, 82, 88);
        vmaskstep = _mm256_setr_epi32(20, 12, 55, 43, 37, 19, 7, 78);
        vmasknums = _mm256_load_si256((__m256i *)(ddata->offsets + 16));

        for (k = 0; k < FLAGSIZE >> 8; k++)
        {
            __m256i vsieve = _mm256_load_si256((__m256i *)(&flagblock64[k * 4]));
            _mm256_store_si256((__m256i *)t, vmasknums);

            vmask = _mm256_load_si256((__m256i *)(&_59_MASKS_256[t[0] * 4]));
            vsieve = _mm256_and_si256(vmask, vsieve);

            vmask = _mm256_load_si256((__m256i *)(&_61_MASKS_256[t[1] * 4]));
            vsieve = _mm256_and_si256(vmask, vsieve);

            vmask = _mm256_load_si256((__m256i *)(&_67_MASKS_256[t[2] * 4]));
            vsieve = _mm256_and_si256(vmask, vsieve);

            vmask = _mm256_load_si256((__m256i *)(&_71_MASKS_256[t[3] * 4]));
            vsieve = _mm256_and_si256(vmask, vsieve);

            vmask = _mm256_load_si256((__m256i *)(&_73_MASKS_256[t[4] * 4]));
            vsieve = _mm256_and_si256(vmask, vsieve);

            vmask = _mm256_load_si256((__m256i *)(&_79_MASKS_256[t[5] * 4]));
            vsieve = _mm256_and_si256(vmask, vsieve);

            vmask = _mm256_load_si256((__m256i *)(&_83_MASKS_256[t[6] * 4]));
            vsieve = _mm256_and_si256(vmask, vsieve);

            vmask = _mm256_load_si256((__m256i *)(&_89_MASKS_256[t[7] * 4]));
            vsieve = _mm256_and_si256(vmask, vsieve);

            _mm256_store_si256((__m256i *)(&flagblock64[k * 4]), vsieve);

            vmasknums = _mm256_sub_epi32(vmasknums, vmaskstep);
            vmasknums = _mm256_add_epi32(vmasknums, vmaskp);
            vmask = _mm256_cmpgt_epi32(vmasknums, vmaskp1);
            vmask = _mm256_and_si256(vmask, vmaskp);
            vmasknums = _mm256_sub_epi32(vmasknums, vmask);
        }

        _mm256_store_si256((__m256i *)(ddata->offsets + 16), vmasknums);

        // continue presieving 16 more primes whose mask info has been
        // dynamically generated (to avoid super long constant lists).
        // we can't keep doing this forever because the memory footprint 
        // quickly grows too big and the flag density too sparse, negating 
        // the benefits of this approach.
        for (j = 24; j < sdata->presieve_max_id; j += 8)
        {
            vmaskp = _mm256_load_si256((__m256i *)(&presieve_primes[j - 24]));
            vmaskp1 = _mm256_load_si256((__m256i *)(&presieve_p1[j - 24]));
            vmaskstep = _mm256_load_si256((__m256i *)(&presieve_steps[j - 24]));
            vmasknums = _mm256_load_si256((__m256i *)(ddata->offsets + j));

            for (k = 0; k < FLAGSIZE >> 8; k++)
            {
                __m256i vsieve = _mm256_load_si256((__m256i *)(&flagblock64[k * 4]));
                _mm256_store_si256((__m256i *)t, vmasknums);

                vmask = _mm256_load_si256((__m256i *)(&presieve_largemasks[j - 24 + 0][t[0]]));
                vsieve = _mm256_and_si256(vmask, vsieve);

                vmask = _mm256_load_si256((__m256i *)(&presieve_largemasks[j - 24 + 1][t[1]]));
                vsieve = _mm256_and_si256(vmask, vsieve);

                vmask = _mm256_load_si256((__m256i *)(&presieve_largemasks[j - 24 + 2][t[2]]));
                vsieve = _mm256_and_si256(vmask, vsieve);

                vmask = _mm256_load_si256((__m256i *)(&presieve_largemasks[j - 24 + 3][t[3]]));
                vsieve = _mm256_and_si256(vmask, vsieve);

                vmask = _mm256_load_si256((__m256i *)(&presieve_largemasks[j - 24 + 4][t[4]]));
                vsieve = _mm256_and_si256(vmask, vsieve);

                vmask = _mm256_load_si256((__m256i *)(&presieve_largemasks[j - 24 + 5][t[5]]));
                vsieve = _mm256_and_si256(vmask, vsieve);

                vmask = _mm256_load_si256((__m256i *)(&presieve_largemasks[j - 24 + 6][t[6]]));
                vsieve = _mm256_and_si256(vmask, vsieve);

                vmask = _mm256_load_si256((__m256i *)(&presieve_largemasks[j - 24 + 7][t[7]]));
                vsieve = _mm256_and_si256(vmask, vsieve);

                _mm256_store_si256((__m256i *)(&flagblock64[k * 4]), vsieve);

                vmasknums = _mm256_sub_epi32(vmasknums, vmaskstep);
                vmasknums = _mm256_add_epi32(vmasknums, vmaskp);
                vmask = _mm256_cmpgt_epi32(vmasknums, vmaskp1);
                vmask = _mm256_and_si256(vmask, vmaskp);
                vmasknums = _mm256_sub_epi32(vmasknums, vmask);
            }

            _mm256_store_si256((__m256i *)(ddata->offsets + j), vmasknums);
        }

    }

    return;
}

#endif

void pre_sieve(soe_dynamicdata_t *ddata, soe_staticdata_t *sdata, uint8 *flagblock)
{
    uint64 startprime = sdata->startprime;
    uint64 k;
    int i, j;
    int mask_step, mask_step2, mask_step3, mask_step4;
    int mask_num, mask_num2, mask_num3, mask_num4;
    uint64 *flagblock64;

    // flagblock is always a multiple of 8 bytes
    flagblock64 = (uint64 *)flagblock;

    // do the smallest primes in predetermined batches.
    // first handle the different startprime cases, going up through prime 13.
    if (startprime == 2)
    {
        for (k = 0, mask_step = _64_MOD_P[0], mask_num = ddata->offsets[2],
            mask_step2 = _64_MOD_P[1], mask_num2 = ddata->offsets[3];
            k<FLAGSIZE >> 6; k++)
        {
            flagblock64[k] &= (_5_MASKS[mask_num] & _7_MASKS[mask_num2]);
            mask_num -= mask_step;
            if (mask_num < 0) mask_num = 5 + mask_num;
            mask_num2 -= mask_step2;
            if (mask_num2 < 0) mask_num2 = 7 + mask_num2;
        }
        ddata->offsets[2] = (uint32)mask_num;
        ddata->offsets[3] = (uint32)mask_num2;

        for (k = 0, mask_step = _64_MOD_P[2], mask_num = ddata->offsets[4],
            mask_step2 = _64_MOD_P[3], mask_num2 = ddata->offsets[5];
            k<FLAGSIZE >> 6; k++)
        {
            flagblock64[k] &= (_11_MASKS[mask_num] & _13_MASKS[mask_num2]);
            mask_num -= mask_step;
            if (mask_num < 0) mask_num = 11 + mask_num;
            mask_num2 -= mask_step2;
            if (mask_num2 < 0) mask_num2 = 13 + mask_num2;
        }
        ddata->offsets[4] = (uint32)mask_num;
        ddata->offsets[5] = (uint32)mask_num2;

    }
    else if (startprime == 3)
    {
        for (k = 0, mask_step = _64_MOD_P[1], mask_num = ddata->offsets[3];
            k<FLAGSIZE >> 6; k++)
        {
            flagblock64[k] &= _7_MASKS[mask_num];
            mask_num -= mask_step;
            if (mask_num < 0) mask_num = 7 + mask_num;
        }
        ddata->offsets[3] = (uint32)mask_num;

        for (k = 0, mask_step = _64_MOD_P[2], mask_num = ddata->offsets[4],
            mask_step2 = _64_MOD_P[3], mask_num2 = ddata->offsets[5];
            k<FLAGSIZE >> 6; k++)
        {
            flagblock64[k] &= (_11_MASKS[mask_num] & _13_MASKS[mask_num2]);
            mask_num -= mask_step;
            if (mask_num < 0) mask_num = 11 + mask_num;
            mask_num2 -= mask_step2;
            if (mask_num2 < 0) mask_num2 = 13 + mask_num2;
        }
        ddata->offsets[4] = (uint32)mask_num;
        ddata->offsets[5] = (uint32)mask_num2;

    }
    else if (startprime == 4)
    {
        for (k = 0, mask_step = _64_MOD_P[2], mask_num = ddata->offsets[4],
            mask_step2 = _64_MOD_P[3], mask_num2 = ddata->offsets[5];
            k<FLAGSIZE >> 6; k++)
        {
            flagblock64[k] &= (_11_MASKS[mask_num] & _13_MASKS[mask_num2]);
            mask_num -= mask_step;
            if (mask_num < 0) mask_num = 11 + mask_num;
            mask_num2 -= mask_step2;
            if (mask_num2 < 0) mask_num2 = 13 + mask_num2;
        }
        ddata->offsets[4] = (uint32)mask_num;
        ddata->offsets[5] = (uint32)mask_num2;

    }
    else if (startprime == 5)
    {
        for (k = 0, mask_step = _64_MOD_P[3], mask_num = ddata->offsets[5];
            k<FLAGSIZE >> 6; k++)
        {
            flagblock64[k] &= _13_MASKS[mask_num];
            mask_num -= mask_step;
            if (mask_num < 0) mask_num = 13 + mask_num;
        }
        ddata->offsets[5] = (uint32)mask_num;

    }

    // then do some more presieving, starting at prime 17.
    for (k = 0, mask_step = _64_MOD_P[4], mask_num = ddata->offsets[6],
        mask_step2 = _64_MOD_P[5], mask_num2 = ddata->offsets[7];
        k<FLAGSIZE >> 6; k++)
    {
        flagblock64[k] &= (_17_MASKS[mask_num] & _19_MASKS[mask_num2]);
        mask_num -= mask_step;
        if (mask_num < 0) mask_num = 17 + mask_num;
        mask_num2 -= mask_step2;
        if (mask_num2 < 0) mask_num2 = 19 + mask_num2;
    }
    ddata->offsets[6] = (uint32)mask_num;
    ddata->offsets[7] = (uint32)mask_num2;

    for (k = 0, mask_step = _64_MOD_P[6], mask_num = ddata->offsets[8],
        mask_step2 = _64_MOD_P[7], mask_num2 = ddata->offsets[9];
        k<FLAGSIZE >> 6; k++)
    {
        flagblock64[k] &= (_23_MASKS[mask_num] & _29_MASKS[mask_num2]);
        mask_num -= mask_step;
        if (mask_num < 0) mask_num = 23 + mask_num;
        mask_num2 -= mask_step2;
        if (mask_num2 < 0) mask_num2 = 29 + mask_num2;
    }
    ddata->offsets[8] = (uint32)mask_num;
    ddata->offsets[9] = (uint32)mask_num2;

    return;
}