# Copyright (C) 2001,2002 Jens Franke
# This file is part of gnfs4linux, distributed under the terms of the 
# GNU General Public Licence and WITHOUT ANY WARRANTY.

# You should have received a copy of the GNU General Public License along
# with this program; see the file COPYING.  If not, write to the Free
# Software Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
# 02111-1307, USA.


#gcc on osx
#CFLAGS= -m64 -Ofast -march=native -I/users/buhrow/src/c/gmp_install/gmp-6.2.0/include -Wa,-q
#CFLAGS= -m64 -Ofast -march=native -I/usr/local/Cellar/gmp4/4.3.2/include -Wa,-q

#windows
#CFLAGS= -Wall  -Wno-unused-variable -Wno-unused-function -Wno-unused-but-set-variable  -Ofast -fomit-frame-pointer -march=corei7 -mtune=corei7 -I. -Ic:/users/steve/progs/local/include
#CFLAGS= -Wimplicit-function-declaration -m64 -Ofast -march=native -Ic:/users/steve/progs/local/include

# mingw
# BRB: ULL_NO_UL didn't quite do the job by itself, because gmp ./configure set ulong to be 4 bytes
# but ulong is set to 8 bytes in siever_config.h when compiling lasieve.  Changing the siever_config
# definition breaks the montgomery asm code.  An easier solution is to set the following new definition:
# BITS_PER_GMP_ULONG, and a few relevant lines of code in mpz-ull.c
CFLAGS= -m64 -g -Ofast -march=native -I. -I../../../../gmp-install/mingw/include -DULL_NO_UL -DBITS_PER_GMP_ULONG=32 -fcommon



#clang
#CFLAGS= -Wall -Wno-unused-variable -Wno-unused-function -m64 -Ofast -march=native -I/Users/searle/progs/ensc-dependencies/include
#CFLAGS= -m64 -Ofast -march=corei7 -I/Users/searle/progs/ensc-dependencies/include 
#clang for profiling
#CFLAGS= -O1 -g -march=native -I/Users/searle/progs/ensc-dependencies/include

#linux
#CFLAGS= -m64 -g -Ofast -march=native -I. -I../../../../gmp-install/6.2.1/include -fcommon

include ../paths

CC=gcc $(CFLAGS)
#CC=gcc-4.9 $(CFLAGS)
#CC=clang $(CFLAGS)

ifeq ($(AVX512_LASCHED),1)
	CFLAGS += -DAVX512_LASCHED
endif

ifeq ($(AVX512_TD),1)
	CFLAGS += -DAVX512_TD
endif

CTANGLE=../../../../cweb-master/ctangle
#CTANGLE='c:/progra~2/cweb/bin/ctangle'

SRCFILES=gcd32.c 32bit.h modinv32.c psp.c Makefile \
	64bit.h asm_64bit.asm \
	siever-config.w lasched.w medsched.w \
	mpz-trialdiv.w asm-zeit.asm zeit.[ch] \
	ls-defs.asm lasched0.asm lasched1.asm lasched3.asm \
	lasieve-asmprepn2.asm lasieve-asmprepn.asm \
	lasieve-asmprepn64.asm lasieve-asmprepn2_64.asm \
	translate.asm rescale.asm \
	medsched0.asm MMX-TD.c MMX-TD1.asm MMX-TD2.asm modinv1000.asm \
	modinv1000b.asm mpz-td.asm pt64.asm ri-aux.asm \
	schedsieve.asm search0.asm \
	slinie3A.asm slinie2A.asm slinie1A.asm slinieG.asm \
	tdslinie1.asm tdslinie2.asm tdslinie3.asm tdslinie.asm \
	tdsieve-from-sched.asm \
	mpqs-config.h mpqs3_sieve.asm mpqs_eval.asm mpqs_gauss.asm \
	mpqs_nextpol.asm mpqs_sieve.asm mpqs_sieveinit.asm mpqs_td.asm \
	mpqs3arith.c invtab.c mpqs3_nextpol.asm mpqs3_td.asm \
	mpqs3_gauss.asm \
	noasm64.c basemath.c montgomery_mul.[ch] \
	asm_arith64.asm asm_arith128.asm asm_arith192.asm

%.c: %.w %.ch
	$(CTANGLE) $*.w $*.ch

%.c: %.w
	$(CTANGLE) $*.w

%.h: %.w %.ch
	$(CTANGLE) $*.w $*.ch

%.h: %.w
	$(CTANGLE) $*.w

%.tex: %.w %.ch
	cweave $*.w $*.ch

%.tex: %.w
	cweave $*.w

%.dvi: %.tex
	tex $<

%.S: %.asm
	m4 -g  -Dn_i_bits=0 ls-defs.asm $^ > $@

%.o: %.S
	$(CC) -c $^

%.o: %.s
	$(CC) -c $^

.SECONDARY:  *.s *.S

.PHONY: bup

all: liblasieve.a liblasieveI11.a liblasieveI12.a liblasieveI13.a liblasieveI14.a liblasieveI14.a liblasieveI15.a liblasieveI16.a

../gmp-aux.h:
	$(CTANGLE) ../gmp-aux.w

../if.h:
	$(CTANGLE) ../if.w

gcd32.o mpqs.o psp.o MMX-TD.o lasched.o: siever-config.h


liblasieve.a: gcd32.o psp.o pt64.o modinv1000.o ri-aux.o mpz-trialdiv.o \
	asm_64bit.o \
	modinv1000b.o schedsieve.o mpz-td.o MMX-TD.o MMX-TD1.o MMX-TD2.o \
	tdsieve-from-sched.o lasieve-asmprepn.o zeit.o asm-zeit.o \
	lasieve-asmprepn64.o lasieve-asmprepn2_64.o translate.o rescale.o \
	mpqs_td.o invtab.o mpqs_sieve.o mpqs_eval.o mpqs_sieveinit.o \
	mpqs_nextpol.o mpqs_gauss.o lasieve-asmprepn2.o \
	noasm64.o basemath.o montgomery_mul.o asm_arith64.o asm_arith128.o \
	asm_arith192.o mpqs3_sieve.o mpqs3_nextpol.o mpqs3_td.o \
	mpqs3_gauss.o
	ar rcs $@ $^

lasched0I%.s: ls-defs.asm lasched0.asm
	m4 -g -Dn_i_bits=`expr $* - 1` $^ > $@

lasched1I%.s: ls-defs.asm lasched1.asm
	m4 -g -Dn_i_bits=`expr $* - 1` -Dot=1 $^ > $@

lasched2I%.s: ls-defs.asm lasched1.asm
	m4 -g -Dn_i_bits=`expr $* - 1` -Dot=2 $^ > $@

lasched3I%.s: ls-defs.asm lasched1.asm
	m4 -g -Dn_i_bits=`expr $* - 1` -Dot=3 $^ > $@

lasched0ntI%.s: ls-defs.asm lasched0.asm
	m4 -g -Dn_i_bits=`expr $* - 1` -Dnt_sched=1 $^ > $@

lasched1ntI%.s: ls-defs.asm lasched1.asm
	m4 -g -Dn_i_bits=`expr $* - 1` -Dot=1 -Dnt_sched=1 $^ > $@

lasched2ntI%.s: ls-defs.asm lasched1.asm
	m4 -g -Dn_i_bits=`expr $* - 1` -Dot=2 -Dnt_sched=1 $^ > $@

lasched3ntI%.s: ls-defs.asm lasched1.asm
	m4 -g -Dn_i_bits=`expr $* - 1` -Dot=3 -Dnt_sched=1 $^ > $@

medsched0I%.s: ls-defs.asm medsched0.asm
	m4 -g -Dn_i_bits=`expr $* - 1` $^ > $@

laschedI%.o: lasched.c lasched.h siever-config.h ../if.h
	$(CC) $(CFLAGS) -DI_bits=$* -c -o $@ $<

medschedI%.o: medsched.c medsched.h siever-config.h ../if.h
	$(CC) $(CFLAGS) -DI_bits=$* -c -o $@ $<

slinieI%.s: ls-defs.asm slinieG.asm
	m4 -g -Dn_i_bits=`expr $* - 1` $^ > $@

slinie3I%.s: ls-defs.asm slinie3A.asm
	m4 -g -Dn_i_bits=`expr $* - 1` $^ > $@

slinie2I%.s: ls-defs.asm slinie2A.asm
	m4 -g -Dn_i_bits=`expr $* - 1` $^ > $@

slinie1I%.s: ls-defs.asm slinie1A.asm
	m4 -g -Dn_i_bits=`expr $* - 1` $^ > $@

tdslinieI%.s: ls-defs.asm tdslinie.asm
	m4 -g -Dn_i_bits=`expr $* - 1` $^ > $@

tdslinie1I%.s: ls-defs.asm tdslinie1.asm
	m4 -g -Dn_i_bits=`expr $* - 1` $^ > $@

tdslinie2I%.s: ls-defs.asm tdslinie2.asm
	m4 -g -Dn_i_bits=`expr $* - 1` $^ > $@

tdslinie3I%.s: ls-defs.asm tdslinie3.asm
	m4 -g -Dn_i_bits=`expr $* - 1` $^ > $@

search0I%.s: ls-defs.asm search0.asm
	m4 -g -Dn_i_bits=`expr $* - 1` $^ > $@

liblasieveI%.a:	laschedI%.o lasched0I%.o lasched1I%.o lasched2I%.o \
	lasched3I%.o lasched0ntI%.o lasched1ntI%.o lasched2ntI%.o \
	lasched3ntI%.o medschedI%.o medsched0I%.o search0I%.o \
	slinieI%.o slinie3I%.o slinie2I%.o slinie1I%.o \
	tdslinieI%.o tdslinie1I%.o tdslinie2I%.o tdslinie3I%.o 
	ar rcs $@ $^

bup:
	echo $(SRCFILES)

clean:
	rm *.o *.a
