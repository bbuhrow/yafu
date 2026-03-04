# Copyright (C) 2001,2002 Jens Franke
# This file is part of gnfs4linux, distributed under the terms of the 
# GNU General Public Licence and WITHOUT ANY WARRANTY.

# You should have received a copy of the GNU General Public License along
# with this program; see the file COPYING.  If not, write to the Free
# Software Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
# 02111-1307, USA.


# for people that have custom gmp locations or are evaluating a new
# gmp build in a custom location, fill in these fields to point to
# the include and lib locations.  Relative locations should be relative
# to this makefile location.  I have several to support various build
# environments.  Just comment them out or remove them as you wish.
CUSTOM_GMP_INC=../../../gmp-install/6.2.0-gcc/include
CUSTOM_GMP_LIB=../../../gmp-install/6.2.0-gcc/lib
#CUSTOM_GMP_INC=../../../gmp-install/6.2.0-gcc/include
#CUSTOM_GMP_LIB=../../../gmp-install/6.2.0-gcc/lib
#CUSTOM_GMP_INC=../../../gmp-install/mingw/include
#CUSTOM_GMP_LIB=../../../gmp-install/mingw/lib

CC=gcc

# gcc linux
CFLAGS=-Wall -Wno-unused-variable -Wno-unused-function -Ofast -fomit-frame-pointer -march=native -funroll-loops  \
	-I$(CUSTOM_GMP_INC) -L$(CUSTOM_GMP_LIB) -I../../include -fcommon

# mingw windows
# BRB: ULL_NO_UL didn't quite do the job by itself, because gmp ./configure set ulong to be 4 bytes
# but ulong is set to 8 bytes in siever_config.h when compiling lasieve.  Changing the siever_config
# definition breaks the montgomery asm code.  An easier solution is to set the following new definition:
# BITS_PER_GMP_ULONG, and a few relevant lines of code in mpz-ull.c
#CFLAGS=-Wall -Wno-unused-variable -Wno-unused-function -O2 -Ofast -DULL_NO_UL -DBITS_PER_GMP_ULONG=32 -fomit-frame-pointer -march=native -funroll-loops  \
#	-I. -I$(CUSTOM_GMP_INC) -L$(CUSTOM_GMP_LIB) -I../../include -fcommon
	
# clang linux
#CFLAGS=-Wall -Wno-unused-variable -Wno-unused-function -Ofast -fomit-frame-pointer -march=native -ffp-model=precise -funroll-loops  \
#	-I$(CUSTOM_GMP_INC) -L$(CUSTOM_GMP_LIB)  -fcommon

GMP_LIB=$(CUSTOM_GMP_LIB)/libgmp.a

include paths


#CTANGLE=ctangle
CTANGLE=~/src/c/cweb-master/ctangle


ifeq ($(AVX512_ECM),1)
	CFLAGS += -DAVX512_ECM
endif

ifeq ($(AVX512_LASETUP),1)
	CFLAGS += -DAVX512_LASIEVE_SETUP
endif

ifeq ($(AVX512_LASCHED),1)
	CFLAGS += -DAVX512_LASCHED
	ASM_FLAGS += AVX512_LASCHED=1
endif

ifeq ($(AVX512_TD),1)
	ifeq ($(CC),icc)
		CFLAGS += -DAVX512_TD
		ASM_FLAGS += AVX512_TD=1
	endif
endif

ifeq ($(AVX512_TDSCHED),1)
	CFLAGS += -DAVX512_TDSCHED
endif

ifeq ($(AVX512_SIEVE1),1)
	CFLAGS += -DAVX512_SIEVE1
endif

ifeq ($(AVX512_TDS0),1)
	CFLAGS += -DAVX512_TDS0
endif

ifeq ($(AVX512_SIEVE_SEARCH),1)
	CFLAGS += -DAVX512_SIEVE_SEARCH
endif

ifeq ($(LASTATS),1)
	CFLAGS += -DGATHER_STAT
endif

ifeq ($(MPQS),1)
	CFLAGS += -DGGNFS_MPQS
endif

ifeq ($(AVX512_ALL),1)
	ifeq ($(CC),icc)
		CFLAGS += -DUSE_AVX512F -DAVX512_LASIEVE_SETUP -DAVX512_LASCHED -DAVX512_TD -DAVX512_ECM -DAVX512_SIEVE1 -DAVX512_SIEVE_SEARCH -DAVX512_TDS0 -DAVX512_TDSCHED 
		ASM_FLAGS += AVX512_TD=1
	else
		CFLAGS += -DAVX512_LASIEVE_SETUP -DAVX512_ECM -DAVX512_LASCHED -DAVX512_SIEVE1 -DAVX512_SIEVE_SEARCH -DAVX512_TDS0 -DAVX512_TDSCHED 
		# -DAVX512_ECM 
	endif
	ASM_FLAGS += AVX512_LASCHED=1
	
endif


.SUFFIXES:

.SECONDARY: *.c *.o *.a

SRCFILES=fbgen.c fbgen.h lasieve-prepn.c if.c la-cs.c gmp-aux.c mpz-ull.c \
	fbgen64.c fbgen64.h \
	real-poly-aux.c real-poly-aux.h recurrence6.c redu2.c input-poly.c \
	input-poly-orig.c \
	primgen32.c primgen64.c gnfs-lasieve4e.c gnfs-lasieve4g.w \
	td.[ch] \
	ecm.[ch] strategy.c ecmtest.c \
	ecmstat.c mpqsstat.c mpqstest.c mpqs3test.c mpqs.c mpqs3.c tinyecm.c microecm.c \
	micropm1.c pprime_p.c pm1.[ch] pm1test.c pm1stat.c gnfs-lasieve4f.w


ASMDIRS=asm

list_asm_files = \
	$(foreach file, $(shell $(MAKE) -s -C $(asm_dir) bup),$(asm_dir)/$(file))

ASMFILES=$(foreach asm_dir,$(ASMDIRS),$(list_asm_files))

asm/%:	force
	$(MAKE) $(ASM_FLAGS) -C asm $*

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

%.s: %.asm
	gasp -c '#' -o $@ $^

%.o: %.c if.h asm/siever-config.h
	$(CC)  $(CFLAGS) -c -o $@ $<

%.o: %.s
	$(CC) $(CFLAGS) -c $^

.PHONY:	force

alle: gnfs-lasieve4I11e  gnfs-lasieve4I12e  gnfs-lasieve4I13e  gnfs-lasieve4I14e   gnfs-lasieve4I15e  gnfs-lasieve4I16e

allf: gnfs-lasieve4I11f  gnfs-lasieve4I12f  gnfs-lasieve4I13f  gnfs-lasieve4I14f   gnfs-lasieve4I15f  gnfs-lasieve4I16f

allg: gnfs-lasieve4I11g  gnfs-lasieve4I12g  gnfs-lasieve4I13g  gnfs-lasieve4I14g   gnfs-lasieve4I15g  gnfs-lasieve4I16g

input-poly.o: input-poly.h

input-poly-orig.o: input-poly-orig.h

fbgen.o: gmp-aux.h

fbgen.o lasieve-prepn.o: asm/32bit.h

recurrence6I%.o: recurrence6.c recurrence6.h if.h asm/siever-config.h
	$(CC)  $(CFLAGS) -c -DI_bits=$* -o $@ $<

redu2.o: redu2.h

primgen32.o fbgen.o: primgen32.h

gnfs-lasieve4e.c: la-cs.c

gnfs-lasieve4eI%.o: gnfs-lasieve4e.c if.h primgen32.h asm/32bit.h redu2.h \
	recurrence6.h fbgen.h real-poly-aux.h gmp-aux.h asm/medsched.h \
	asm/siever-config.h lasieve-prepn.h input-poly.h asm/lasched.h \
	strategy.h
	$(CC)  $(CFLAGS) -c -DI_bits=$* -o $@ $<

gnfs-lasieve4fI%.o: gnfs-lasieve4f.c if.h primgen32.h asm/32bit.h redu2.h \
	recurrence6.h fbgen.h real-poly-aux.h gmp-aux.h asm/medsched.h \
	asm/siever-config.h lasieve-prepn.h input-poly.h asm/lasched.h \
	strategy.h
	$(CC) $(CFLAGS) -c -DI_bits=$* -o $@ $<

gnfs-lasieve4gI%.o: gnfs-lasieve4g.c if.h primgen32.h asm/32bit.h redu2.h \
	recurrence6.h fbgen.h real-poly-aux.h gmp-aux.h asm/medsched.h \
	asm/siever-config.h lasieve-prepn.h input-poly.h asm/lasched.h \
	strategy.h
	$(CC) $(CFLAGS) -c -DI_bits=$* -o $@ $<

libgmp-aux.a: gmp-aux.o mpz-ull.o
	$(AR) rcs $@ $^

lasieve-prepn.o: lasieve-prepn.h recurrence6.h

mpqs3.o: asm/mpqs-config.h

mpqs.o: mpqs.c asm/mpqs-config.h asm/siever-config.h
#	gcc -c -O2 -o $@ $<
	$(CC) $(CFLAGS) -c -o $@ $<

mpqs3.o: mpqs3.c asm/mpqs-config.h asm/siever-config.h
#	gcc -c -O2 -o $@ $<
	$(CC) $(CFLAGS) -c -o $@ $<

gnfs-lasieve4I%e: gnfs-lasieve4eI%.o if.o input-poly.o libgmp-aux.a redu2.o \
	recurrence6I%.o fbgen.o fbgen64.o real-poly-aux.o mpqs.o \
	primgen32.o lasieve-prepn.o pprime_p.o \
	strategy.o ecm.o mpqs3.o pm1.o microecm.o micropm1.o tinyecm.o \
	asm/liblasieve.a asm/liblasieveI%.a
	$(CC) $(CFLAGS) -o $@ $^ $(GMP_LIB) -lm

gnfs-lasieve4I%f: gnfs-lasieve4fI%.o if.o input-poly-orig.o libgmp-aux.a redu2.o \
	recurrence6.o fbgen.o fbgen64.o real-poly-aux.o mpqs.o \
	primgen32.o lasieve-prepn.o pprime_p.o \
	strategy.o ecm.o mpqs3.o pm1.o \
	asm/liblasieve.a asm/liblasieveI%.a
	$(CC) $(CFLAGS) -o $@ $^ $(GMP_LIB) -lm

gnfs-lasieve4I%g: gnfs-lasieve4gI%.o if.o input-poly-orig.o libgmp-aux.a redu2.o \
	recurrence6.o fbgen.o fbgen64.o real-poly-aux.o mpqs.o \
	primgen32.o primgen64.o lasieve-prepn.o pprime_p.o \
	strategy.o ecm.o mpqs3.o pm1.o td.o \
	asm/liblasieve.a asm/liblasieveI%.a
	$(CC) $(CFLAGS) -o $@ $^ $(GMP_LIB) -lm

mpqsz.o: mpqs.c asm/mpqs-config.h
	$(CC) $(CFLAGS) -O2 -DMPQS_STAT -DMPQS_ZEIT -c -o $@ $<
#	gcc -g -DMPQS_STAT -DMPQS_ZEIT -c -o $@ $<

mpqstest: mpqstest.o mpqsz.o if.o mpz-ull.o asm/liblasieve.a
	$(CC) $(CFLAGS) -o $@ $^ $(GMP_LIB) -lm

mpqst: mpqst.o mpqsz.o if.o asm/liblasieve.a
	$(CC) $(CFLAGS) -o $@ $^ $(GMP_LIB) -lm

mpqszt.o: mpqs.c  asm/mpqs-config.h
	$(CC) $(CFLAGS) -O2 -DTOTAL_STAT -DMPQS_ZEIT -c -o $@ $<

tmpqs: mpqst.o mpqszt.o if.o asm/liblasieve.a
	$(CC) $(CFLAGS) -o $@ $^ $(GMP_LIB) -lm

mpqs3z.o: mpqs3.c asm/mpqs-config.h
	$(CC) $(CFLAGS) -O2 -DMPQS3_STAT -DMPQS3_ZEIT -c -o $@ $<
#	gcc -g -DMPQS3_STAT -DMPQS3_ZEIT -c -o $@ $<

mpqs3test: mpqs3test.o mpqs3z.o if.o libgmp-aux.a asm/liblasieve.a
	$(CC) $(CFLAGS) -o $@ $^ $(GMP_LIB) -lm

mpqs3z2.o: mpqs3.c asm/mpqs-config2.h
	$(CC) $(CFLAGS) -O2 -DMPQS3_STAT -DMPQS3_ZEIT -DVARIANT2 -c -o $@ $<
#	gcc -g -DMPQS3_STAT -DMPQS3_ZEIT -c -o $@ $<

mpqs3test2: mpqs3test.o mpqs3z2.o if.o libgmp-aux.a asm/liblasieve.a
	$(CC) $(CFLAGS) -o $@ $^ $(GMP_LIB) -lm

mpqsstat: mpqsstat.o mpqs.o mpqs3.o if.o mpz-ull.o libgmp-aux.a \
	asm/liblasieve.a
	$(CC) $(CFLAGS) -o $@ $^ $(GMP_LIB) -lm

ecmz.o: ecm.c asm/siever-config.h
	$(CC) $(CFLAGS) -DECM_STAT -DECM_ZEIT -c -o $@ $<
	
ecmtest: ecmtest.o ecmz.o if.o asm/liblasieve.a
	$(CC) $(CFLAGS) -o $@ $^ $(GMP_LIB) -lm

ecmstat: ecmstat.o ecmz.o if.o asm/liblasieve.a
	$(CC) $(CFLAGS) -o $@ $^ $(GMP_LIB) -lm

pm1z.o: pm1.c asm/siever-config.h
	$(CC) $(CFLAGS) -DPM1_STAT -DPM1_ZEIT -c -o $@ $<

pm1test: pm1test.o pm1z.o if.o asm/liblasieve.a
	$(CC) $(CFLAGS) -o $@ $^ $(GMP_LIB) -lm

pm1stat: pm1stat.o pm1z.o if.o asm/liblasieve.a
	$(CC) $(CFLAGS) -o $@ $^ $(GMP_LIB) -lm

lasieve5.tgz: Makefile paths INSTALL.and.USE COPYING $(SRCFILES) $(ASMFILES)
	tar cvO $^ | gzip --best --stdout > $@

clean:
	rm -f *.o *.a asm/*.o asm/*.a asm/*.s asm/*.S

binclean: realclean
	rm -f gnfs-lasieve4I1* forumexs/*.out readmeex/*.gz
