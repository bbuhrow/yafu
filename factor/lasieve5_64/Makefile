# Copyright (C) 2001,2002 Jens Franke
# This file is part of gnfs4linux, distributed under the terms of the 
# GNU General Public Licence and WITHOUT ANY WARRANTY.

# You should have received a copy of the GNU General Public License along
# with this program; see the file COPYING.  If not, write to the Free
# Software Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
# 02111-1307, USA.

#ifdef DEBUG
#CFLAGS= -DGATHER_STAT -DDEBUG -g
#else

# SMJS various ones I've used

# For gcc installed using brew on osx 10.9 (note -Wa,-q which specifies using clang as rather than system one)
#CFLAGS= -m64 -Ofast -march=native -fomit-frame-pointer -funroll-loops -I/Users/searle/progs/ensc-dependencies/include -Wa,-q
#CFLAGS= -m64 -Ofast -march=corei7 -mtune=corei7 -I/Users/searle/progs/ensc-dependencies/include -Wa,-q

# windows one
#CFLAGS=-Wall -Wno-unused-variable -Wno-unused-function -Wno-unused-but-set-variable -Ofast -fomit-frame-pointer -march=corei7 -mtune=corei7 -funroll-loops -Ic:/users/steve/progs/local/include 

# linux
CFLAGS=-Wall -Wno-unused-variable -Wno-unused-function -Ofast -fomit-frame-pointer -march=native -funroll-loops  \
	-I../../../gmp_install/gmp-6.2.0/include -I../../include -fcommon

#clang
#CFLAGS= -Wall -Wno-unused-variable -Wno-unused-function -Ofast -march=native -I/Users/searle/progs/ensc-dependencies/include
#clang for profiling
#CFLAGS= -O1 -g -march=native -I/Users/searle/progs/ensc-dependencies/include

#endif


GMP_LIB=-L../../../gmp_install/gmp-6.2.0/include/lib -lgmp
#For windows
#GMP_LIB=-Lc:/Users/steve/progs/local/lib -lgmp -lws2_32
#Using brew installed gmp
#GMP_LIB=-L/usr/local/Cellar/gmp4/4.3.2/lib -lgmp
#orig
#GMP_LIB=/home/franke/itanium-bin/libgmp.a

include paths

CC=gcc  $(CFLAGS)
#CC=gcc-4.9  $(CFLAGS)
#CC=clang  $(CFLAGS)

#CTANGLE=ctangle
CTANGLE=../../../cweb-master/ctangle


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
		CFLAGS += -DAVX512_LASIEVE_SETUP -DAVX512_LASCHED -DAVX512_TD -DAVX512_ECM -DAVX512_SIEVE1 -DAVX512_SIEVE_SEARCH -DAVX512_TDS0 -DAVX512_TDSCHED 
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
	ecmstat.c mpqsstat.c mpqstest.c mpqs3test.c mpqs.c mpqs3.c tinyecm.c ../gmp-ecm/microecm.c \
	../gmp-ecm/micropm1.c pprime_p.c pm1.[ch] pm1test.c pm1stat.c gnfs-lasieve4f.w


ASMDIRS=athlon athlon64 piii generic xeon64

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
	$(CC) -c $^

.PHONY:	force

alle: gnfs-lasieve4I11e  gnfs-lasieve4I12e  gnfs-lasieve4I13e  gnfs-lasieve4I14e   gnfs-lasieve4I15e  gnfs-lasieve4I16e

allf: gnfs-lasieve4I11f  gnfs-lasieve4I12f  gnfs-lasieve4I13f  gnfs-lasieve4I14f   gnfs-lasieve4I15f  gnfs-lasieve4I16f

allg: gnfs-lasieve4I11g  gnfs-lasieve4I12g  gnfs-lasieve4I13g  gnfs-lasieve4I14g   gnfs-lasieve4I15g  gnfs-lasieve4I16g

input-poly.o: input-poly.h

input-poly-orig.o: input-poly-orig.h

fbgen.o: gmp-aux.h

fbgen.o lasieve-prepn.o: asm/32bit.h

recurrence6.o: recurrence6.c recurrence6.h if.h asm/siever-config.h
	$(CC)  $(CFLAGS) -c -o $@ $<

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
	$(CC) -c -DI_bits=$* -o $@ $<

gnfs-lasieve4gI%.o: gnfs-lasieve4g.c if.h primgen32.h asm/32bit.h redu2.h \
	recurrence6.h fbgen.h real-poly-aux.h gmp-aux.h asm/medsched.h \
	asm/siever-config.h lasieve-prepn.h input-poly.h asm/lasched.h \
	strategy.h
	$(CC) -c -DI_bits=$* -o $@ $<

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
	recurrence6.o fbgen.o fbgen64.o real-poly-aux.o mpqs.o \
	primgen32.o lasieve-prepn.o pprime_p.o \
	strategy.o ecm.o mpqs3.o pm1.o ../gmp-ecm/microecm.o ../gmp-ecm/micropm1.o tinyecm.o \
	asm/liblasieve.a asm/liblasieveI%.a
	$(CC) -o $@ $^ $(GMP_LIB) -lm

gnfs-lasieve4I%f: gnfs-lasieve4fI%.o if.o input-poly-orig.o libgmp-aux.a redu2.o \
	recurrence6.o fbgen.o fbgen64.o real-poly-aux.o mpqs.o \
	primgen32.o lasieve-prepn.o pprime_p.o \
	strategy.o ecm.o mpqs3.o pm1.o \
	asm/liblasieve.a asm/liblasieveI%.a
	$(CC) -o $@ $^ $(GMP_LIB) -lm

gnfs-lasieve4I%g: gnfs-lasieve4gI%.o if.o input-poly-orig.o libgmp-aux.a redu2.o \
	recurrence6.o fbgen.o fbgen64.o real-poly-aux.o mpqs.o \
	primgen32.o primgen64.o lasieve-prepn.o pprime_p.o \
	strategy.o ecm.o mpqs3.o pm1.o td.o \
	asm/liblasieve.a asm/liblasieveI%.a
	$(CC) -o $@ $^ $(GMP_LIB) -lm

mpqsz.o: mpqs.c asm/mpqs-config.h
	$(CC) -O2 -DMPQS_STAT -DMPQS_ZEIT -c -o $@ $<
#	gcc -g -DMPQS_STAT -DMPQS_ZEIT -c -o $@ $<

mpqstest: mpqstest.o mpqsz.o if.o mpz-ull.o asm/liblasieve.a
	$(CC) -o $@ $^ $(GMP_LIB) -lm

mpqst: mpqst.o mpqsz.o if.o asm/liblasieve.a
	$(CC) -o $@ $^ $(GMP_LIB) -lm

mpqszt.o: mpqs.c  asm/mpqs-config.h
	$(CC) -O2 -DTOTAL_STAT -DMPQS_ZEIT -c -o $@ $<

tmpqs: mpqst.o mpqszt.o if.o asm/liblasieve.a
	$(CC) -o $@ $^ $(GMP_LIB) -lm

mpqs3z.o: mpqs3.c asm/mpqs-config.h
	$(CC) -O2 -DMPQS3_STAT -DMPQS3_ZEIT -c -o $@ $<
#	gcc -g -DMPQS3_STAT -DMPQS3_ZEIT -c -o $@ $<

mpqs3test: mpqs3test.o mpqs3z.o if.o libgmp-aux.a asm/liblasieve.a
	$(CC) -o $@ $^ $(GMP_LIB) -lm

mpqs3z2.o: mpqs3.c asm/mpqs-config2.h
	$(CC) -O2 -DMPQS3_STAT -DMPQS3_ZEIT -DVARIANT2 -c -o $@ $<
#	gcc -g -DMPQS3_STAT -DMPQS3_ZEIT -c -o $@ $<

mpqs3test2: mpqs3test.o mpqs3z2.o if.o libgmp-aux.a asm/liblasieve.a
	$(CC) -o $@ $^ $(GMP_LIB) -lm

mpqsstat: mpqsstat.o mpqs.o mpqs3.o if.o mpz-ull.o libgmp-aux.a \
	asm/liblasieve.a
	$(CC) -o $@ $^ $(GMP_LIB) -lm

ecmz.o: ecm.c asm/siever-config.h
	$(CC) -DECM_STAT -DECM_ZEIT -c -o $@ $<
	
ecmtest: ecmtest.o ecmz.o if.o asm/liblasieve.a
	$(CC) -o $@ $^ $(GMP_LIB) -lm

ecmstat: ecmstat.o ecmz.o if.o asm/liblasieve.a
	$(CC) -o $@ $^ $(GMP_LIB) -lm

pm1z.o: pm1.c asm/siever-config.h
	$(CC) -DPM1_STAT -DPM1_ZEIT -c -o $@ $<

pm1test: pm1test.o pm1z.o if.o asm/liblasieve.a
	$(CC) -o $@ $^ $(GMP_LIB) -lm

pm1stat: pm1stat.o pm1z.o if.o asm/liblasieve.a
	$(CC) -o $@ $^ $(GMP_LIB) -lm

lasieve5.tgz: Makefile paths INSTALL.and.USE COPYING $(SRCFILES) $(ASMFILES)
	tar cvO $^ | gzip --best --stdout > $@

clean:
	rm -f *.o *.a asm/*.o asm/*.a asm/*.s asm/*.S

binclean: realclean
	rm -f gnfs-lasieve4I1* forumexs/*.out readmeex/*.gz
