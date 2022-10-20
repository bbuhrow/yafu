# Copyright (C) 2001,2002 Jens Franke
# This file is part of gnfs4linux, distributed under the terms of the 
# GNU General Public Licence and WITHOUT ANY WARRANTY.

# You should have received a copy of the GNU General Public License along
# with this program; see the file COPYING.  If not, write to the Free
# Software Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
# 02111-1307, USA.

ifdef DEBUG
CFLAGS= -DGATHER_STAT -DDEBUG -g
else
CFLAGS= -O3 -DGATHER_STAT
endif
#CFLAGS= -DGATHER_STAT -DDEBUG -g

GMP_LIB=-lgmp
# GMP_LIB=/home/franke/itanium-bin/libgmp.a

include paths

CC=gcc -static $(CFLAGS)
TANGLE=ctangle

.SUFFIXES:

.SECONDARY: *.c *.o *.a

SRCFILES=fbgen.c fbgen.h lasieve-prepn.w la-cs.w if.w gmp-aux.w mpz-ull.w \
	fbgen64.c fbgen64.h \
	real-poly-aux.c real-poly-aux.h recurrence6.w redu2.w input-poly.w \
	primgen32.w primgen64.w gnfs-lasieve4e.w gnfs-lasieve4g.w \
	td.[ch] \
	ecm.[ch] strategy.w ecmtest.c \
	ecmstat.c mpqsstat.c mpqstest.c mpqs3test.c mpqs.c mpqs3.c \
	pprime_p.c pm1.[ch] pm1test.c pm1stat.c gnfs-lasieve4f.w


ASMDIRS=athlon athlon64 piii generic xeon64

list_asm_files = \
	$(foreach file, $(shell $(MAKE) -s -C $(asm_dir) bup),$(asm_dir)/$(file))

ASMFILES=$(foreach asm_dir,$(ASMDIRS),$(list_asm_files))

asm/%:	force
	$(MAKE) -C asm $*

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
	$(CC) -c -o $@ $<

%.o: %.s
	$(CC) -c $^

.PHONY:	force

clean:
	rm *.o *.a asm/*.o asm/*.a

input-poly.o: input-poly.h

fbgen.o: gmp-aux.h

fbgen.o lasieve-prepn.o: asm/32bit.h

recurrence6.o: recurrence6.c recurrence6.h if.h asm/siever-config.h
	$(CC) -c -o $@ $<

redu2.o: redu2.h

primgen32.o fbgen.o: primgen32.h

gnfs-lasieve4e.c: la-cs.w

gnfs-lasieve4eI%.o: gnfs-lasieve4e.c if.h primgen32.h asm/32bit.h redu2.h \
	recurrence6.h fbgen.h real-poly-aux.h gmp-aux.h asm/medsched.h \
	asm/siever-config.h lasieve-prepn.h input-poly.h asm/lasched.h \
	strategy.h
	$(CC) -c -DI_bits=$* -o $@ $<

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
	$(CC) -c -o $@ $<

mpqs3.o: mpqs3.c asm/mpqs-config.h asm/siever-config.h
#	gcc -c -O2 -o $@ $<
	$(CC) -c -o $@ $<

gnfs-lasieve4I%e: gnfs-lasieve4eI%.o if.o input-poly.o redu2.o \
	recurrence6.o fbgen.o fbgen64.o real-poly-aux.o mpqs.o \
	primgen32.o libgmp-aux.a lasieve-prepn.o pprime_p.o \
	strategy.o ecm.o mpqs3.o pm1.o \
	asm/liblasieve.a asm/liblasieveI%.a
	$(CC) -o $@ $^ $(GMP_LIB) -lm

gnfs-lasieve4I%f: gnfs-lasieve4fI%.o if.o input-poly.o redu2.o \
	recurrence6.o fbgen.o fbgen64.o real-poly-aux.o mpqs.o \
	primgen32.o libgmp-aux.a lasieve-prepn.o pprime_p.o \
	strategy.o ecm.o mpqs3.o pm1.o \
	asm/liblasieve.a asm/liblasieveI%.a
	$(CC) -o $@ $^ $(GMP_LIB) -lm

gnfs-lasieve4I%g: gnfs-lasieve4gI%.o if.o input-poly.o redu2.o \
	recurrence6.o fbgen.o fbgen64.o real-poly-aux.o mpqs.o \
	primgen32.o primgen64.o libgmp-aux.a lasieve-prepn.o pprime_p.o \
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

