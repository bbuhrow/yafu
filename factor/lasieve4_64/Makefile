# Copyright (C) 2001,2002 Jens Franke
# This file is part of gnfs4linux, distributed under the terms of the
# GNU General Public Licence and WITHOUT ANY WARRANTY.

# You should have received a copy of the GNU General Public License along
# with this program; see the file COPYING.  If not, write to the Free
# Software Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
# 02111-1307, USA.
#
# 6/13/04: Hacked up for use in GGNFS by Chris Monico.
# 11/11/22: Ben Buhrow, added AVX512 functionality, compile options, and fast ecm cofactoring

INC=-I. -I./asm  -I../../../gmp_install/gmp-6.2.0/include -I../../include
LIBFLAGS=-L.  -L../../../gmp_install/gmp-6.2.0/lib
#CFLAGS=-Os -march=amdfam10 -mtune=amdfam10 -funroll-loops
CFLAGS=-O3 -g -march=native -mtune=native -funroll-loops -static-intel
# -DCONTIGUOUS_RI

ifeq ($(AVX512_ECM),1)
	CFLAGS += -DAVX512_ECM
endif

ifeq ($(AVX512_LASETUP),1)
	CFLAGS += -DAVX512_LASIEVE_SETUP
endif

ifeq ($(AVX512_LASCHED),1)
	CFLAGS += -DAVX512_LASCHED
endif

ifeq ($(AVX512_TD),1)
	CFLAGS += -DAVX512_TD
endif

ifeq ($(AVX512_SIEVE1),1)
	CFLAGS += -DAVX512_SIEVE1
endif

ifeq ($(LASTATS),1)
	CFLAGS += -DGATHER_STAT
endif

ifeq ($(MPQS),1)
	CFLAGS += -DGGNFS_MPQS
endif

ifeq ($(AVX512_ALL),1)
	CFLAGS += -DAVX512_LASIEVE_SETUP -DAVX512_LASCHED -DAVX512_TD -DAVX512_SIEVE1 -DAVX512_ECM 
endif


.SUFFIXES:

.SECONDARY: *.c *.o

.PHONY: all clean

SRCFILES=fbgen.c gnfs-lasieve4e.c input-poly.c mpz-ull.c mpqs.c \
         real-poly-aux.c redu2.c gmp-aux.c if.c lasieve-prepn.c ../gmp-ecm/microecm.c tinyecm.c \
	 primgen32.c recurrence6.c asm/siever-config.c asm/siever-config.h 

OBJS=if.o input-poly.o redu2.o recurrence6.o fbgen.o ../gmp-ecm/microecm.o tinyecm.o \
     real-poly-aux.o primgen32.o lasieve-prepn.o mpqs.o

ifeq ($(AVX512_LASCHED),1)
	SRCFILES += asm/lasched.c asm/medsched.c
endif


LIBS=-lgmp-aux -lgmp -lm


all : gnfs-lasieve4I11e gnfs-lasieve4I12e gnfs-lasieve4I13e \
      gnfs-lasieve4I14e gnfs-lasieve4I15e gnfs-lasieve4I16e


%.o: %.c asm/siever-config.h
	$(CC) $(CFLAGS) $(INC) -c -o $@ $<

%.o: %.s
	$(CC) $(INC) -c $^

libgmp-aux.a: gmp-aux.o mpz-ull.o
	$(AR) rcs $@ $^

gnfs-lasieve4eI%.o: gnfs-lasieve4e.c 
	$(CC) $(CFLAGS) $(INC) -c -DI_bits=$* -o $@ $<

asm/lasched%.o: asm/lasched.c
	$(CC) $(CFLAGS) $(INC) -c -DI_bits=$* -o $@ $<

asm/medsched%.o: asm/medsched.c
	$(CC) $(CFLAGS) $(INC) -c -DI_bits=$* -o $@ $<
	
gnfs-lasieve4I11e: gnfs-lasieve4eI11.o asm/lasched11.o asm/medsched11.o $(OBJS) libgmp-aux.a \
                             liblasieve.a liblasieveI11.a 
	$(CC) $(CFLAGS) $(INC) $(LIBFLAGS) -o $@ $^ $(LIBS)
gnfs-lasieve4I12e: gnfs-lasieve4eI12.o asm/lasched12.o asm/medsched12.o $(OBJS) libgmp-aux.a \
                             liblasieve.a liblasieveI12.a 
	$(CC) $(CFLAGS) $(INC) $(LIBFLAGS) -o $@ $^ $(LIBS)
gnfs-lasieve4I13e: gnfs-lasieve4eI13.o asm/lasched13.o asm/medsched13.o $(OBJS) libgmp-aux.a \
                             liblasieve.a liblasieveI13.a 
	$(CC) $(CFLAGS) $(INC) $(LIBFLAGS) -o $@ $^ $(LIBS)
gnfs-lasieve4I14e: gnfs-lasieve4eI14.o asm/lasched14.o asm/medsched14.o $(OBJS) libgmp-aux.a \
                             liblasieve.a liblasieveI14.a 
	$(CC) $(CFLAGS) $(INC) $(LIBFLAGS) -o $@ $^ $(LIBS)
gnfs-lasieve4I15e: gnfs-lasieve4eI15.o asm/lasched15.o asm/medsched15.o $(OBJS) libgmp-aux.a \
                             liblasieve.a liblasieveI15.a 
	$(CC) $(CFLAGS) $(INC) $(LIBFLAGS) -o $@ $^ $(LIBS)
gnfs-lasieve4I16e: gnfs-lasieve4eI16.o asm/lasched16.o asm/medsched16.o $(OBJS) libgmp-aux.a \
                             liblasieve.a liblasieveI16.a 
	$(CC) $(CFLAGS) $(INC) $(LIBFLAGS) -o $@ $^ $(LIBS)


clean:
ifeq ($(OS),Windows_NT)
	-rm -f *.o *.a $(BINDIR)/gnfs-lasieve4I1?e.exe
else
	-rm -f *.o *.a $(BINDIR)/gnfs-lasieve4I1?e
endif
	(test -d asm && $(MAKE) -C asm clean) || exit 0
		

