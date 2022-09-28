# Copyright (C) 2001,2002 Jens Franke
# This file is part of gnfs4linux, distributed under the terms of the
# GNU General Public Licence and WITHOUT ANY WARRANTY.

# You should have received a copy of the GNU General Public License along
# with this program; see the file COPYING.  If not, write to the Free
# Software Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
# 02111-1307, USA.
#
# 6/13/04: Hacked up for use in GGNFS by Chris Monico.

INC=-I. -I./asm  -I../../../../gmp_install/gmp-6.2.0/include
LIBFLAGS=-L.  -L../../../../gmp_install/gmp-6.2.0/lib
#CFLAGS=-Os -march=amdfam10 -mtune=amdfam10 -funroll-loops
CFLAGS=-O3 -march=native -mtune=native -funroll-loops

.SUFFIXES:

.SECONDARY: *.c *.o

.PHONY: all clean

SRCFILES=fbgen.c gnfs-lasieve4e.c input-poly.c mpz-ull.c mpqs.c \
         real-poly-aux.c redu2.c gmp-aux.c if.c lasieve-prepn.c ../gmp-ecm/microecm.c tinyecm.c \
	 primgen32.c recurrence6.c asm/siever-config.c asm/siever-config.h 

OBJS=if.o input-poly.o redu2.o recurrence6.o fbgen.o ../gmp-ecm/microecm.o tinyecm.o \
     real-poly-aux.o primgen32.o lasieve-prepn.o mpqs.o

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

gnfs-lasieve4I11e: gnfs-lasieve4eI11.o $(OBJS) libgmp-aux.a \
                             liblasieve.a liblasieveI11.a 
	$(CC) $(CFLAGS) $(INC) $(LIBFLAGS) -o $@ $^ $(LIBS)
gnfs-lasieve4I12e: gnfs-lasieve4eI12.o $(OBJS) libgmp-aux.a \
                             liblasieve.a liblasieveI12.a 
	$(CC) $(CFLAGS) $(INC) $(LIBFLAGS) -o $@ $^ $(LIBS)
gnfs-lasieve4I13e: gnfs-lasieve4eI13.o $(OBJS) libgmp-aux.a \
                             liblasieve.a liblasieveI13.a 
	$(CC) $(CFLAGS) $(INC) $(LIBFLAGS) -o $@ $^ $(LIBS)
gnfs-lasieve4I14e: gnfs-lasieve4eI14.o $(OBJS) libgmp-aux.a \
                             liblasieve.a liblasieveI14.a 
	$(CC) $(CFLAGS) $(INC) $(LIBFLAGS) -o $@ $^ $(LIBS)
gnfs-lasieve4I15e: gnfs-lasieve4eI15.o $(OBJS) libgmp-aux.a \
                             liblasieve.a liblasieveI15.a 
	$(CC) $(CFLAGS) $(INC) $(LIBFLAGS) -o $@ $^ $(LIBS)
gnfs-lasieve4I16e: gnfs-lasieve4eI16.o $(OBJS) libgmp-aux.a \
                             liblasieve.a liblasieveI16.a 
	$(CC) $(CFLAGS) $(INC) $(LIBFLAGS) -o $@ $^ $(LIBS)


clean:
ifeq ($(OS),Windows_NT)
	-rm -f *.o *.a $(BINDIR)/gnfs-lasieve4I1?e.exe
else
	-rm -f *.o *.a $(BINDIR)/gnfs-lasieve4I1?e
endif
	(test -d asm && $(MAKE) -C asm clean) || exit 0
		

