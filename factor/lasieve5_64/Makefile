# Copyright (C) 2001,2002 Jens Franke
# This file is part of gnfs4linux, distributed under the terms of the
# GNU General Public Licence and WITHOUT ANY WARRANTY.
#
# You should have received a copy of the GNU General Public License along
# with this program; see the file COPYING.  If not, write to the Free
# Software Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
# 02111-1307, USA.
#
# -----------------------------------------------------------------------------
# factor/lasieve5_64/Makefile
#
# Builds the NFS sieve executables:
#   gnfs-lasieve4I{11..16}e   (alle)
#   gnfs-lasieve4I{11..16}f   (allf)
#   gnfs-lasieve4I{11..16}g   (allg)
#
# Assembly libraries (liblasieve.a, liblasieveI%.a) are built in the asm/
# subdirectory and pulled in automatically via the asm/%: force delegation.
#
# Variables accepted from the top-level yafu Makefile (all have sensible
# defaults so this file can also be used standalone):
#
#   CC            Compiler: gcc, clang, or icc (default: gcc)
#   GMP_INCDIR    Path to gmp.h               (default: /usr/include)
#   GMP_LIBDIR    Path to libgmp.a            (default: /usr/lib)
#   DETECTED_OS   Windows or Linux            (auto-detected if absent)
#   DEBUG         Set to 1 for a debug build
#
# AVX-512 feature flags (set individually or use AVX512_ALL=1 for all):
#   AVX512_ALL        Enable all AVX-512 sieve paths
#   AVX512_ECM        AVX-512 ECM
#   AVX512_LASETUP    AVX-512 lasieve setup
#   AVX512_LASCHED    AVX-512 scheduler  (also forwarded to asm/)
#   AVX512_TD         AVX-512 trial division  (ICC only)
#   AVX512_TDSCHED    AVX-512 TD scheduler
#   AVX512_SIEVE1     AVX-512 sieve pass 1
#   AVX512_TDS0       AVX-512 TDS0
#   AVX512_SIEVE_SEARCH  AVX-512 sieve search
#   LASTATS           Gather statistics
#   MPQS              Enable MPQS
#
# Standalone usage examples:
#   make alle
#   make alle AVX512_ALL=1
#   make alle GMP_INCDIR=/path/to/include GMP_LIBDIR=/path/to/lib
#   make clean
# -----------------------------------------------------------------------------


# --- Standalone fallbacks (overridden when called from top-level) ------------
CC         ?= gcc
GMP_INCDIR ?= /usr/include
GMP_LIBDIR ?= /usr/lib

ifndef DETECTED_OS
    ifeq ($(OS),Windows_NT)
        DETECTED_OS := Windows
    else
        DETECTED_OS := $(shell uname -s)
    endif
endif

# Compiler family — mirrors the top-level detection so ifeq ($(CC),icc) checks
# work correctly even when CC is a full path or versioned name.
CC_BASENAME := $(notdir $(CC))
ifneq (,$(findstring icx,$(CC_BASENAME)))
    COMPILER_FAMILY := icc
else ifneq (,$(findstring icc,$(CC_BASENAME)))
    COMPILER_FAMILY := icc
else ifneq (,$(findstring clang,$(CC_BASENAME)))
    COMPILER_FAMILY := clang
else
    COMPILER_FAMILY := gcc
endif


# --- GMP static library ------------------------------------------------------
# Lasieve links statically against GMP so the sieve binaries are self-contained.
GMP_LIB := $(GMP_LIBDIR)/libgmp.a


# --- CFLAGS: platform and compiler selection ---------------------------------
#
# All three variants share a common base; platform and compiler add extras.
#
# gcc linux  : base flags only
# clang linux: base + -ffp-model=precise (required for reproducible FP results)
# mingw win  : base + -I. + -DULL_NO_UL -DBITS_PER_GMP_ULONG=32
#              (GMP's ./configure sets ulong=32-bit on Windows, but lasieve
#              assumes 64-bit; the defines reconcile this without changing
#              siever-config.h or the montgomery asm code)
#
CFLAGS_BASE := \
    -Wall -Wno-unused-variable -Wno-unused-function \
    -Ofast -fomit-frame-pointer -march=native -funroll-loops \
    -I$(GMP_INCDIR) -L$(GMP_LIBDIR) -I../../include -fcommon

ifeq ($(DETECTED_OS),Windows)
    CFLAGS = -I. $(CFLAGS_BASE) -DULL_NO_UL -DBITS_PER_GMP_ULONG=32
else ifeq ($(COMPILER_FAMILY),clang)
    CFLAGS = $(CFLAGS_BASE) -ffp-model=precise
else
    CFLAGS = $(CFLAGS_BASE)
endif

ifeq ($(DEBUG),1)
    CFLAGS := $(filter-out -Ofast -fomit-frame-pointer -march=native -funroll-loops,$(CFLAGS))
    CFLAGS += -O0 -g
endif


# --- ASM_FLAGS: variables forwarded to the asm/ sub-make --------------------
# Always forward the compiler, GMP path, and OS so the asm Makefile doesn't
# need its own hardcoded values.  AVX512 flags are appended below as needed.
ASM_FLAGS := CC=$(CC) GMP_INCDIR=$(GMP_INCDIR) DETECTED_OS=$(DETECTED_OS) \
             $(if $(DEBUG),DEBUG=$(DEBUG))


# --- AVX-512 feature flags ---------------------------------------------------
# Individual flags can be set on the command line or in config.mk (via yafu's
# top-level Makefile which passes AVX512_ALL=1 when USE_AVX512=1 is set).
# AVX512_TD is ICC-only: that compiler generates faster AVX-512 TD code.

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
    ifeq ($(COMPILER_FAMILY),icc)
        CFLAGS    += -DAVX512_TD
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
    ifeq ($(COMPILER_FAMILY),icc)
        CFLAGS    += -DUSE_AVX512F -DAVX512_LASIEVE_SETUP -DAVX512_LASCHED \
                     -DAVX512_TD -DAVX512_ECM -DAVX512_SIEVE1 \
                     -DAVX512_SIEVE_SEARCH -DAVX512_TDS0 -DAVX512_TDSCHED
        ASM_FLAGS += AVX512_TD=1
    else
        # gcc / clang: all AVX-512 paths except TD (which requires ICC codegen)
        CFLAGS    += -DAVX512_LASIEVE_SETUP -DAVX512_ECM -DAVX512_LASCHED \
                     -DAVX512_SIEVE1 -DAVX512_SIEVE_SEARCH \
                     -DAVX512_TDS0 -DAVX512_TDSCHED
    endif
    ASM_FLAGS += AVX512_LASCHED=1
endif


# --- Build rules -------------------------------------------------------------

.SUFFIXES:
.SECONDARY: *.c *.o *.a

ASMDIRS = asm

list_asm_files = \
    $(foreach file,$(shell $(MAKE) -s -C $(asm_dir) bup),$(asm_dir)/$(file))

ASMFILES = $(foreach asm_dir,$(ASMDIRS),$(list_asm_files))

# Delegate any asm/* target to the asm sub-make, forwarding ASM_FLAGS.
asm/%: force
	$(MAKE) $(ASM_FLAGS) -C asm $*

# Generic compile rules
%.o: %.c if.h asm/siever-config.h
	$(CC) $(CFLAGS) -c -o $@ $<

%.o: %.s
	$(CC) $(CFLAGS) -c $^

.PHONY: force alle allf allg clean


# --- Targets -----------------------------------------------------------------

alle: gnfs-lasieve4I11e gnfs-lasieve4I12e gnfs-lasieve4I13e \
      gnfs-lasieve4I14e gnfs-lasieve4I15e gnfs-lasieve4I16e

allf: gnfs-lasieve4I11f gnfs-lasieve4I12f gnfs-lasieve4I13f \
      gnfs-lasieve4I14f gnfs-lasieve4I15f gnfs-lasieve4I16f

allg: gnfs-lasieve4I11g gnfs-lasieve4I12g gnfs-lasieve4I13g \
      gnfs-lasieve4I14g gnfs-lasieve4I15g gnfs-lasieve4I16g


# --- Explicit header dependencies --------------------------------------------

input-poly.o: input-poly.h
input-poly-orig.o: input-poly-orig.h
fbgen.o: gmp-aux.h
fbgen.o lasieve-prepn.o: asm/32bit.h
redu2.o: redu2.h
primgen32.o fbgen.o: primgen32.h
gnfs-lasieve4e.c: la-cs.c
lasieve-prepn.o: lasieve-prepn.h recurrence6.h
mpqs3.o: asm/mpqs-config.h


# --- Per-I-value compile rules -----------------------------------------------

recurrence6I%.o: recurrence6.c recurrence6.h if.h asm/siever-config.h
	$(CC) $(CFLAGS) -c -DI_bits=$* -o $@ $<

gnfs-lasieve4eI%.o: gnfs-lasieve4e.c if.h primgen32.h asm/32bit.h redu2.h \
	recurrence6.h fbgen.h real-poly-aux.h gmp-aux.h asm/medsched.h \
	asm/siever-config.h lasieve-prepn.h input-poly.h asm/lasched.h \
	strategy.h
	$(CC) $(CFLAGS) -c -DI_bits=$* -o $@ $<

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


# --- Library rules -----------------------------------------------------------

libgmp-aux.a: gmp-aux.o mpz-ull.o
	$(AR) rcs $@ $^

mpqs.o: mpqs.c asm/mpqs-config.h asm/siever-config.h
	$(CC) $(CFLAGS) -c -o $@ $<

mpqs3.o: mpqs3.c asm/mpqs-config.h asm/siever-config.h
	$(CC) $(CFLAGS) -c -o $@ $<


# --- Executable link rules ---------------------------------------------------

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


# --- Test / stat targets -----------------------------------------------------

mpqsz.o: mpqs.c asm/mpqs-config.h
	$(CC) $(CFLAGS) -O2 -DMPQS_STAT -DMPQS_ZEIT -c -o $@ $<

mpqstest: mpqstest.o mpqsz.o if.o mpz-ull.o asm/liblasieve.a
	$(CC) $(CFLAGS) -o $@ $^ $(GMP_LIB) -lm

mpqst: mpqst.o mpqsz.o if.o asm/liblasieve.a
	$(CC) $(CFLAGS) -o $@ $^ $(GMP_LIB) -lm

mpqszt.o: mpqs.c asm/mpqs-config.h
	$(CC) $(CFLAGS) -O2 -DTOTAL_STAT -DMPQS_ZEIT -c -o $@ $<

tmpqs: mpqst.o mpqszt.o if.o asm/liblasieve.a
	$(CC) $(CFLAGS) -o $@ $^ $(GMP_LIB) -lm

mpqs3z.o: mpqs3.c asm/mpqs-config.h
	$(CC) $(CFLAGS) -O2 -DMPQS3_STAT -DMPQS3_ZEIT -c -o $@ $<

mpqs3test: mpqs3test.o mpqs3z.o if.o libgmp-aux.a asm/liblasieve.a
	$(CC) $(CFLAGS) -o $@ $^ $(GMP_LIB) -lm

mpqs3z2.o: mpqs3.c asm/mpqs-config2.h
	$(CC) $(CFLAGS) -O2 -DMPQS3_STAT -DMPQS3_ZEIT -DVARIANT2 -c -o $@ $<

mpqs3test2: mpqs3test.o mpqs3z2.o if.o libgmp-aux.a asm/liblasieve.a
	$(CC) $(CFLAGS) -o $@ $^ $(GMP_LIB) -lm

mpqsstat: mpqsstat.o mpqs.o mpqs3.o if.o mpz-ull.o libgmp-aux.a asm/liblasieve.a
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


# --- Clean -------------------------------------------------------------------

clean:
	-rm -f *.o *.a asm/*.o asm/*.a asm/*.s asm/*.S gnfs-lasieve4I1*
