# Standalone build of vecpoly (8-lane polynomials over Z/NZ for AVX-ECM's
# FFT stage 2) and its tests.
#   make                IFMA (-march=native)
#   make NOIFMA=1       AVX-512F only: FP-emulated 52-bit products
#   make AVX2=1         AVX2 + FMA only
#   ... INTNTT=1        integer NTT instead of the native FP one (A/B)
#   ./vp_test [bits ...] [-l logL] [-b] [-k] [-s school_max]

GMP_PREFIX ?= ../../../gmp-install/6.2.1
CC ?= gcc
ifeq ($(AVX2),1)
ARCH = -mavx2 -mfma -mbmi2 -mtune=native
ISA = -DUSE_AVX2 -DUSE_SSE41 -DUSE_BMI2
else ifeq ($(NOIFMA),1)
ARCH = -march=skylake-avx512
ISA = -DUSE_AVX512F -DUSE_AVX512BW -DSKYLAKEX
else
ARCH ?= -march=native
ISA = -DUSE_AVX512F -DUSE_AVX512BW -DSKYLAKEX -DUSE_AVX512IFMA -DIFMA
endif
# INTNTT=1: integer NTT (emulated 52-bit products) on non-IFMA targets too
ifeq ($(INTNTT),1)
ISA += -DVP_INT_NTT
endif
CFLAGS ?= -O3 -g -Wall
CPPFLAGS += -I$(GMP_PREFIX)/include $(ISA)
LDLIBS = $(GMP_PREFIX)/lib/libgmp.a -lm

HDRS = vecpoly.h ../rns/rns.h ../rns/rns_vec52.h ../vecmont/vecmont.h
OBJS = vp_setup.o vp_kern.o vp_poly.o vm_tmpl.o

all: vp_test

vp_test: vp_test.o $(OBJS)
	$(CC) $(ARCH) $(CFLAGS) -o $@ $^ $(LDLIBS)

# FP kernels may run with MXCSR set to round toward zero (AVX2)
vp_kern.o: vp_kern.c $(HDRS)
	$(CC) $(CPPFLAGS) $(ARCH) $(CFLAGS) -frounding-math -c -o $@ $<

vm_tmpl.o: ../vecmont/vecarith52_tmpl.c ../vecmont/vecmont.h
	$(CC) $(CPPFLAGS) $(ARCH) $(CFLAGS) -frounding-math -w -c -o $@ $<

%.o: %.c $(HDRS)
	$(CC) $(CPPFLAGS) $(ARCH) $(CFLAGS) -c -o $@ $<

clean:
	rm -f *.o vp_test

.PHONY: all clean
