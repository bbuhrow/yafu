# AVX2 vector Montgomery experiments (handoff 4.3).
#   make [GMP_PREFIX=...] [ARCH=-march=haswell]
#   ./a2m_bench -t          correctness only
#   ./a2m_bench -b          benchmarks (vs mpz_powm)
# ARCH defaults to AVX2 + FMA only, so the AVX-512 machines here measure
# the AVX2 code (no AVX-512 instructions are generated).

GMP_PREFIX ?= ../../../gmp-install/6.2.1
CC ?= gcc
ARCH ?= -mavx2 -mfma -mbmi2 -mtune=native
CFLAGS ?= -O3 -g -Wall
CPPFLAGS += -I$(GMP_PREFIX)/include
LDLIBS = $(GMP_PREFIX)/lib/libgmp.a -lm

OBJS = a2m_fp52.o a2m_int.o a2m_bench.o

all: a2m_bench

a2m_bench: $(OBJS)
	$(CC) $(ARCH) $(CFLAGS) -o $@ $^ $(LDLIBS)

# the FP52 kernel changes MXCSR: keep the FP ops between the writes
a2m_fp52.o: a2m_fp52.c a2m.h
	$(CC) $(CPPFLAGS) $(ARCH) $(CFLAGS) -frounding-math -c -o $@ $<

%.o: %.c a2m.h
	$(CC) $(CPPFLAGS) $(ARCH) $(CFLAGS) -c -o $@ $<

clean:
	rm -f *.o a2m_bench

.PHONY: all clean
