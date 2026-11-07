# Standalone test/benchmark of AVX-ECM's vector Montgomery arithmetic.
#   make GMP_PREFIX=/path/to/gmp [ARCH=-march=icelake-client] [TMPL_FLAGS=...]
#   ./ae_bench -t          correctness only
#   ./ae_bench -b -p       benchmarks (with powm)

GMP_PREFIX ?= ../../../../gmp-install/6.2.1
CC ?= gcc
ARCH ?= -march=native
CFLAGS ?= -O3 -g
CPPFLAGS += -I$(GMP_PREFIX)/include
LDLIBS = $(GMP_PREFIX)/lib/libgmp.a -lm

AE_DIR = ..
VM_DIR = ../../../arith/vecmont
AE_SRCS = vec_common vecarith52 vecarith52_common vecarith52_special
AE_OBJS = $(addprefix ae_,$(addsuffix .o,$(AE_SRCS)))
# make NOIFMA=1: AVX-512F only (FP-emulated products), e.g. for Skylake-X;
# runs on IFMA machines too
ifeq ($(NOIFMA),1)
AE_ISA = -DUSE_AVX512F -DUSE_AVX512BW -DSKYLAKEX
ARCH = -march=skylake-avx512
else
AE_ISA = -DUSE_AVX512F -DUSE_AVX512BW -DSKYLAKEX -DUSE_AVX512IFMA -DIFMA
endif
AE_FLAGS = $(AE_ISA) -I../../../include -I../../../ytools -I$(AE_DIR)
TMPL_FLAGS ?=

all: ae_bench

ae_bench: ae_bench.o ae_tmpl.o $(AE_OBJS)
	$(CC) $(ARCH) $(CFLAGS) -o $@ $^ $(LDLIBS)

ae_bench.o: ae_bench.c
	$(CC) $(CPPFLAGS) $(AE_FLAGS) $(ARCH) $(CFLAGS) -w -c -o $@ $<

ae_tmpl.o: $(VM_DIR)/vecarith52_tmpl.c $(VM_DIR)/vecmont.h
	$(CC) $(CPPFLAGS) $(AE_FLAGS) $(TMPL_FLAGS) $(ARCH) $(CFLAGS) -Wall -Wno-comment -c -o $@ $<

ae_%.o: $(AE_DIR)/%.c
	$(CC) $(CPPFLAGS) $(AE_FLAGS) $(ARCH) $(CFLAGS) -w -c -o $@ $<

clean:
	rm -f *.o ae_bench

.PHONY: all clean
