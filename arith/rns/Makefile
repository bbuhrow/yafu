# Standalone build of the RNS modmul library, tests and benchmark.
#   make GMP_PREFIX=/path/to/gmp [ARCH=-march=icelake-client]
#   make NOIFMA=1   AVX-512F only: FP-emulated 52-bit products (Skylake-X);
#                   runs on IFMA machines too
#   make AVX2=1     AVX2 + FMA only (no AVX-512 instructions): the batch Yao
#                   kernel (rnsyb_*_avx2) and the references; bench compares
#                   it with the AVX2 template cores
# The AVX-512 kernels are built when ARCH enables avx512f (+ avx512ifma for
# the IFMA products).  "make bench" also builds AVX-ECM's vecarith52 cores
# (old fixed/BFIPS and the constant-S template cores) for comparison.

GMP_PREFIX ?= ../../../gmp-install/6.2.1
CC ?= gcc
ifeq ($(AVX2),1)
ARCH = -mavx2 -mfma -mbmi2 -mtune=native
AE_ISA = -DUSE_AVX2 -DUSE_SSE41 -DUSE_BMI2
else ifeq ($(NOIFMA),1)
ARCH = -march=skylake-avx512
AE_ISA = -DUSE_AVX512F -DUSE_AVX512BW -DSKYLAKEX
else
ARCH ?= -march=native
AE_ISA = -DUSE_AVX512F -DUSE_AVX512BW -DSKYLAKEX -DUSE_AVX512IFMA -DIFMA
endif
CFLAGS ?= -O3 -g -Wall
CPPFLAGS += -I$(GMP_PREFIX)/include
LDLIBS = $(GMP_PREFIX)/lib/libgmp.a -lm

OBJS = rns_ref.o rns_avx512.o
YOBJS = rns_yao_ref.o rns_yao_avx512.o rns_yaob_avx512.o rns_yaob_avx2.o

AE_DIR = ../../factor/avx-ecm
VM_DIR = ../vecmont
AE_SRCS = vec_common vecarith52 vecarith52_oldcores vecarith52_common vecarith52_special
VM_SRCS = vecarith52_tmpl
AE_OBJS = $(addprefix ae_,$(addsuffix .o,$(AE_SRCS) $(VM_SRCS)))
AE_FLAGS = $(AE_ISA) -I../../include -I../../ytools -I$(AE_DIR)

all: rns_test rns_yao_test

bench: rns_bench

rns_test: rns_test.o $(OBJS)
	$(CC) $(ARCH) $(CFLAGS) -o $@ $^ $(LDLIBS)

rns_yao_test: rns_yao_test.o $(YOBJS) $(OBJS)
	$(CC) $(ARCH) $(CFLAGS) -o $@ $^ $(LDLIBS)

rns_bench: rns_bench.o $(YOBJS) $(OBJS) $(AE_OBJS)
	$(CC) $(ARCH) $(CFLAGS) -o $@ $^ $(LDLIBS)

# above 16400 bits (handoff 4.7): karatsuba REDC vs template vs large-s Yao (AVX-512)
big_bench: big_bench.o $(YOBJS) $(OBJS) $(AE_OBJS)
	$(CC) $(ARCH) $(CFLAGS) -o $@ $^ $(LDLIBS)

big_bench.o: big_bench.c rns.h rns_yao.h
	$(CC) $(CPPFLAGS) $(AE_FLAGS) $(ARCH) $(CFLAGS) -Wno-comment -c -o $@ $<

rns_bench.o: rns_bench.c rns.h rns_yao.h rns_vec52.h
	$(CC) $(CPPFLAGS) $(AE_FLAGS) $(ARCH) $(CFLAGS) -Wno-comment -c -o $@ $<

ae_%.o: $(AE_DIR)/%.c
	$(CC) $(CPPFLAGS) $(AE_FLAGS) $(ARCH) $(CFLAGS) -w -c -o $@ $<

ae_%.o: $(VM_DIR)/%.c $(VM_DIR)/vecmont.h
	$(CC) $(CPPFLAGS) $(AE_FLAGS) $(ARCH) $(CFLAGS) -frounding-math -w -c -o $@ $<

HDRS = rns.h rns_yao.h rns_yao_vec.h rns_vec52.h rns_yaob_kern.h

# the AVX2 kernel switches MXCSR to round toward zero around each call
rns_yaob_avx2.o: rns_yaob_avx2.c $(HDRS)
	$(CC) $(CPPFLAGS) $(ARCH) $(CFLAGS) -frounding-math -c -o $@ $<

%.o: %.c $(HDRS)
	$(CC) $(CPPFLAGS) $(ARCH) $(CFLAGS) -c -o $@ $<

clean:
	rm -f *.o rns_test rns_bench rns_yao_test big_bench

.PHONY: all bench clean
