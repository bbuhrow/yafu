# Standalone build of the RNS modmul library, tests and benchmark.
#   make GMP_PREFIX=/path/to/gmp [ARCH=-march=icelake-client]
# The AVX-512 kernels are built when ARCH enables avx512f + avx512ifma.
# "make bench" also builds AVX-ECM's vecarith52 cores (old fixed/BFIPS and the
# constant-S template cores, IFMA) for comparison.

GMP_PREFIX ?= ../../../gmp-install/6.2.1
CC ?= gcc
ARCH ?= -march=native
CFLAGS ?= -O3 -g -Wall
CPPFLAGS += -I$(GMP_PREFIX)/include
LDLIBS = $(GMP_PREFIX)/lib/libgmp.a -lm

OBJS = rns_ref.o rns_avx512.o
YOBJS = rns_yao_ref.o rns_yao_avx512.o rns_yaob_avx512.o

AE_DIR = ../../factor/avx-ecm
VM_DIR = ../vecmont
AE_SRCS = vec_common vecarith52 vecarith52_common vecarith52_special
VM_SRCS = vecarith52_tmpl
AE_OBJS = $(addprefix ae_,$(addsuffix .o,$(AE_SRCS) $(VM_SRCS)))
AE_FLAGS = -DUSE_AVX512F -DUSE_AVX512BW -DSKYLAKEX -DUSE_AVX512IFMA -DIFMA \
	-I../../include -I../../ytools -I$(AE_DIR)

all: rns_test rns_yao_test

bench: rns_bench

rns_test: rns_test.o $(OBJS)
	$(CC) $(ARCH) $(CFLAGS) -o $@ $^ $(LDLIBS)

rns_yao_test: rns_yao_test.o $(YOBJS) $(OBJS)
	$(CC) $(ARCH) $(CFLAGS) -o $@ $^ $(LDLIBS)

rns_bench: rns_bench.o $(YOBJS) $(OBJS) $(AE_OBJS)
	$(CC) $(ARCH) $(CFLAGS) -o $@ $^ $(LDLIBS)

rns_bench.o: rns_bench.c rns.h rns_yao.h
	$(CC) $(CPPFLAGS) $(AE_FLAGS) $(ARCH) $(CFLAGS) -Wno-comment -c -o $@ $<

ae_%.o: $(AE_DIR)/%.c
	$(CC) $(CPPFLAGS) $(AE_FLAGS) $(ARCH) $(CFLAGS) -w -c -o $@ $<

ae_%.o: $(VM_DIR)/%.c $(VM_DIR)/vecmont.h
	$(CC) $(CPPFLAGS) $(AE_FLAGS) $(ARCH) $(CFLAGS) -w -c -o $@ $<

%.o: %.c rns.h rns_yao.h rns_yao_vec.h
	$(CC) $(CPPFLAGS) $(ARCH) $(CFLAGS) -c -o $@ $<

clean:
	rm -f *.o rns_test rns_bench rns_yao_test

.PHONY: all bench clean
