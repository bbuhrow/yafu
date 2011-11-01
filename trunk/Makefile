# This source distribution is placed in the public domain by its author,
# Ben Buhrow. You may use it for any purpose, free of charge,
# without having to notify anyone. I disclaim any responsibility for any
# errors.
# 
# Optionally, please be nice and tell me if you find this source to be
# useful. Again optionally, if you add to the functionality present here
# please consider making those additions public too, so that others may 
# benefit from your work.	
# 
# Some parts of the code (and also this header), included in this 
# distribution have been reused from other sources. In particular I 
# have benefitted greatly from the work of Jason Papadopoulos's msieve @ 
# www.boo.net/~jasonp, Scott Contini's mpqs implementation, and Tom St. 
# Denis Tom's Fast Math library.  Many thanks to their kind donation of 
# code to the public domain.
#        				   --bbuhrow@gmail.com 7/28/09
# ----------------------------------------------------------------------*/

CC = gcc
#CC = x86_64-w64-mingw32-gcc-4.5.1
#CFLAGS = -march=core2 -mtune=core2
CFLAGS = -g
WARN_FLAGS = -Wall #-W -Wconversion
OPT_FLAGS = -O3
INC = -I. -Iinclude
LIBS = -lecm -lgmp

# INC += -I/sppdg/scratch/buhrow/gmp-4.2.3/install/include/
INC += -I../gmp/include
# LIBS += -L/sppdg/scratch/buhrow/gmp-4.2.3/install/lib/
LIBS += -L../gmp/lib/linux

# INC += -I/sppdg/scratch/buhrow/ecm-6.2.3/install/include/
INC += -I../gmp-ecm/include
# LIBS += -L/sppdg/scratch/buhrow/ecm-6.2.3/install/lib/
LIBS += -L../gmp-ecm/lib/linux

ifeq ($(STATIC),1)
	CFLAGS += -static
endif

ifeq ($(BLOCK),64)
	CFLAGS += -DYAFU_64K
endif

ifeq ($(PROFILE),1)
	CFLAGS += -pg 
	CFLAGS += -DPROFILING
else
	OPT_FLAGS += -fomit-frame-pointer
endif

ifeq ($(OPT_DEBUG),1)
	CFLAGS += -DOPT_DEBUG
endif

ifeq ($(TIMING),1)
	CFLAGS += -DQS_TIMING
endif

ifeq ($(NFS),1)
	CFLAGS += -DUSE_NFS
	LIBS += -L../msieve/lib/linux -lmsieve
endif

#MINGW builds don't need -pthread
#LIBS += -lm 
LIBS += -lm -lpthread

ifeq ($(FORCE_MODERN),1)
	CFLAGS += -DFORCE_MODERN
endif

ifeq ($(CC),icc)
	CFLAGS += -mtune=core2 -march=core2
endif

CFLAGS += $(OPT_FLAGS) $(WARN_FLAGS) $(INC)
	
x86: CFLAGS += -m32

	
#---------------------------Msieve file lists -------------------------
MSIEVE_SRCS = \
	factor/qs/msieve/lanczos.c \
	factor/qs/msieve/lanczos_matmul0.c \
	factor/qs/msieve/lanczos_matmul1.c \
	factor/qs/msieve/lanczos_matmul2.c \
	factor/qs/msieve/lanczos_pre.c \
	factor/qs/msieve/sqrt.c \
	factor/qs/msieve/savefile.c \
	factor/qs/msieve/gf2.c
	
MSIEVE_OBJS = $(MSIEVE_SRCS:.c=.o)
	
#---------------------------YAFU file lists -------------------------
YAFU_SRCS = \
	top/driver.c \
	top/utils.c \
	top/stack.c \
	top/calc.c \
	top/test.c \
	factor/factor_common.c \
	factor/rho.c \
	factor/squfof.c \
	factor/trialdiv.c \
	factor/tune.c \
	factor/qs/filter.c \
	factor/qs/tdiv.c \
	factor/qs/tdiv_small.c \
	factor/qs/tdiv_med.c \
	factor/qs/tdiv_large.c \
	factor/qs/tdiv_scan.c \
	factor/qs/sieve.c \
	factor/qs/new_poly.c \
	factor/qs/poly_roots.c \
	factor/qs/update_poly_roots.c \
	factor/qs/siqs_test.c \
	factor/qs/siqs_aux.c \
	factor/qs/smallmpqs.c \
	factor/qs/SIQS.c \
	factor/gmp-ecm/ecm.c \
	factor/gmp-ecm/pp1.c \
	factor/gmp-ecm/pm1.c \
	factor/nfs/nfs.c \
	factor/nfs/nfs_sieving.c \
	factor/nfs/nfs_poly.c \
	factor/nfs/nfs_postproc.c \
	factor/nfs/nfs_filemanip.c \
	factor/nfs/nfs_threading.c \
	arith/arith0.c \
	arith/arith1.c \
	arith/arith2.c \
	arith/arith3.c \
	top/eratosthenes/count.c \
	top/eratosthenes/offsets.c \
	top/eratosthenes/primes.c \
	top/eratosthenes/roots.c \
	top/eratosthenes/linesieve.c \
	top/eratosthenes/soe.c \
	top/eratosthenes/tiny.c \
	top/eratosthenes/worker.c \
	top/eratosthenes/soe_util.c \
	top/eratosthenes/wrapper.c
	
YAFU_OBJS = $(YAFU_SRCS:.c=.o)

#---------------------------Header file lists -------------------------
HEAD = include/yafu.h  \
	include/qs.h  \
	include/lanczos.h  \
	include/types.h  \
	include/calc.h  \
	include/common.h  \
	include/factor.h  \
	include/soe.h  \
	include/util.h  \
	include/types.h \
	include/yafu_string.h  \
	include/arith.h  \
	include/msieve.h  \
	include/yafu_stack.h  \
	include/yafu_ecm.h \
	include/gmp_xface.h \
	include/nfs.h

#---------------------------Make Targets -------------------------

all:
	@echo "pick a target:"
	@echo "x86       32-bit Intel/AMD systems (required if gcc used)"
	@echo "x86_64    64-bit Intel/AMD systems (required if gcc used)"
	@echo "add 'BLOCK=64' to make with 64kB QS blocksize "
	@echo "add 'TIMING=1' to make with expanded QS timing info (slower) "
	@echo "add 'PROFILE=1' to make with profiling enabled (slower) "

x86: $(MSIEVE_OBJS) $(YAFU_OBJS)
	$(CC) -m32 $(CFLAGS) $(MSIEVE_OBJS) $(YAFU_OBJS) -o yafu $(LIBS)
	
	
x86_64: $(MSIEVE_OBJS) $(YAFU_OBJS)
	$(CC) $(CFLAGS) $(MSIEVE_OBJS) $(YAFU_OBJS) -o yafu $(LIBS)
	
clean:
	rm -f $(MSIEVE_OBJS) $(YAFU_OBJS) 
	
#---------------------------Build Rules -------------------------

	
%.o: %.c $(HEAD)
	$(CC) $(CFLAGS) -c -o $@ $<





