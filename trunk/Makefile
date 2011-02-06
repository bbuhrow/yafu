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
#CFLAGS = -g -march=native -mtune=native
#CFLAGS =
WARN_FLAGS = -Wall #-W -Wconversion
OPT_FLAGS = -O3
INC = -I. -Iinclude

#MINGW builds don't need -pthread
#LIBS = -lm 
LIBS = -lm -pthread

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

ifeq ($(GMPECM),1)
	CFLAGS += -DHAVE_GMP_ECM
	CFLAGS += -DHAVE_GMP
	INC += -I../gmp/include/
	INC += -I../gmp-ecm/include/
	LIBS += -L../gmp/lib/ -lgmp
	LIBS += -L../gmp-ecm/lib/ -lecm
endif

ifeq ($(NFS),1)	
	LIBS += -L../msieve/ -lmsieve

	# NFS builds require GMP
	ifneq ($(GMPECM),1)
		CFLAGS += -DHAVE_GMP
		# INC += -I/sppdg/scratch/buhrow/gmp-4.2.3/install/include/
		INC += -I../gmp/include/
		# LIBS += -L/sppdg/scratch/buhrow/gmp-4.2.3/install/lib/ -lgmp
		LIBS += -L../gmp/lib/ -lgmp
	endif

	# NOTE: if the included msieve library was build with ECM, then 
	# ECM must also be enabled here
	CFLAGS += -DHAVE_GMP_ECM
	INC += -I../gmp-ecm/include/
	# INC += -I/sppdg/scratch/buhrow/ecm-6.2.3/install/include/
	LIBS += -L../gmp-ecm/lib/ -lecm
	# LIBS += -L/sppdg/scratch/buhrow/ecm-6.2.3/install/lib/ -lecm
endif

ifeq ($(FORCE_MODERN),1)
	CFLAGS += -DFORCE_MODERN
endif

ifeq ($(CC),icc)
	CFLAGS += -mtune=core2 -march=core2
endif

CFLAGS += $(OPT_FLAGS) $(WARN_FLAGS) $(INC)
	
x86: CFLAGS += -m32

#---------------------------Tom's Fast Math file lists -------------------------
TFM_SRCS = \
	tfm/fp_mul_comba.c \
	tfm/fp_mul_comba_small_set.c \
	tfm/fp_sqr_comba.c \
	tfm/fp_sqr_comba_small_set.c \
	tfm/fp_sqr_comba_generic.c \
	tfm/fp_2expt.c \
	tfm/fp_cmp_mag.c \
	tfm/fp_mont_small.c \
	tfm/fp_montgomery_calc_normalization.c \
	tfm/fp_montgomery_reduce.c \
	tfm/fp_montgomery_setup.c \
	tfm/fp_mul_2.c \
	tfm/s_fp_sub.c
	
TFM_OBJS = $(TFM_SRCS:.c=.o)
	
#---------------------------Msieve file lists -------------------------
MSIEVE_SRCS = \
	msieve/lanczos.c \
	msieve/lanczos_matmul0.c \
	msieve/lanczos_matmul1.c \
	msieve/lanczos_matmul2.c \
	msieve/lanczos_pre.c \
	msieve/mpqs_gf2.c \
	msieve/sqrt.c \
	msieve/savefile.c 
#	msieve/squfof_jp.c
	
MSIEVE_OBJS = $(MSIEVE_SRCS:.c=.o)
	
#---------------------------YAFU file lists -------------------------
YAFU_SRCS = \
	top/driver.c \
	top/utils.c \
	top/stack.c \
	top/calc.c \
	top/test.c \
	factor/factor_common.c \
	factor/algfactor.c \
	factor/relation.c \
	factor/sieve.c \
	factor/poly.c \
	factor/siqs_test.c \
	factor/siqs_aux.c \
	factor/smallmpqs.c \
	factor/ecm.c \
	factor/pp1.c \
	factor/pm1.c \
	factor/rho.c \
	factor/squfof.c \
	factor/trialdiv.c \
	factor/MPQS.c \
	factor/gf2.c \
	factor/pQS.c \
	factor/SIQS.c \
	factor/tune.c \
	arith/arith0.c \
	arith/monty.c \
	arith/arith1.c \
	arith/arith2.c \
	arith/arith3.c \
	eratosthenes/count.c \
	eratosthenes/offsets.c \
	eratosthenes/primes.c \
	eratosthenes/roots.c \
	eratosthenes/linesieve.c \
	eratosthenes/soe.c \
	eratosthenes/tiny.c \
	eratosthenes/worker.c \
	eratosthenes/wrapper.c \
	nfs/nfs.c

YAFU_OBJS = $(YAFU_SRCS:.c=.o)

#---------------------------Header file lists -------------------------
HEAD = include/yafu.h  \
	include/qs.h  \
	include/lanczos.h  \
	include/types.h  \
	include/tfm.h  \
	include/ecm.h \
	include/calc.h  \
	include/common.h  \
	include/factor.h  \
	include/monty.h  \
	include/soe.h  \
	include/util.h  \
	include/types.h \
	include/yafu_string.h  \
	include/arith.h  \
	include/msieve.h  \
	include/yafu_stack.h  \
	include/gmp_xface.h

#---------------------------Make Targets -------------------------

all:
	@echo "pick a target:"
	@echo "x86       32-bit Intel/AMD systems (required if gcc used)"
	@echo "x86_64    64-bit Intel/AMD systems (required if gcc used)"
	@echo "add 'BLOCK=64' to make with 64kB QS blocksize "
	@echo "add 'TIMING=1' to make with expanded QS timing info (slower) "
	@echo "add 'PROFILE=1' to make with profiling enabled (slower) "

x86: $(TFM_OBJS) $(MSIEVE_OBJS) $(YAFU_OBJS)
	$(CC) -m32 $(CFLAGS) $(TFM_OBJS) $(MSIEVE_OBJS) $(YAFU_OBJS) -o yafu $(LIBS)
	
	
x86_64: $(TFM_OBJS) $(MSIEVE_OBJS) $(YAFU_OBJS)
	$(CC) $(CFLAGS) $(TFM_OBJS) $(MSIEVE_OBJS) $(YAFU_OBJS) -o yafu $(LIBS)
	
clean:
	rm -f $(TFM_OBJS) $(MSIEVE_OBJS) $(YAFU_OBJS) 
	
#---------------------------Build Rules -------------------------

	
%.o: %.c $(HEAD)
	$(CC) $(CFLAGS) -c -o $@ $<





