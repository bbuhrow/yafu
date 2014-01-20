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
WARN_FLAGS = -Wall # -Wconversion
OPT_FLAGS = -O3
INC = -I. -Iinclude
BINNAME = yafu

# modify these for your particular gmp/gmp-ecm installation
INC += -I../gmp/include
LIBS += -L../gmp/lib/linux/x86_64

INC += -I../gmp-ecm/include/linux
LIBS += -L../gmp-ecm/lib/linux/x86_64

ifeq ($(COMPILER),icc)
	CC = icc
endif

ifeq ($(PROFILE),1)
	CFLAGS += -pg 
	CFLAGS += -DPROFILING
	BINNAME := ${BINNAME:%=%_prof}
else
	OPT_FLAGS += -fomit-frame-pointer
endif

# if this option is specified then compile *both* the sse2 and sse4.1 versions of the
# appropriate files.  The executable will then choose between them based on the runtime
# capability of the user's cpu.  In other words, sse4.1 capability is required on the
# host cpu in order to compile the fat binary, but once it is compiled it should run
# to the capability of the target user cpu.
ifeq ($(USE_SSE41),1)
	CFLAGS += -DUSE_SSE41
endif

ifeq ($(OPT_DEBUG),1)
	CFLAGS += -DOPT_DEBUG
endif

ifeq ($(TIMING),1)
	CFLAGS += -DQS_TIMING
endif

ifeq ($(NFS),1)
	CFLAGS += -DUSE_NFS
#	modify the following line for your particular msieve installation
	LIBS += -L../msieve/lib/linux/x86_64 
	LIBS += -lmsieve
endif

# modify these for your particular cuda installation
ifeq ($(CUDA),1)
	CFLAGS += -DHAVE_CUDA
#	INC += -I/users/buhrow/NVIDIA_GPU_Computing_SDK/C/common/inc
	INC += -I/usr/local/cuda/include/
	LIBS += -L/usr/lib64 -lcuda
#	LIBS += /users/buhrow/NVIDIA_GPU_Computing_SDK/C/lib/libcutil_x86_64.a
endif

ifeq ($(FORCE_MODERN),1)
	CFLAGS += -DFORCE_MODERN
endif

ifeq ($(CC),icc)
#	CFLAGS += -mtune=core2 -march=core2 
	INC += -L/usr/lib/gcc/x86_64-redhat-linux/4.4.4
endif

ifeq ($(MIC),1)
	CFLAGS += -mmic -DTARGET_MIC
	BINNAME := ${BINNAME:%=%_mic}
	OBJ_EXT = .mo
	
	INC += -I../msieve/zlib
else
	OBJ_EXT = .o
endif

LIBS += -lecm -lgmp

# attempt to get static builds to work... unsuccessful so far
ifeq ($(STATIC),1)
	CFLAGS += -static
#	LIBS += -Wl,-Bstatic -lm -Wl,Bdynamic -pthread
	LIBS += -L/usr/lib/x86_64-redhat-linux5E/lib64/ -lpthread -lm
else
	LIBS += -lpthread -lm -ldl
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
	top/mpz_prp_prime.c \
	factor/factor_common.c \
	factor/rho.c \
	factor/squfof.c \
	factor/trialdiv.c \
	factor/tune.c \
	factor/qs/filter.c \
	factor/qs/tdiv.c \
	factor/qs/tdiv_small.c \
	factor/qs/tdiv_large.c \
	factor/qs/tdiv_scan.c \
	factor/qs/large_sieve.c \
	factor/qs/new_poly.c \
	factor/qs/siqs_test.c \
	factor/tinyqs/tinySIQS.c \
	factor/qs/siqs_aux.c \
	factor/qs/smallmpqs.c \
	factor/qs/SIQS.c \
	factor/gmp-ecm/ecm.c \
	factor/gmp-ecm/pp1.c \
	factor/gmp-ecm/pm1.c \
	factor/nfs/nfs.c \
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
	
ifeq ($(MIC),1)
# just add 32k versions of these files
	YAFU_SRCS += factor/qs/tdiv_med_32k.c \
		factor/qs/tdiv_resieve_32k.c \
		factor/qs/med_sieve_32k.c \
		factor/qs/poly_roots_32k.c \
		factor/qs/update_poly_roots_32k.c
else
# add both versions
	YAFU_SRCS += factor/qs/tdiv_med_32k.c \
		factor/qs/tdiv_med_64k.c \
		factor/qs/tdiv_resieve_32k.c \
		factor/qs/tdiv_resieve_64k.c \
		factor/qs/med_sieve_32k.c \
		factor/qs/med_sieve_64k.c \
		factor/qs/poly_roots_32k.c \
		factor/qs/poly_roots_64k.c \
		factor/qs/update_poly_roots_32k.c \
		factor/qs/update_poly_roots_64k.c
endif
		
ifeq ($(USE_SSE41),1)
# these files require SSE4.1 to compile
	YAFU_SRCS += factor/qs/update_poly_roots_32k_sse4.1.c
	YAFU_SRCS += factor/qs/med_sieve_32k_sse4.1.c
endif

YAFU_OBJS = $(YAFU_SRCS:.c=$(OBJ_EXT))

#---------------------------YAFU NFS file lists -----------------------
ifeq ($(NFS),1)

YAFU_NFS_SRCS = \
	factor/nfs/nfs_sieving.c \
	factor/nfs/nfs_poly.c \
	factor/nfs/nfs_postproc.c \
	factor/nfs/nfs_filemanip.c \
	factor/nfs/nfs_threading.c \
	factor/nfs/snfs.c

YAFU_NFS_OBJS = $(YAFU_NFS_SRCS:.c=$(OBJ_EXT))

else

YAFU_NFS_OBJS =

endif

#---------------------------Header file lists -------------------------
HEAD = include/yafu.h  \
	include/qs.h  \
	factor/qs/poly_macros_32k.h \
	factor/qs/poly_macros_64k.h \
	factor/qs/poly_macros_common.h \
	factor/qs/sieve_macros_32k.h \
	factor/qs/sieve_macros_64k.h \
	factor/qs/tdiv_macros_32k.h \
	factor/qs/tdiv_macros_64k.h \
	factor/qs/tdiv_macros_common.h \
	include/lanczos.h  \
	include/types.h  \
	include/calc.h  \
	include/common.h  \
	include/factor.h  \
	include/soe.h  \
	include/util.h  \
	include/types.h \
	include/yafu_string.h  \
	include/mpz_prp_prime.h \
	include/jacobi_sum.h \
	include/arith.h  \
	include/msieve.h  \
	include/yafu_stack.h  \
	include/yafu_ecm.h \
	include/gmp_xface.h \
	include/nfs.h

ifeq ($(USE_SSE41),1)
# these files require SSE4.1 to compile
	HEAD += factor/qs/poly_macros_common_sse4.1.h
	HEAD += factor/qs/sieve_macros_32k_sse4.1.h
endif

#---------------------------Make Targets -------------------------

all:
	@echo "pick a target:"
	@echo "x86       32-bit Intel/AMD systems (required if gcc used)"
	@echo "x86_64    64-bit Intel/AMD systems (required if gcc used)"
	@echo "add 'TIMING=1' to make with expanded QS timing info (slower) "
	@echo "add 'PROFILE=1' to make with profiling enabled (slower) "

x86: $(MSIEVE_OBJS) $(YAFU_OBJS) $(YAFU_NFS_OBJS)
	$(CC) -m32 $(CFLAGS) $(MSIEVE_OBJS) $(YAFU_OBJS) $(YAFU_NFS_OBJS) -o $(BINNAME) $(LIBS)

x86_64: $(MSIEVE_OBJS) $(YAFU_OBJS) $(YAFU_NFS_OBJS)
	$(CC) $(CFLAGS) $(MSIEVE_OBJS) $(YAFU_OBJS) $(YAFU_NFS_OBJS) -o $(BINNAME) $(LIBS)


clean:
	rm -f $(MSIEVE_OBJS) $(YAFU_OBJS) $(YAFU_NFS_OBJS)

#---------------------------Build Rules -------------------------


%$(OBJ_EXT): %.c $(HEAD)
	$(CC) $(CFLAGS) -c -o $@ $<





