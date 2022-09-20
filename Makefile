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
CFLAGS = -g -m64 -std=gnu99 -DUSE_SSE2
#CFLAGS += -march=core2 -mtune=core2
WARN_FLAGS = -Wall # -Wconversion
OPT_FLAGS = -O2

BINNAME = yafu
OBJ_EXT = .o

# ===================== path options =============================

# standard search directories for headers/libraries within yafu.
# These should normally not be modified.
INC = -I. -Iinclude -Itop/aprcl -Itop/cmdParser -Itop/  -Ifactor/gmp-ecm
LIBS = -L.

# we require additional search directories for msieve, zlib, gmp, gmp-ecm, 
# ytools, and ysieve for libraries and headers.  By default, we look
# adjacent to the yafu folder (i.e., ../ysieve, ../ytools, etc.).  Change
# these if your installation locations differ.
INC += -I../ysieve -I../ytools
LIBS += -L../ysieve/ -L../ytools/

INC += -I../gmp_install/gmp-6.2.0/include
LIBS += -L../gmp_install/gmp-6.2.0/lib

INC += -I../ecm_install/include/
LIBS += -L../ecm_install/lib/

INC += -I../msieve/zlib 
LIBS += -L../msieve/

# ===================== compiler options =========================
ifeq ($(COMPILER),icc)
	CC = icc
	INC += -L/usr/lib/gcc/x86_64-redhat-linux/4.4.4
	CFLAGS += -qopt-report=5
endif

# ===================== architecture options =========================
# if this option is specified then compile *both* the sse2 and sse4.1 versions of the
# appropriate files.  The executable will then choose between them based on the runtime
# capability of the user's cpu.  In other words, sse4.1 capability is required on the
# host cpu in order to compile the fat binary, but once it is compiled it should run
# to the capability of the target user cpu.


ifeq ($(USE_SSE41),1)
	CFLAGS += -DUSE_SSE41 -msse4.1
endif

ifeq ($(USE_AVX2),1)
	USE_SSE41=1
	CFLAGS += -DUSE_AVX2 -DUSE_SSE41 

ifeq ($(COMPILER),icc)
	CFLAGS += -march=core-avx2  
else
	CFLAGS += -mavx2 
endif

endif

ifeq ($(ICELAKE),1)
	CFLAGS += -DUSE_BMI2 -DUSE_AVX2 -DUSE_AVX512F -DUSE_AVX512BW -DSKYLAKEX -DIFMA -march=icelake-client
	SKYLAKEX = 1
else

ifeq ($(SKYLAKEX),1)
	CFLAGS += -DUSE_BMI2 -DUSE_AVX2 -DUSE_AVX512F -DUSE_AVX512BW -DSKYLAKEX -march=skylake-avx512 
endif
	
endif

ifeq ($(USE_BMI2),1)
# -mbmi enables _blsr_u64 and -mbmi2 enables _pdep_u64 when using gcc
  CFLAGS += -mbmi2 -mbmi -DUSE_BMI2
endif

ifeq ($(KNL),1)
ifneq ($(USE_AVX2),1)
    CFLAGS += -DUSE_AVX2 -DUSE_SSE41 
endif

CFLAGS += -DTARGET_KNL -DUSE_AVX512F -DUSE_AVX512PF -DSMALL_SIQS_INTERVALS 
BINNAME = yafu_knl

ifeq ($(COMPILER),icc)
    CFLAGS += -xMIC-AVX512 
else
    CFLAGS += -march=knl
endif

endif

# ===================== feature options =========================
ifeq ($(NO_ZLIB),1)
  CFLAGS += -DNO_ZLIB
endif

ifeq ($(SMALLINT),1)
	CFLAGS += -DSMALL_SIQS_INTERVALS
endif

ifeq ($(PROFILE),1)
	CFLAGS += -pg 
	CFLAGS += -DPROFILING
	BINNAME := ${BINNAME:%=%_prof}
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
	LIBS += -lmsieve
endif

ifeq ($(FORCE_GENERIC),1)
	CFLAGS += -DFORCE_GENERIC
endif

# make sure we get the correct libgmp linked by using an absolute path
LIBS += -lecm /users/buhrow/src/c/gmp_install/gmp-6.2.0/lib/libgmp.a -lytools -lysieve
#LIBS += -lecm -lgmp -lytools -lysieve

ifeq ($(SKYLAKEX),1)
    # define KNL now for skylakex, after handling an actual command line KNL
    KNL=1
endif

# attempt to get static builds to work... unsuccessful so far
ifeq ($(STATIC),1)
# https://software.intel.com/en-us/articles/error-ld-cannot-find-lm
	CFLAGS += -static-intel -static
	LIBS += -L/usr/lib/x86_64-redhat-linux6E/lib64/ -lpthread -lm
else
	LIBS += -lpthread -lm
endif

ifeq ($(MINGW),1)
# not needed with mingw
#	-ldl
else
	LIBS += -ldl
endif

ifeq ($(COMPILER),icc)
	LIBS +=  -lsvml
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

MSIEVE_OBJS = $(MSIEVE_SRCS:.c=$(OBJ_EXT))

#---------------------------YAFU file lists -------------------------
YAFU_SRCS = \
	top/driver.c \
	top/test.c \
	factor/tune.c \
	factor/autofactor.c \
	top/cmdParser/cmdOptions.c \
	top/cmdParser/calc.c
	
COMMON_SRCS = \
	top/aprcl/mpz_aprcl.c \
	factor/factor_common.c \
	factor/rho.c \
	factor/squfof.c \
	factor/trialdiv.c \
	arith/arith.c \
	arith/monty.c \
	factor/gmp-ecm/tinyecm.c \
    factor/gmp-ecm/microecm.c
	
ECM_SRCS = \
	factor/gmp-ecm/ecm.c \
	factor/gmp-ecm/pp1.c \
	factor/gmp-ecm/pm1.c \
    factor/avx-ecm/avxecm.c \
	factor/avx-ecm/avxppm1.c \
    factor/avx-ecm/avx_ecm_main.c \
    factor/avx-ecm/vec_common.c \
    factor/avx-ecm/vecarith.c \
    factor/avx-ecm/vecarith52.c \
	factor/avx-ecm/vecarith52_special.c \
	factor/avx-ecm/vecarith52_common.c

YAFU_SIQS_SRCS = \
	factor/qs/filter.c \
	factor/qs/tdiv.c \
	factor/qs/tdiv_small.c \
	factor/qs/tdiv_large.c \
	factor/qs/tdiv_scan.c \
	factor/qs/large_sieve.c \
	factor/qs/new_poly.c \
	factor/qs/siqs_test.c \
	factor/qs/siqs_aux.c \
	factor/qs/smallmpqs.c \
	factor/qs/SIQS.c \
	factor/qs/med_sieve_32k.c \
	factor/qs/poly_roots_32k.c \
	factor/prime_sieve.c \
    factor/batch_factor.c \
    factor/qs/cofactorize_siqs.c
	
SIQS_BIN_SRCS = factor/qs/qs_demo/siqs_demo.c \
	factor/qs/qs_demo/cmdOptions.c \
	factor/qs/qs_demo/calc.c

ECM_BIN_SRCS = factor/ecm_demo/ecm_demo.c \
	factor/ecm_demo/cmdOptions.c \
	factor/ecm_demo/calc.c

ifeq ($(USE_AVX2),1)

    YAFU_SIQS_SRCS += factor/qs/tdiv_med_32k_avx2.c 
    YAFU_SIQS_SRCS += factor/qs/update_poly_roots_32k_avx2.c
    YAFU_SIQS_SRCS += factor/qs/med_sieve_32k_avx2.c
    YAFU_SIQS_SRCS += factor/qs/tdiv_resieve_32k_avx2.c

endif

ifeq ($(USE_SSE41),1)

    # these files require SSE4.1 to compile
	YAFU_SIQS_SRCS += factor/qs/update_poly_roots_32k_sse4.1.c
	YAFU_SIQS_SRCS += factor/qs/med_sieve_32k_sse4.1.c    

endif


ifeq ($(KNL),1)

	YAFU_SIQS_SRCS += factor/qs/update_poly_roots_32k_knl.c

endif

YAFU_SIQS_SRCS += factor/qs/update_poly_roots_32k.c
YAFU_SIQS_SRCS += factor/qs/tdiv_med_32k.c
YAFU_SIQS_SRCS += factor/qs/tdiv_resieve_32k.c

YAFU_OBJS = $(YAFU_SRCS:.c=$(OBJ_EXT))
YAFU_SIQS_OBJS = $(YAFU_SIQS_SRCS:.c=$(OBJ_EXT))
YAFU_ECM_OBJS = $(ECM_SRCS:.c=$(OBJ_EXT))
YAFU_COMMON_OBJS = $(COMMON_SRCS:.c=$(OBJ_EXT))
SIQS_BIN_OBJS = $(SIQS_BIN_SRCS:.c=$(OBJ_EXT))

#---------------------------YAFU NFS file lists -----------------------
ifeq ($(NFS),1)

YAFU_NFS_SRCS = \
	factor/nfs/nfs_sieving.c \
	factor/nfs/nfs_poly.c \
	factor/nfs/nfs_postproc.c \
	factor/nfs/nfs_filemanip.c \
	factor/nfs/nfs_threading.c \
	factor/nfs/snfs.c \
	factor/nfs/nfs.c

YAFU_NFS_OBJS = $(YAFU_NFS_SRCS:.c=$(OBJ_EXT))

else

YAFU_NFS_OBJS =

endif

#---------------------------Header file lists -------------------------
COMMON_HEAD = include/gmp_xface.h \
	include/monty.h \
	include/arith.h  \
	include/common.h  \
	top/aprcl/jacobi_sum.h \
	top/aprcl/mpz_aprcl.h \
	include/factor.h
	
ECM_HEAD = include/yafu_ecm.h \
	factor/avx-ecm/avx_ecm.h \
	factor/gmp-ecm/microecm.h

SIQS_HEAD = include/qs.h  \
	include/qs_impl.h \
	factor/qs/poly_macros_32k.h \
	factor/qs/poly_macros_common.h \
	factor/qs/sieve_macros_32k.h \
	factor/qs/tdiv_macros_common.h \
	include/lanczos.h  \
	include/prime_sieve.h \
    include/batch_factor.h \
    include/cofactorize.h
	
SIQS_BIN_HEAD = factor/qs_demo/siqs_demo.h \
	factor/qs_demo/calc.h \
	factor/qs_demo/cmdOptions.h

YAFU_HEAD = include/yafu.h \
	top/cmdParser/cmdOptions.h \
	top/cmdParser/calc.h
	
NFS_HEAD = include/nfs.h \
	include/nfs_impl.h \
	include/msieve_common.h
	

ifeq ($(USE_AVX2),1)

	SIQS_HEAD += factor/qs/poly_macros_common_avx2.h

else
	ifeq ($(USE_SSE41),1)

		# these files require SSE4.1 to compile
		SIQS_HEAD += factor/qs/poly_macros_common_sse4.1.h
		SIQS_HEAD += factor/qs/sieve_macros_32k_sse4.1.h
    
	endif
endif

#---------------------------Make Targets -------------------------

yafu: $(MSIEVE_OBJS) $(YAFU_SIQS_OBJS) $(YAFU_OBJS) $(YAFU_NFS_OBJS) $(YAFU_ECM_OBJS) $(YAFU_COMMON_OBJS)
	rm -f libysiqs.a
	ar r libysiqs.a $(YAFU_SIQS_OBJS) $(YAFU_COMMON_OBJS) $(MSIEVE_OBJS) 
	ranlib libysiqs.a
	rm -f libyecm.a
	ar r libyecm.a $(YAFU_ECM_OBJS) $(YAFU_COMMON_OBJS)
	ranlib libyecm.a
	rm -f libynfs.a
	ar r libynfs.a $(YAFU_NFS_OBJS) $(YAFU_COMMON_OBJS)
	ranlib libynfs.a
	$(CC) $(CFLAGS) $(YAFU_OBJS) -o $(BINNAME) -lysiqs  -lyecm  -lynfs $(LIBS) 

siqs: $(MSIEVE_OBJS) $(YAFU_SIQS_OBJS) $(YAFU_COMMON_OBJS) $(SIQS_BIN_OBJS) 
	rm -f libysiqs.a
	ar r libysiqs.a $(YAFU_SIQS_OBJS) $(YAFU_COMMON_OBJS) $(MSIEVE_OBJS) 
	ranlib libysiqs.a
	$(CC) $(CFLAGS) $(SIQS_BIN_OBJS) -o siqs_demo -lysiqs  $(LIBS)
	
ecm: $(YAFU_ECM_OBJS) $(YAFU_COMMON_OBJS) $(ECM_BIN_OBJS) 
	rm -f libyecm.a
	ar r libyecm.a $(YAFU_ECM_OBJS) $(YAFU_COMMON_OBJS)
	ranlib libyecm.a
	$(CC) $(CFLAGS) $(ECM_BIN_OBJS) -o ecm_demo -libyecm.a  $(LIBS)

clean:
	rm -f $(MSIEVE_OBJS) $(YAFU_OBJS) $(YAFU_NFS_OBJS) $(YAFU_SIQS_OBJS) $(YAFU_ECM_OBJS) $(YAFU_COMMON_OBJS)

#---------------------------Build Rules -------------------------

#%$(OBJ_EXT): %.c $(SIQS_BIN_HEAD) $(COMMON_HEAD) $(NFS_HEAD) $(ECM_HEAD) $(YAFU_HEAD) $(SIQS_HEAD)
#	$(CC) $(CFLAGS) -c -o $@ $<
#
#

%$(OBJ_EXT): %.c $(COMMON_HEAD)
	$(CC) $(CFLAGS) -c -o $@ $<
	
%$(OBJ_EXT): %.c $(NFS_HEAD)
	$(CC) $(CFLAGS) -c -o $@ $<
	
%$(OBJ_EXT): %.c $(ECM_HEAD)
	$(CC) $(CFLAGS) -c -o $@ $<
	
%$(OBJ_EXT): %.c $(YAFU_HEAD)
	$(CC) $(CFLAGS) -c -o $@ $<
	
%$(OBJ_EXT): %.c $(SIQS_HEAD)
	$(CC) $(CFLAGS) -c -o $@ $<

%$(OBJ_EXT): %.c $(SIQS_BIN_HEAD)
	$(CC) $(CFLAGS) -c -o $@ $<
	

