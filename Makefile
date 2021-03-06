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

CC = gcc-7.3.0
#CC = x86_64-w64-mingw32-gcc-4.5.1
#CFLAGS = -march=core2 -mtune=core2
CFLAGS = -g
WARN_FLAGS = -Wall # -Wconversion
OPT_FLAGS = -O2
INC = -I. -Iinclude -Itop/aprcl -Itop/cmdParser -Itop/ -I../msieve/zlib -I../../ysieve.git/trunk -I../../ytools.git/trunk
LIBS = -L../../ysieve.git/trunk -L../../ytools.git/trunk
BINNAME = yafu


# ===================== compiler options =========================
ifeq ($(COMPILER),icc)
	CC = icc
	INC += -L/usr/lib/gcc/x86_64-redhat-linux/4.4.4
	CFLAGS += -qopt-report=5
endif

ifeq ($(COMPILER),gcc)
	CC = gcc
endif

# ===================== architecture options =========================
# if this option is specified then compile *both* the sse2 and sse4.1 versions of the
# appropriate files.  The executable will then choose between them based on the runtime
# capability of the user's cpu.  In other words, sse4.1 capability is required on the
# host cpu in order to compile the fat binary, but once it is compiled it should run
# to the capability of the target user cpu.
ifeq ($(SKYLAKEX),1)
	CFLAGS += -DUSE_BMI2 -DUSE_AVX2 -DUSE_AVX512F -DUSE_AVX512BW -march=skylake-avx512 
endif

ifeq ($(SMALLINT),1)
	CFLAGS += -DSMALL_SIQS_INTERVALS
endif

ifeq ($(USE_BMI2),1)
# -mbmi enables _blsr_u64 and -mbmi2 enables _pdep_u64 when using gcc
  CFLAGS += -mbmi2 -mbmi -DUSE_BMI2
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

ifeq ($(NO_ZLIB),1)
  CFLAGS += -DNO_ZLIB
endif

ifeq ($(USE_SSE41),1)
	CFLAGS += -DUSE_SSE41 -m64 -msse4.1
endif

ifeq ($(KNL),1)
    ifneq ($(USE_AVX2),1)
        CFLAGS += -DUSE_AVX2 -DUSE_SSE41 
    endif

    ifeq ($(COMPILER),icc)
        CFLAGS += -DTARGET_KNL -DUSE_AVX512F -DUSE_AVX512PF -DSMALL_SIQS_INTERVALS -xMIC-AVX512 
        BINNAME = yafu_knl
    else
        CFLAGS += -DTARGET_KNL -DUSE_AVX512F -DUSE_AVX512PF -DSMALL_SIQS_INTERVALS -march=knl
        BINNAME = yafu_knl_gcc
    endif
  #-openmp
endif


ifeq ($(KNC),1)
	CFLAGS += -mmic -DTARGET_KNC -vec-report3
	BINNAME := ${BINNAME:%=%_knc}
	OBJ_EXT = .mo

	INC += -I../msieve/zlib

	INC += -I../gmp/include
	LIBS += -L../gmp/lib/linux/phi

	INC += -I../gmp-ecm/include/linux
	LIBS += -L../gmp-ecm/lib/phi/lib
else
	OBJ_EXT = .o
  
    ifeq ($(SKYLAKEX),1)
        INC += -I../../gmp_install/gmp-6.2.0/include
        LIBS += -L../../gmp_install/gmp-6.2.0/lib/
        INC += -I../../ecm_install/include/
        LIBS += -L../../ecm_install/lib/
    else
        ifeq ($(KNL),1)
			OBJ_EXT = .ko
            INC += -I../../gmp_install/gmp-6.2.0-knl/include
            LIBS += -L../../gmp_install/gmp-6.2.0-knl/lib/
            INC += -I../../ecm_install_gmp620_knl/include/
            LIBS += -L../../ecm_install_gmp620_knl/lib/
        else
            # for non avx512 systems
            INC += -I../gmp/include
            LIBS += -L../gmp/lib/
            INC += -I../gmp-ecm/include/
            LIBS += -L../gmp-ecm/lib/
        endif
    endif

	
endif

# ===================== feature options =========================
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
#	modify the following line for your particular msieve installation

	ifeq ($(KNC),1)
		LIBS += -L../msieve/
	else
        ifeq ($(COMPILER),icc)
            LIBS += -L../../msieve/lib/linux
        else
            LIBS += -L../../msieve/lib/linux/gcc73/
        endif
	endif
	LIBS += -lmsieve
endif

# modify these for your particular cuda installation
ifeq ($(CUDA),1)
	INC += -I/usr/local/cuda/include/
	LIBS += -L/usr/lib64 -lcuda
endif

ifeq ($(FORCE_GENERIC),1)
	CFLAGS += -DFORCE_GENERIC
endif

ifeq ($(SKYLAKEX),1)
    LIBS += -lecm /users/buhrow/src/c/gmp_install/gmp-6.2.0/lib/libgmp.a -lytools -lysieve
else
    LIBS += -lecm -lgmp -lytools -lysieve
endif


ifeq ($(SKYLAKEX),1)
    # define KNL now for skylakex, after handling an actual command line KNL
    KNL=1
endif




# attempt to get static builds to work... unsuccessful so far
ifeq ($(STATIC),1)
# https://software.intel.com/en-us/articles/error-ld-cannot-find-lm
	CFLAGS += -static-intel -static
#	LIBS += -Wl,-Bstatic -lm -Wl,Bdynamic -pthread
  LIBS += -L/usr/lib/x86_64-redhat-linux6E/lib64/ -lpthread -lm
else
	LIBS += -lpthread -lm -ldl
endif

ifeq ($(COMPILER),icc)
  LIBS +=  -lsvml
# -L/apps/intel/parallel_studio_xe/2017_U1/compilers_and_libraries_2017.1.132/linux/compiler/lib/intel64/ 
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
	top/aprcl/mpz_aprcl.c \
	top/cmdParser/cmdOptions.c \
	top/cmdParser/calc.c \
	factor/factor_common.c \
	factor/rho.c \
	factor/squfof.c \
	factor/trialdiv.c \
	factor/tune.c \
	factor/qs/filter.c \
	factor/qs/tdiv.c \
	factor/qs/tdiv_small.c \
	factor/qs/tdiv_large.c \
	factor/qs/large_sieve.c \
	factor/qs/new_poly.c \
	factor/qs/siqs_test.c \
	factor/tinyqs/tinySIQS.c \
	factor/qs/siqs_aux.c \
	factor/qs/smallmpqs.c \
	factor/qs/SIQS.c \
	factor/qs/med_sieve_32k.c \
	factor/qs/poly_roots_32k.c \
	factor/gmp-ecm/ecm.c \
	factor/gmp-ecm/pp1.c \
	factor/gmp-ecm/pm1.c \
    factor/gmp-ecm/tinyecm.c \
    factor/gmp-ecm/microecm.c \
	factor/nfs/nfs.c \
	arith/arith.c \
	arith/monty.c \
    factor/prime_sieve.c \
    factor/batch_factor.c \
    factor/qs/cofactorize_siqs.c \
    factor/avx-ecm/avxecm.c \
    factor/avx-ecm/avx_ecm_main.c \
    factor/avx-ecm/vec_common.c \
    factor/avx-ecm/vecarith.c \
    factor/avx-ecm/vecarith52.c
	
#top/threadpool.c \
#top/queue.c \
#top/utils.c \

ifeq ($(USE_AVX2),1)

    YAFU_SRCS += factor/qs/tdiv_med_32k_avx2.c 
    YAFU_SRCS += factor/qs/update_poly_roots_32k_avx2.c
    YAFU_SRCS += factor/qs/med_sieve_32k_avx2.c
    YAFU_SRCS += factor/qs/tdiv_resieve_32k_avx2.c

endif

ifeq ($(USE_SSE41),1)

    # these files require SSE4.1 to compile
	YAFU_SRCS += factor/qs/update_poly_roots_32k_sse4.1.c
	YAFU_SRCS += factor/qs/med_sieve_32k_sse4.1.c    

endif

ifeq ($(KNC),1)

    # these files target running on KNC hardware
	YAFU_SRCS += factor/qs/update_poly_roots_32k_knc.c
    YAFU_SRCS += factor/qs/tdiv_med_32k_knc.c
	YAFU_SRCS += factor/qs/tdiv_resieve_32k_knc.c
	YAFU_SRCS += factor/qs/tdiv_scan_knc.c

else

    ifeq ($(KNL),1)

        YAFU_SRCS += factor/qs/tdiv_scan_knl.c
        YAFU_SRCS += factor/qs/update_poly_roots_32k_knl.c
        #YAFU_SRCS += factor/qs/tdiv_resieve_32k_knl.c 

    else

        YAFU_SRCS += factor/qs/tdiv_scan.c

    endif

    # won't build with KNC
    YAFU_SRCS += factor/qs/update_poly_roots_32k.c
    YAFU_SRCS += factor/qs/tdiv_med_32k.c
    YAFU_SRCS += factor/qs/tdiv_resieve_32k.c

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
	factor/qs/poly_macros_common.h \
	factor/qs/sieve_macros_32k.h \
	factor/qs/tdiv_macros_common.h \
	include/lanczos.h  \
	include/types.h  \
	include/common.h  \
	include/factor.h  \
	top/aprcl/mpz_aprcl.h \
	top/aprcl/jacobi_sum.h \
	include/arith.h  \
	include/msieve.h  \
	include/yafu_ecm.h \
	include/gmp_xface.h \
    include/monty.h \
	include/nfs.h \
    include/prime_sieve.h \
    include/batch_factor.h \
    include/cofactorize.h \
    factor/avx-ecm/avx_ecm.h \
	top/cmdParser/cmdOptions.h \
	top/cmdParser/calc.h
	
#include/util.h  \
#top/threadpool.h \

ifeq ($(USE_AVX2),1)

	HEAD += factor/qs/poly_macros_common_avx2.h
	HEAD += factor/qs/sieve_macros_32k_avx2.h

else
  ifeq ($(USE_SSE41),1)

  # these files require SSE4.1 to compile
    HEAD += factor/qs/poly_macros_common_sse4.1.h
    HEAD += factor/qs/sieve_macros_32k_sse4.1.h
    
  endif
endif

#---------------------------Make Targets -------------------------

all: $(MSIEVE_OBJS) $(YAFU_OBJS) $(YAFU_NFS_OBJS)
	$(CC) $(CFLAGS) $(MSIEVE_OBJS) $(YAFU_OBJS) $(YAFU_NFS_OBJS) -o $(BINNAME) $(LIBS)


clean:
	rm -f $(MSIEVE_OBJS) $(YAFU_OBJS) $(YAFU_NFS_OBJS)

#---------------------------Build Rules -------------------------


%$(OBJ_EXT): %.c $(HEAD)
	$(CC) $(CFLAGS) -c -o $@ $<





