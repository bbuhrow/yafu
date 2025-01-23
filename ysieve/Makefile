# MIT License
# 
# Copyright (c) 2021 Ben Buhrow
# 
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
# 
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
# 
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

CC = gcc
CFLAGS = -g
WARN_FLAGS = -Wall # -Wconversion
OPT_FLAGS = -O3
INC = -I. -I../ytools/
LIBS = -L../ytools/
BINNAME = ysieve
OBJ_EXT = .o

# ===================== path options =============================

# we require additional search directories for gmp and ytools.
# By default, we look adjacent to the ysieve folder (i.e., ../ytools, etc.).  Change
# these if your installation locations differ.

INC += -I../../gmp-install/6.2.1/include
LIBS += -L../../gmp-install/6.2.1/lib


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


# ===================== feature options =========================
ifeq ($(PROFILE),1)
	CFLAGS += -pg 
	CFLAGS += -DPROFILING
	BINNAME := ${BINNAME:%=%_prof}
else
	OPT_FLAGS += -fomit-frame-pointer
endif

# attempt to get static builds to work... unsuccessful so far
ifeq ($(STATIC),1)
# https://software.intel.com/en-us/articles/error-ld-cannot-find-lm
	CFLAGS += -static-intel -static
#	LIBS += -Wl,-Bstatic -lm -Wl,Bdynamic -pthread
  LIBS += -L/usr/lib/x86_64-redhat-linux6E/lib64/  -lgmp -lpthread -lm -lytools
else
	LIBS +=  -lm -lgmp -lpthread -lytools
#	-ldl
endif

CFLAGS += $(OPT_FLAGS) $(WARN_FLAGS) $(INC)



#--------------------------- file list -------------------------
SRCS = \
	presieve.c \
	count.c \
	offsets.c \
	primes.c \
	roots.c \
	linesieve.c \
	soe.c \
	tiny.c \
	worker.c \
	soe_util.c \
	wrapper.c \
	demo/calc.c \
    demo/cmdOptions.c \
    demo/ysieve.c

OBJS = $(SRCS:.c=$(OBJ_EXT))

#---------------------------Header file lists -------------------------
HEAD = calc.h  \
	soe.h  \
    soe_impl.h \
	threadpool.h \
	ytools.h \
    demo/cmdOptions.h


#---------------------------Make Targets -------------------------

all: $(OBJS)
	rm -f libysieve.a
	ar r libysieve.a $(OBJS)
	ranlib libysieve.a
	$(CC) $(CFLAGS) $(OBJS) -o $(BINNAME) libysieve.a $(LIBS)


clean:
	rm -f $(OBJS)

#---------------------------Build Rules -------------------------


%$(OBJ_EXT): %.c $(HEAD)
	$(CC) $(CFLAGS) -c -o $@ $<





