
WIN = 0
WIN64 = 0

ifeq ($(WIN),1)
	CUDA_ROOT = $(shell echo $$CUDA_PATH)
	NVCC = "$(CUDA_ROOT)/bin/nvcc"
	EXT = dll
ifeq ($(WIN64),1)
	NVCCFLAGS += -m64
else
	NVCCFLAGS += -m32
endif
	NVCCFLAGS += -Xptxas -v -Xcudafe -\# -shared -Xcompiler /fp:strict
else
	NVCC = "$(shell which nvcc)"
	CUDA_ROOT = $(shell dirname $(NVCC))/../
	EXT = so
	NVCCFLAGS += -std=c++11 -Xptxas -v -Xcudafe -\# -shared -Xcompiler -ffloat-store \
			-Xcompiler -fPIC -Xcompiler -fvisibility=hidden
endif

COMMA = ,
ifdef sm
	SM_ARCH = $(subst $(COMMA),-,$(sm))
else 
    SM_ARCH = 600
endif

ifeq (900, $(findstring 900, $(SM_ARCH)))
    SM_TARGETS  += -gencode=arch=compute_90,code=\"sm_90,compute_90\"
    SM_DEF              += -DSM900
endif
ifeq (860, $(findstring 860, $(SM_ARCH)))
    SM_TARGETS  += -gencode=arch=compute_86,code=\"sm_86,compute_86\"
    SM_DEF              += -DSM800
endif
ifeq (800, $(findstring 800, $(SM_ARCH)))
    SM_TARGETS  += -gencode=arch=compute_80,code=\"sm_80,compute_80\"
    SM_DEF              += -DSM800
endif
ifeq (750, $(findstring 750, $(SM_ARCH)))
    SM_TARGETS  += -gencode=arch=compute_75,code=\"sm_75,compute_75\"
    SM_DEF              += -DSM700
endif
ifeq (700, $(findstring 700, $(SM_ARCH)))
    SM_TARGETS  += -gencode=arch=compute_70,code=\"sm_70,compute_70\"
    SM_DEF              += -DSM700
endif
ifeq (620, $(findstring 620, $(SM_ARCH)))
    SM_TARGETS  += -gencode=arch=compute_62,code=\"sm_62,compute_62\"
    SM_DEF              += -DSM620
endif
ifeq (610, $(findstring 610, $(SM_ARCH)))
    SM_TARGETS  += -gencode=arch=compute_61,code=\"sm_61,compute_61\"
    SM_DEF              += -DSM610
endif
ifeq (600, $(findstring 600, $(SM_ARCH)))
    SM_TARGETS  += -gencode=arch=compute_60,code=\"sm_60,compute_60\"
    SM_DEF              += -DSM600
endif
ifeq (520, $(findstring 520, $(SM_ARCH)))
    SM_TARGETS 	+= -gencode=arch=compute_52,code=\"sm_52,compute_52\"
    SM_DEF 		+= -DSM520
endif
ifeq (370, $(findstring 370, $(SM_ARCH)))
    SM_TARGETS 	+= -gencode=arch=compute_37,code=\"sm_37,compute_37\"
    SM_DEF 		+= -DSM370
endif
ifeq (350, $(findstring 350, $(SM_ARCH)))
    SM_TARGETS 	+= -gencode=arch=compute_35,code=\"sm_35,compute_35\"
    SM_DEF 		+= -DSM350
endif
ifeq (300, $(findstring 300, $(SM_ARCH)))
    SM_TARGETS 	+= -gencode=arch=compute_30,code=\"sm_30,compute_30\"
    SM_DEF 		+= -DSM300
endif
ifeq (210, $(findstring 210, $(SM_ARCH)))
    SM_TARGETS 	+= -gencode=arch=compute_20,code=\"sm_21,compute_20\"
    SM_DEF 		+= -DSM210
endif
ifeq (200, $(findstring 200, $(SM_ARCH)))
    SM_TARGETS 	+= -gencode=arch=compute_20,code=\"sm_20,compute_20\"
    SM_DEF 		+= -DSM200
endif

rwildcard=$(foreach d,$(wildcard $1*),$(call rwildcard,$d/,$2) $(filter $(subst *,%,$2),$d))

CUB_DEPS = 	$(call rwildcard, cub,*.cuh)

INC = -I. -I"$(CUDA_ROOT)/include"
#INC = -I"$(CUDA_ROOT)/include" -I.

DEPS = ./Makefile \
	sort_engine.cu \
	sort_engine.h \
	spmv_engine.cu \
	spmv_engine.h \
	$(CUB_DEPS)

LIBNAME = sort_engine
LIBNAME2 = spmv_engine

all: $(LIBNAME).$(EXT) $(LIBNAME2).$(EXT)
	touch built

clean :
	rm -f  *.dll *.so *.lib *.exp built

$(LIBNAME).$(EXT) : $(DEPS)
	$(NVCC) $(SM_TARGETS) $(SM_DEF) -o $@ sort_engine.cu $(NVCCFLAGS) $(INC) -O3 -DTHRUST_IGNORE_CUB_VERSION_CHECK

$(LIBNAME2).$(EXT) : $(DEPS)
	$(NVCC) $(SM_TARGETS) $(SM_DEF) -o $@ spmv_engine.cu $(NVCCFLAGS) -DVBITS=$(VBITS) $(INC) -O3 -DTHRUST_IGNORE_CUB_VERSION_CHECK

