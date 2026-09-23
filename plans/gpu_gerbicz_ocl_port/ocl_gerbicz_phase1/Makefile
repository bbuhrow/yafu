# Makefile for the OCL_GERBICZ Phase 1 (OpenCL foundation) deliverable.
#
# This builds ONLY the new Phase 1 files. It does not touch or rebuild
# anything else in the YAFU tree.
#
# PROJDIR must point at a directory containing ocl_xface.h and
# collision_bucket.h. Defaults to ./patched -- a copy of those two
# project headers with ocl_xface.h.patch already applied, shipped ONLY
# so this Makefile builds standalone out of the box. Project sources
# are not modified in place: point PROJDIR at the real project tree
# instead once ocl_xface.h.patch has been applied there, e.g.:
#   patch -p1 -d /path/to/yafu/gnfs/poly/stage1 < ocl_xface.h.patch
#   make PROJDIR=/path/to/yafu/gnfs/poly/stage1 smoke
PROJDIR   ?= patched
STD       ?= c99
CFLAGS    := -std=$(STD) -Wall -Wextra -pedantic -DHAVE_OCL_POLY -I. -I$(PROJDIR)
LDLIBS    := -lOpenCL

SRCS      := ocl_shared.c ocl_gerbicz_ctx.c ocl_xface_stub.c test_ocl_foundation.c
OBJS      := $(SRCS:.c=.o)
BIN       := test_ocl_foundation

.PHONY: all smoke clean mingw-check

all: $(BIN)

$(BIN): $(SRCS)
	$(CC) $(CFLAGS) $(SRCS) $(LDLIBS) -o $@

smoke: $(BIN)
	./$(BIN)

# Compile-only cross-check under mingw-w64 (a GNU toolchain targeting
# Windows). This is NOT a substitute for a real MSVC (cl.exe) build --
# no MSVC toolchain was available in the Phase 1 sandbox (see STATUS).
# It catches gross portability breaks (missing headers, GNU-only
# extensions) but not MSVC-specific ones (pre-C99 declaration rules on
# old MSVC, snprintf/sscanf conformance quirks, etc.).
# CL_HEADERS_DIR must point at a directory containing a CL/cl.h (the
# portable Khronos headers, e.g. from the opencl-c-headers package or
# the Windows OpenCL SDK) -- do NOT point this at /usr/include directly,
# it will pull in glibc's own stdio.h ahead of mingw's.
mingw-check:
	x86_64-w64-mingw32-gcc -std=$(STD) -Wall -Wextra -pedantic -DHAVE_OCL_POLY \
		-I. -I$(PROJDIR) -I$(CL_HEADERS_DIR) -c ocl_shared.c ocl_gerbicz_ctx.c
	rm -f ocl_shared.o ocl_gerbicz_ctx.o

clean:
	rm -f $(BIN) *.o gerbicz_smoke_*.bin
