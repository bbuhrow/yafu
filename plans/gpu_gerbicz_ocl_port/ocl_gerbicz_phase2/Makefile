# Makefile for the OCL_GERBICZ Phase 2 (primitives) deliverable.
# Builds on Phase 1's foundation (ocl_shared.*, ocl_gerbicz_ctx.*), carried
# forward into this same directory per the user's Phase 2 answer #3.
#
# PROJDIR must point at a directory containing ocl_xface.h and
# collision_bucket.h. Defaults to ./patched (see that Makefile note,
# same as Phase 1) for a standalone build; point it at the real project
# tree once ocl_xface.h.patch is applied there.
#
# Usage:
#   make smoke           # gcc, correctness tests + benchmarks
#   make CC=clang smoke
#   make mingw-check CL_HEADERS_DIR=<dir-with-only-CL/cl.h>
#   make clean

CC        ?= gcc
STD       ?= c99
PROJDIR   ?= patched
CFLAGS    := -std=$(STD) -Wall -Wextra -pedantic -DHAVE_OCL_POLY -I. -I$(PROJDIR)
LDLIBS    := -lOpenCL

SRCS      := ocl_shared.c ocl_gerbicz_ctx.c ocl_primitives.c ocl_xface_stub.c test_ocl_primitives.c
BIN       := test_ocl_primitives

.PHONY: all smoke clean mingw-check

all: $(BIN)

$(BIN): $(SRCS) ocl_primitives.h ocl_primitives_kernels.cl
	$(CC) $(CFLAGS) $(SRCS) $(LDLIBS) -o $@

smoke: $(BIN)
	./$(BIN)

# Compile-only cross-check under mingw-w64 -- see Phase 1's Makefile note
# for the caveats (not a substitute for real MSVC).
mingw-check:
	x86_64-w64-mingw32-gcc -std=$(STD) -Wall -Wextra -pedantic -DHAVE_OCL_POLY \
		-I. -I$(PROJDIR) -I$(CL_HEADERS_DIR) -c ocl_shared.c ocl_gerbicz_ctx.c ocl_primitives.c
	rm -f ocl_shared.o ocl_gerbicz_ctx.o ocl_primitives.o

clean:
	rm -f $(BIN) *.o gerbicz_smoke_*.bin gerbicz_primitives_*.bin
