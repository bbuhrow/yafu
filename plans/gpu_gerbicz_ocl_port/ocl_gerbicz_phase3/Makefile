# Makefile for the OCL_GERBICZ Phase 3 (collision engine) deliverable.
# Builds on Phase 1 (ocl_shared.*, ocl_gerbicz_ctx.*) and Phase 2
# (ocl_primitives.*), carried forward into this directory per the
# user's earlier answer (alongside for now).
#
# Also builds and links directly against the REAL Phase 0 collharness
# source (harness/*.c) for validation (task 3.5) -- not a reimplemented
# comparator, the actual gen_generate/cpuref_exact/cpuref_stats/cmp_run.
#
# PROJDIR defaults to ./patched (Phase 1/2's guard-fixed, arg-bumped
# copy of ocl_xface.h + unmodified collision_bucket.h). Point it at the
# real project tree once ocl_xface_phase3.h.patch is applied there
# (on top of Phase 1's ocl_xface.h.patch, already applied per STATUS).
#
# Usage:
#   make smoke                    # gcc, ./patched, full validation suite
#   make CC=clang smoke
#   make PROJDIR=/real/path smoke
#   make mingw-check CL_HEADERS_DIR=<dir-with-only-CL/cl.h>
#   make clean

CC        ?= gcc
STD       ?= c99
PROJDIR   ?= patched
CFLAGS    := -std=$(STD) -Wall -Wextra -pedantic -DHAVE_OCL_POLY -I. -I$(PROJDIR) -Iharness
LDLIBS    := -lOpenCL -lm

SRCS      := ocl_shared.c ocl_gerbicz_ctx.c ocl_primitives.c ocl_collision.c \
             ocl_xface_stub.c test_ocl_collision.c
HARNESS_SRCS := harness/collcase_io.c harness/gen.c harness/cpuref_exact.c \
                harness/cpuref_stats.c harness/cmp.c
BIN       := test_ocl_collision

.PHONY: all smoke clean mingw-check

all: $(BIN)

$(BIN): $(SRCS) $(HARNESS_SRCS) ocl_collision.h ocl_collision_kernels.cl
	$(CC) $(CFLAGS) $(SRCS) $(HARNESS_SRCS) $(LDLIBS) -o $@

smoke: $(BIN)
	./$(BIN)

# Compile-only cross-check under mingw-w64 -- see Phase 1's Makefile
# note for the caveats (not a substitute for real MSVC).
mingw-check:
	x86_64-w64-mingw32-gcc -std=$(STD) -Wall -Wextra -pedantic -DHAVE_OCL_POLY \
		-I. -I$(PROJDIR) -I$(CL_HEADERS_DIR) \
		-c ocl_shared.c ocl_gerbicz_ctx.c ocl_primitives.c ocl_collision.c
	rm -f ocl_shared.o ocl_gerbicz_ctx.o ocl_primitives.o ocl_collision.o

clean:
	rm -f $(BIN) *.o gerbicz_*.bin
	rm -f /tmp/occ_*.collcase
