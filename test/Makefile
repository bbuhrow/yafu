# YAFU modular test system -- Makefile
#
# Two builds:
#   make subset   framework + Layer 0 (primitives). Needs only your include
#                 path (for mp_platform.h). No GMP, no YAFU objects. This is
#                 the quick "does the harness run" build.  [default]
#   make full     framework + Layer 0 + Layer 1. Links GMP and the YAFU
#                 objects that provide arith/monty/aprcl/tinyprp.
#
# Point the variables below at your tree, then `make full`. Everything is
# overridable on the command line, e.g.:  make full YAFU_ROOT=/path/to/yafu
# ---------------------------------------------------------------------------

CC      ?= gcc
CSTD    ?= -std=c11
OPT     ?= -O2
WARN    ?= -Wall -Wextra
CFLAGS  ?= $(CSTD) $(OPT) $(WARN)

# Root of your YAFU checkout (this test/ dir is assumed to live under it).
YAFU_ROOT ?= /sppdg/scratch/buhrow/projects/factoring/yafu.git

# Include paths. Layer 0 needs mp_platform.h; adjust if it lives elsewhere.
INCLUDES  ?= -I. \
             -I$(YAFU_ROOT)/include \
             -I$(YAFU_ROOT)/top/aprcl \
             -I$(YAFU_ROOT)/ysieve

# GMP (Layer 1 only).
GMP_INC   ?= -I$(YAFU_ROOT)/../gmp-install/6.2.0-aocc/include 
GMP_LIBS  ?= -I$(YAFU_ROOT)/../gmp-install/6.2.0-aocc/lib -lgmp

# YAFU objects to link for Layer 1. The simplest option is to build YAFU
# normally and point this at every object except the one containing main().
# Adjust to match your build (object dir, or a static lib if you make one).
#   e.g.  YAFU_OBJS = $(filter-out %/driver.o, $(wildcard $(YAFU_ROOT)/*.o))
YAFU_OBJS ?=

# ---------------------------------------------------------------------------
FRAMEWORK = testkit.c test_data.c test_main.c
LAYER0    = layer0/test_mp_arith.c layer0/test_mp_bitscan.c
LAYER1    = layer1/test_sp_arith.c layer1/test_modular.c layer1/test_primality.c

BIN       ?= yafu_test

.PHONY: all subset full run list clean

all: subset

subset:
	$(CC) $(CFLAGS) -DTK_NO_LAYER1 $(INCLUDES) \
	    $(FRAMEWORK) $(LAYER0) -o $(BIN)

full:
	@if [ -z "$(YAFU_OBJS)" ]; then \
	  echo "error: set YAFU_OBJS to the YAFU objects to link (see Makefile header)"; \
	  exit 1; \
	fi
	$(CC) $(CFLAGS) $(INCLUDES) $(GMP_INC) \
	    $(FRAMEWORK) $(LAYER0) $(LAYER1) \
	    $(YAFU_OBJS) $(GMP_LIBS) -o $(BIN)

run: $(BIN)
	./$(BIN)

list: $(BIN)
	./$(BIN) --list

clean:
	rm -f $(BIN)
