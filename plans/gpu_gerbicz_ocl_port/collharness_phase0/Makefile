# Portable C99 build for gcc or clang. For MSVC, see README.md's
# "Building on Windows / MSVC" section (a plain cl.exe invocation --
# no VLAs, no __int128, no GNU extensions are used anywhere).

CC      ?= cc
CFLAGS  ?= -std=c99 -O2 -Wall -Wextra -pedantic -Iinclude
LDFLAGS ?=

SRCS = src/main.c src/collcase_io.c src/gen.c src/cpuref_exact.c \
       src/cpuref_stats.c src/cmp.c
OBJS = $(SRCS:.c=.o)
BIN  = collharness

all: $(BIN)

$(BIN): $(OBJS)
	$(CC) $(CFLAGS) -o $@ $(OBJS) $(LDFLAGS) -lm

%.o: %.c
	$(CC) $(CFLAGS) -c -o $@ $<

clean:
	rm -f $(OBJS) $(BIN)
	rm -rf suite_out

.PHONY: all clean test

test: $(BIN)
	./run_tests.sh
