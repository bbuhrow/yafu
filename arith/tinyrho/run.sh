#!/bin/bash
# run rho_test from the yafu root: run.sh ISA MODE [args...]
cd /mnt/c/Users/bbuhr/projects/factoring/yafu-3.1.8
isa=$1; shift
arith/tinyrho/rho_test_$isa "$@"
