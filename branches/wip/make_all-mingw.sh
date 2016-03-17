#!/bin/sh
make clean
make -f Makefile.mingw x86_64 NFS=1
make clean
#make x86
#mv yafu yafu-32k-linux32
#make clean
#make x86 BLOCK=64
#mv yafu yafu-64k-linux32
#make clean




