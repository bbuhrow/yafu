#!/bin/sh
make clean
make -f Makefile.mingw x86_64 NFS=1
mv yafu-x64 yafu-32k-x64
make clean
make -f Makefile.mingw x86_64 BLOCK=64 NFS=1
mv yafu-x64 yafu-64k-x64
make clean
#make x86
#mv yafu yafu-32k-linux32
#make clean
#make x86 BLOCK=64
#mv yafu yafu-64k-linux32
#make clean




