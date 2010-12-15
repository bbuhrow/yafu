#!/bin/sh
make clean
make x86_64 GMPECM=1
mv yafu yafu-32k-linux64
make clean
make x86_64 GMPECM=1 BLOCK=64
mv yafu yafu-64k-linux64
make clean
make x86
mv yafu yafu-32k-linux32
make clean
make x86 BLOCK=64
mv yafu yafu-64k-linux32
make clean




