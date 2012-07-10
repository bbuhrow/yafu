#!/bin/sh
make clean
make x86_64 NFS=1 STATIC=1
mv yafu yafu-linux64
make clean


