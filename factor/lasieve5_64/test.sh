#!/bin/sh

bin/mingw/mpqs/gnfs-lasieve4I11e.exe -c 10000 -f 650000 -v -o nfs300_650k_I11.job.mingw.mpqs.out.1 -a nfs300_650k_I11.job 
bin/mingw/mpqs/gnfs-lasieve4I12e.exe -c 5000 -f 900000 -v -o nfs330_900k_I12.job.mingw.mpqs.out.1 -a nfs330_900k_I12.job 
bin/mingw/mpqs/gnfs-lasieve4I13e.exe -c 1000 -f 1600000 -v -o nfs360_1600k_I13.job.mingw.mpqs.out.1 -a nfs360_1600k_I13.job 
bin/mingw/mpqs/gnfs-lasieve4I13e.exe -c 1000 -f 2700000 -v -o nfs390_2700k_I13.job.mingw.mpqs.out.1 -a nfs390_2700k_I13.job 
bin/mingw/mpqs/gnfs-lasieve4I13e.exe -c 1000 -f 6500000 -v -o nfs420_6500k_I13.job.mingw.mpqs.out.1 -a nfs420_6500k_I13.job 
bin/mingw/mpqs/gnfs-lasieve4I14e.exe -c 1000 -f 15300000 -v -o nfs512_15300k_I14.job.mingw.mpqs.out.1 -a nfs512_15300k_I14.job 
bin/mingw/mpqs/gnfs-lasieve4I16e.exe -c 250 -f 50000000 -v -o nfs631_50000k_I16.job.mingw.mpqs.out.1 -a nfs631_50000k_I16.job 
bin/mingw/mpqs/gnfs-lasieve4I15e.exe -c 250 -f 100000000 -v -o snfs_100000k_I15.job.mingw.mpqs.out.1 -r snfs_100000k_I15.txt
bin/mingw/mpqs/gnfs-lasieve4I16e.exe -c 200 -f 250000000 -v -o R1340L_poly.job.mingw.mpqs.out.1 -a R1340L_poly.txt
bin/mingw/mpqs/gnfs-lasieve4I16e.exe -c 200 -f 325000000 -v -o R942_poly.job.mingw.mpqs.out.1 -r R942_poly.txt

bin/mingw/gnfs-lasieve4I11e.exe -c 10000 -f 650000 -v -o nfs300_650k_I11.job.mingw.out.1 -a nfs300_650k_I11.job 
bin/mingw/gnfs-lasieve4I12e.exe -c 5000 -f 900000 -v -o nfs330_900k_I12.job.mingw.out.1 -a nfs330_900k_I12.job 
bin/mingw/gnfs-lasieve4I13e.exe -c 1000 -f 1600000 -v -o nfs360_1600k_I13.job.mingw.out.1 -a nfs360_1600k_I13.job 
bin/mingw/gnfs-lasieve4I13e.exe -c 1000 -f 2700000 -v -o nfs390_2700k_I13.job.mingw.out.1 -a nfs390_2700k_I13.job 
bin/mingw/gnfs-lasieve4I13e.exe -c 1000 -f 6500000 -v -o nfs420_6500k_I13.job.mingw.out.1 -a nfs420_6500k_I13.job 
bin/mingw/gnfs-lasieve4I14e.exe -c 1000 -f 15300000 -v -o nfs512_15300k_I14.job.mingw.out.1 -a nfs512_15300k_I14.job 
bin/mingw/gnfs-lasieve4I16e.exe -c 250 -f 50000000 -v -o nfs631_50000k_I16.job.mingw.out.1 -a nfs631_50000k_I16.job 
bin/mingw/gnfs-lasieve4I15e.exe -c 250 -f 100000000 -v -o snfs_100000k_I15.job.mingw.out.1 -r snfs_100000k_I15.txt
bin/mingw/gnfs-lasieve4I16e.exe -c 200 -f 250000000 -v -o R1340L_poly.job.mingw.out.1 -a R1340L_poly.txt
bin/mingw/gnfs-lasieve4I16e.exe -c 200 -f 325000000 -v -o R942_poly.job.mingw.out.1 -r R942_poly.txt

bin/mingw/avx512/gnfs-lasieve4I11e.exe -c 10000 -f 650000 -v -o nfs300_650k_I11.job.mingw.avx512.out.1 -a nfs300_650k_I11.job 
bin/mingw/avx512/gnfs-lasieve4I12e.exe -c 5000 -f 900000 -v -o nfs330_900k_I12.job.mingw.avx512.out.1 -a nfs330_900k_I12.job 
bin/mingw/avx512/gnfs-lasieve4I13e.exe -c 1000 -f 1600000 -v -o nfs360_1600k_I13.job.mingw.avx512.out.1 -a nfs360_1600k_I13.job 
bin/mingw/avx512/gnfs-lasieve4I13e.exe -c 1000 -f 2700000 -v -o nfs390_2700k_I13.job.mingw.avx512.out.1 -a nfs390_2700k_I13.job 
bin/mingw/avx512/gnfs-lasieve4I13e.exe -c 1000 -f 6500000 -v -o nfs420_6500k_I13.job.mingw.avx512.out.1 -a nfs420_6500k_I13.job 
bin/mingw/avx512/gnfs-lasieve4I14e.exe -c 1000 -f 15300000 -v -o nfs512_15300k_I14.job.mingw.avx512.out.1 -a nfs512_15300k_I14.job 
bin/mingw/avx512/gnfs-lasieve4I16e.exe -c 250 -f 50000000 -v -o nfs631_50000k_I16.job.mingw.avx512.out.1 -a nfs631_50000k_I16.job 
bin/mingw/avx512/gnfs-lasieve4I15e.exe -c 250 -f 100000000 -v -o snfs_100000k_I15.job.mingw.avx512.out.1 -r snfs_100000k_I15.txt
bin/mingw/avx512/gnfs-lasieve4I16e.exe -c 200 -f 250000000 -v -o R1340L_poly.job.mingw.avx512.out.1 -a R1340L_poly.txt
bin/mingw/avx512/gnfs-lasieve4I16e.exe -c 200 -f 325000000 -v -o R942_poly.job.mingw.avx512.out.1 -r R942_poly.txt