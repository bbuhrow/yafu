#!/bin/sh

bin/mpqs/gnfs-lasieve4I11e -c 10000 -f 650000 -v -o nfs300_650k_I11.job.wsl.mpqs.out.1 -a nfs300_650k_I11.job 
bin/mpqs/gnfs-lasieve4I12e -c 5000 -f 900000 -v -o nfs330_900k_I12.job.wsl.mpqs.out.1 -a nfs330_900k_I12.job 
bin/mpqs/gnfs-lasieve4I13e -c 1000 -f 1600000 -v -o nfs360_1600k_I13.job.wsl.mpqs.out.1 -a nfs360_1600k_I13.job 
bin/mpqs/gnfs-lasieve4I13e -c 1000 -f 2700000 -v -o nfs390_2700k_I13.job.wsl.mpqs.out.1 -a nfs390_2700k_I13.job 
bin/mpqs/gnfs-lasieve4I13e -c 1000 -f 6500000 -v -o nfs420_6500k_I13.job.wsl.mpqs.out.1 -a nfs420_6500k_I13.job 
bin/mpqs/gnfs-lasieve4I14e -c 1000 -f 15300000 -v -o nfs512_15300k_I14.job.wsl.mpqs.out.1 -a nfs512_15300k_I14.job 
bin/mpqs/gnfs-lasieve4I16e -c 250 -f 50000000 -v -o nfs631_50000k_I16.job.wsl.mpqs.out.1 -a nfs631_50000k_I16.job 
bin/mpqs/gnfs-lasieve4I15e -c 250 -f 100000000 -v -o snfs_100000k_I15.job.wsl.mpqs.out.1 -r snfs_100000k_I15.txt
bin/mpqs/gnfs-lasieve4I16e -c 200 -f 250000000 -v -o R1340L_poly.job.wsl.mpqs.out.1 -a R1340L_poly.txt
bin/mpqs/gnfs-lasieve4I16e -c 200 -f 325000000 -v -o R942_poly.job.wsl.mpqs.out.1 -r R942_poly.txt

bin/gnfs-lasieve4I11e -c 10000 -f 650000 -v -o nfs300_650k_I11.job.wsl.out.1 -a nfs300_650k_I11.job 
bin/gnfs-lasieve4I12e -c 5000 -f 900000 -v -o nfs330_900k_I12.job.wsl.out.1 -a nfs330_900k_I12.job 
bin/gnfs-lasieve4I13e -c 1000 -f 1600000 -v -o nfs360_1600k_I13.job.wsl.out.1 -a nfs360_1600k_I13.job 
bin/gnfs-lasieve4I13e -c 1000 -f 2700000 -v -o nfs390_2700k_I13.job.wsl.out.1 -a nfs390_2700k_I13.job 
bin/gnfs-lasieve4I13e -c 1000 -f 6500000 -v -o nfs420_6500k_I13.job.wsl.out.1 -a nfs420_6500k_I13.job 
bin/gnfs-lasieve4I14e -c 1000 -f 15300000 -v -o nfs512_15300k_I14.job.wsl.out.1 -a nfs512_15300k_I14.job 
bin/gnfs-lasieve4I16e -c 250 -f 50000000 -v -o nfs631_50000k_I16.job.wsl.out.1 -a nfs631_50000k_I16.job 
bin/gnfs-lasieve4I15e -c 250 -f 100000000 -v -o snfs_100000k_I15.job.wsl.out.1 -r snfs_100000k_I15.txt
bin/gnfs-lasieve4I16e -c 200 -f 250000000 -v -o R1340L_poly.job.wsl.out.1 -a R1340L_poly.txt
bin/gnfs-lasieve4I16e -c 200 -f 325000000 -v -o R942_poly.job.wsl.out.1 -r R942_poly.txt

bin/avx512/gnfs-lasieve4I11e -c 10000 -f 650000 -v -o nfs300_650k_I11.job.wsl.avx512.out.1 -a nfs300_650k_I11.job 
bin/avx512/gnfs-lasieve4I12e -c 5000 -f 900000 -v -o nfs330_900k_I12.job.wsl.avx512.out.1 -a nfs330_900k_I12.job 
bin/avx512/gnfs-lasieve4I13e -c 1000 -f 1600000 -v -o nfs360_1600k_I13.job.wsl.avx512.out.1 -a nfs360_1600k_I13.job 
bin/avx512/gnfs-lasieve4I13e -c 1000 -f 2700000 -v -o nfs390_2700k_I13.job.wsl.avx512.out.1 -a nfs390_2700k_I13.job 
bin/avx512/gnfs-lasieve4I13e -c 1000 -f 6500000 -v -o nfs420_6500k_I13.job.wsl.avx512.out.1 -a nfs420_6500k_I13.job 
bin/avx512/gnfs-lasieve4I14e -c 1000 -f 15300000 -v -o nfs512_15300k_I14.job.wsl.avx512.out.1 -a nfs512_15300k_I14.job 
bin/avx512/gnfs-lasieve4I16e -c 250 -f 50000000 -v -o nfs631_50000k_I16.job.wsl.avx512.out.1 -a nfs631_50000k_I16.job 
bin/avx512/gnfs-lasieve4I15e -c 250 -f 100000000 -v -o snfs_100000k_I15.job.wsl.avx512.out.1 -r snfs_100000k_I15.txt
bin/avx512/gnfs-lasieve4I16e -c 200 -f 250000000 -v -o R1340L_poly.job.wsl.avx512.out.1 -a R1340L_poly.txt
bin/avx512/gnfs-lasieve4I16e -c 200 -f 325000000 -v -o R942_poly.job.wsl.avx512.out.1 -r R942_poly.txt
