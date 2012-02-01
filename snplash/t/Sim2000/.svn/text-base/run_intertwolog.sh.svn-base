#!/bin/tcsh

#PBS -q x86_64
#PBS -l nodes=1:ppn=4
#PBS -l walltime=960:00:00
#PBS -l cput=400:00:00
#PBS -W group_list=langefeldGrp
#PBS -m a
#PBS -M guyrt7@wfu.edu
#PBS -j oe
#PBS -l mem=3GB
#PBS -l pmem=3GB

cd $PBS_O_WORKDIR

../../src/snplash -geno sim2000.geno -phen sim2000.covphen -map sim2000.map -engine intertwolog -out sim2000.intertwolog_cov -cov cov1 -beg 1 -end 500 

