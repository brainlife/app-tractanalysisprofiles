#!/bin/bash
#PBS -l nodes=1:ppn=4,vmem=28gb,walltime=3:00:00
#PBS -N tractanalysisprofiles
#PBS -V

set -e

singularity exec -e docker://brainlife/mcr:neurodebian1604-r2017a ./compiled/main

count=$(ls images/*.png | wc -l)
NUMFILES=`cat numfiles.txt`
if [ "$count" !=  "$NUMFILES" ];
then
	echo "not all output generated"
	exit 1
fi
