#!/bin/bash

export FILELIST=`pwd`/datalist/skim/all_002_010.txt
export Nlines=5
Ntot=($(wc -l $FILELIST))
Ntot=${Ntot[0]}
Nmax=$(($Ntot / $Nlines))
ls $FILELIST
sbatch --array=0-$Nmax%200 submit_array.sh
