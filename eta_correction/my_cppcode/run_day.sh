#!/bin/bash

y=20$1
m=$(printf "%02d" $2)
d=$(printf "%02d" $3)

export FILELIST=`pwd`/DataList/$y/$m/datalist_$y$m$d.txt
export Nlines=5
Ntot=($(wc -l $FILELIST))
Ntot=${Ntot[0]}
Nmax=$(($Ntot / $Nlines))
ls $FILELIST
sbatch --array=0-$Nmax%200 submit_array.sh
