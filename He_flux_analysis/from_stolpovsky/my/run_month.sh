#!/bin/bash

y=20$1
m=$(printf "%02d" $2)

export FILELIST=`pwd`/DataList/$y/$m/datalist_month_$y$m.txt
export Nlines=5
Ntot=($(wc -l $FILELIST))
Ntot=${Ntot[0]}
Nmax=$(($Ntot / $Nlines))
ls $FILELIST
sbatch --array=0-$Nmax%200 submit_array.sh
