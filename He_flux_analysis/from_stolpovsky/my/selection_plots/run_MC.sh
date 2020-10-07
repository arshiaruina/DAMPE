#!/bin/bash

# Use as ./run_MC.sh <filelist.txt> <Number of lines per job>
# Recomendations to choose the <N lines per job>
# Different MC have different number of events per file.
# Recommended N_lines = -5e-3 * N_ev + 100
# If N_ev is really low (<= 100) N_lines can be 200 or even more
# This script can also be used for the skim data analysis

export FILELIST=$1
export Nlines=$2
Ntot=($(wc -l $FILELIST))
Ntot=${Ntot[0]}
Nmax=$(($Ntot / $Nlines))
sbatch --array=0-$Nmax%200 submit_array.sh
