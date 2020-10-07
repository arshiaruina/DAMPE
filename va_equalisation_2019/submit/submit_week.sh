#!/bin/bash

#for i in /beegfs/dampe/prod/FM/FlightData/2A/20181019/DAMPE_2A_OBS_20181019_*/*.root; do

#for ex in eta_executables/eta_tq*; do

#ex=eta_executables/eta_tq3.exe
#ex=eta.exe
#ex=../va_equalisation

ex=../va_eq_12
#splitfile=../resources/201810_wk1_split/aa
for splitfile in ../resources/201812_wk1_split/*; do
#for splitfile in ../resources/201811_wk1_split/*; do
	sbatch run.sh $ex $splitfile
done

#done
