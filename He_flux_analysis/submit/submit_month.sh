#!/bin/bash

#for i in /beegfs/dampe/prod/FM/FlightData/2A/20181019/DAMPE_2A_OBS_20181019_*/*.root; do

for splitfile in ../resources/201810_split/*; do
	sbatch run.sh $splitfile
done
