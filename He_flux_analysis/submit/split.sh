#!/bin/bash

mkdir -p ../resources/201810_split/

MONTH_FILE="../resources/201810.txt"
LINES=$(( `wc -l < $MONTH_FILE`/50 ))
LOCATION="../resources/201810_split/"

split -l $LINES $MONTH_FILE $LOCATION

PARTS=`ls $LOCATION | wc -l`

echo "$MONTH_FILE splitted into $PARTS parts and stored in $LOCATION"
