#!/bin/bash

mkdir -p ../resources/201812_wk1_split/

WEEK_FILE="../resources/201812_wk1.txt"
LINES=$(( `wc -l < $WEEK_FILE`/90 ))
LOCATION="../resources/201812_wk1_split/"

split -l $LINES $WEEK_FILE $LOCATION

PARTS=`ls $LOCATION | wc -l`

echo "$WEEK_FILE splitted into $PARTS parts and stored in $LOCATION"
