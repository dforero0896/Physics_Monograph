#!/bin/bash

START=$(date +%s.%N)
./run.sh
END=$(date +%s.%N)
DIFF=$(echo "$END - $START" | bc)
echo Linear code: > time_comparison.txt
echo $DIFF >> time_comparison.txt

START=$(date +%s.%N)
./run_parallel.sh
END=$(date +%s.%N)
DIFF=$(echo "$END - $START" | bc)
echo Parallel code: >> time_comparison.txt
echo $DIFF >> time_comparison.txt
