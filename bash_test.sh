#!/bin/bash
#times=(100 200 300 400 500 600 700 800 900 1000 1100 1200 1300 1400 1500 1600)
times=( 100 200 )

for i in "${times[@]}"; do
    echo Processing: $i
    i2=$((i + 200))
    python Tree_analysis.py $i $i2
done