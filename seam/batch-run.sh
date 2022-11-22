#!/bin/bash

for x  in 0hpa1 0hpa2 12hpa1 12hpa2 36hpa1 36hpa2 3dpa1 3dpa2 5dpa1 5dpa2 7dpa1 7dpa2 10dpa1 10dpa2 14dpa1 14dpa2 WT
do
    mkdir -p $x
    cd $x 
    ../scripts/gemc_3d.py $x
    cd -
done
