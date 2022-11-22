#!/bin/bash

python3 ../scripts/CellColocation.py  ../s00.input/NeibourCell.label.xyz.txt WT 30 50 >WT.dot &

function draw_one(){
    python3 ../scripts/CellColocation.py  ../s00.input/NeibourCell.label.xyz.txt $1 30 50 >${1}.dot
    dot -T png -o ${1}.png ${1}.dot
}


for x in 0hpa1 12hpa2 36hpa2 3dpa2 5dpa1 7dpa2 10dpa1 14dpa1 
do
    draw_one $x &
done

wait
dot -T png -o WT.png WT.dot


