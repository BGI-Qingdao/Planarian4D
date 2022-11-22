#!/bin/bash

#for x in WT 0hpa1 0hpa2 12hpa1 12hpa2 36hpa1 36hpa2 3dpa1 3dpa2 5dpa1 5dpa2 7dpa1 7dpa2 10dpa1 10dpa2 14dpa1 14dpa2
#do
#    ./FISH_3c.py -i cell_pos_db/$x".txt" -m "$x".7.density_list.txt -y "$x".2.density_list.txt -g "$x".5.density_list.txt -o "$x"_t3  --xmin '-300' --ymin '-900' &
#done

for x in WT 0hpa1 0hpa2 12hpa1 12hpa2 36hpa1 36hpa2 3dpa1 3dpa2 5dpa1 5dpa2 7dpa1 7dpa2 10dpa1 10dpa2 14dpa1 14dpa2
do
    ./FISH_3c.py -i cell_pos_db/$x".txt" -m "$x".7.density_list.txt  -o "$x"_c7  --xmin '-300' --ymin '-900' &
done
wait
for x in WT 0hpa1 0hpa2 12hpa1 12hpa2 36hpa1 36hpa2 3dpa1 3dpa2 5dpa1 5dpa2 7dpa1 7dpa2 10dpa1 10dpa2 14dpa1 14dpa2
do
    ./FISH_3c.py -i cell_pos_db/$x".txt" -y "$x".2.density_list.txt  -o "$x"_c2  --xmin '-300' --ymin '-900' &
done
wait
for x in WT 0hpa1 0hpa2 12hpa1 12hpa2 36hpa1 36hpa2 3dpa1 3dpa2 5dpa1 5dpa2 7dpa1 7dpa2 10dpa1 10dpa2 14dpa1 14dpa2
do
    ./FISH_3c.py -i cell_pos_db/$x".txt" -g "$x".5.density_list.txt  -o "$x"_c5  --xmin '-300' --ymin '-900' &
done
wait
