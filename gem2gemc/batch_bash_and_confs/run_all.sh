#!/bin/bash

function run_one(){
    sp=$1
    SP=${sp^^}
    sl=$2
    minx=$3
    miny=$4
    prefix="sample-"$sp"-slice-"$sl
    echo "Run $prefix now ..."
    mkdir -p $prefix 
    cd $prefix
    pwd 
    ssdna="slice_"$sp"."$sl".tif"
    #gem="slice_"$sp"."$sl".gem"
    gem="sample_"${sp}"_slice"${sl}".gem" # this is the new gem name of #2 round 

    mask="sample"$SP"_"$sl".mask.txt"
    border="sample"$SP"_"$sl".outlines.txt"

    ln -s  "/dellfsqd2/ST_OCEAN/USER/guolidong/trackline_dev/wochong_reagain/01.ssDNA/slice_"$sp"."$sl".tif"
    #ln -s  "/dellfsqd2/ST_OCEAN/USER/guolidong/trackline_dev/wochong_reagain/00.gems/slice_"$sp"."$sl".gem"
    ln -s  "/dellfsqd2/ST_OCEAN/USER/guolidong/wochong_rna_202201/s02.splited.gem/"$gem #this is the new gem path of #2 round
    ln -s  "/dellfsqd2/ST_OCEAN/USER/guolidong/trackline_dev/gem2cfm/sample"$SP"/sample"$SP"_"$sl".mask.txt"
    ln -s  "/dellfsqd2/ST_OCEAN/USER/guolidong/trackline_dev/gem2cfm/sample"$SP"/sample"$SP"_"$sl".outlines.txt"
    
    spname="sample_"$sp
    slname="slice_"$sl
    /dellfsqd2/ST_OCEAN/USER/guolidong/wochong_rna_202201/scripts_confs_backup/get_roi.sh $spname $slname >roi.json
    python3 /dellfsqd2/ST_OCEAN/USER/guolidong/wochong_rna_202201/scripts_confs_backup/gem_to_cfm.py  -s $ssdna  -g $gem -b $border -m $mask -r roi.json -o $prefix -x $minx -y $miny >log 2>err
    cd -
}


#run_one b 11 &

#run_one b 1 &
#run_one b 12 &
#run_one b 13 &
run_one b 14 0 0&
run_one b 2  10001	10 &
#run_one b 3 &
#run_one b 4 &
#run_one b 5 &
#run_one b 6 &
#run_one b 7 &
#run_one b 8 &
#run_one b 9 &
#run_one b 10 &

#run_one b 15 &
#run_one b 16 &
#run_one b 17 &
#run_one b 18 &
#run_one b 19 &
#run_one b 20 &
#run_one b 21 &

wait

#run_one a 1 &
#run_one a 2 &
#run_one a 3 &
run_one a 4  2	0 &
#run_one a 5 &
#run_one a 6 &
#run_one a 7 &
#run_one a 8 &
#run_one a 9 &
#run_one a 10 &


#run_one a 11 &
#run_one a 12 &
#run_one a 13 &
#run_one a 14 &
#run_one a 15 &
#run_one a 16 &
#run_one a 17 &
#run_one a 18 &
#run_one a 19 &
#run_one a 20 & 

#run_one a 21 &
#run_one a 22 &
#run_one a 23 &
#run_one a 24 &
#run_one a 25 &
#run_one a 26 &
#run_one a 27 &
#run_one a 28 &
#run_one a 29 &
#run_one a 30 &

wait


#run_one a 31 &
#run_one a 32 &
#run_one a 33 &
#run_one a 34 &
#run_one a 35 &


wait
