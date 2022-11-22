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
    #/dellfsqd2/ST_OCEAN/USER/huangzhi/wt/gem2cfm/get_roi.sh $spname $slname >roi.json
    /dellfsqd2/ST_OCEAN/USER/guolidong/wochong_rna_202201/scripts_confs_backup/get_roi.sh  $spname $slname >roi.json
    python3 /dellfsqd2/ST_OCEAN/USER/guolidong/wochong_rna_202201/scripts_confs_backup/gem_to_cfm300.py  -s $ssdna  -g $gem -b $border -m $mask -r roi.json -o $prefix -x $minx -y $miny >log 2>err
    cd -
}





run_one b 1   0	6801   		&
#run_one b 2  10001	10   	&
run_one b 3   15	10   	&
run_one b 4   7001	0   	&
run_one b 5   1775	9001   	&
run_one b 6   7256	0   	&
run_one b 7   8501	0   	&
run_one b 8   5	0     		&
run_one b 9   9501	0   	&
run_one b 10  0	0     		&
run_one b 11  0	0     		&
run_one b 12  1775	2025   	&
run_one b 13  0	0     		&
#run_one b 14 0	0     		&
run_one b 15  0	7501   		&
run_one b 16  0	1     		&
run_one b 17  2025	10   	&
run_one b 18  15	7501   	&
run_one b 19  0	0     		&
run_one b 20  2426	1975   	&
run_one b 21  3375	1976   	&
wait 
#run_one a 1          &
run_one a 2    12	2175      &
run_one a 3    15	10       &
#run_one a 4   2	0         &
run_one a 5    2	0         &
run_one a 6    9001	1775      &
run_one a 7    10	0        &
run_one a 8    8501	2225      &
run_one a 9    1825	9001      &
run_one a 10   5879	9251      &
run_one a 11   8351	2175      &
run_one a 12   1	7001      &
run_one a 13   0	4         &
run_one a 14   2425	8001      &
run_one a 15   1	7001      &
run_one a 16   0	1         &
run_one a 17   0	6501      &
run_one a 18   7501	1375      &
run_one a 19   7501	1      &
run_one a 20   0	0         & 
run_one a 21   3375	1976      &
wait 
run_one a 22   7256	0      &
run_one a 23   2225	2375      &
run_one a 24   8001	0      &
run_one a 25   2875	1575      &
run_one a 26   1025	1775      &
run_one a 27   925	2226      &
run_one a 28   1837	2425      &
run_one a 29   5875	1775      &
run_one a 30   10501	1975      &
run_one a 31   2225	10001      &
run_one a 32   0	0         &
run_one a 33   2875	10001      &
run_one a 34   0	0         &
run_one a 35   9001	0      &
wait 
