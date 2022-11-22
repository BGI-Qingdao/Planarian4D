

function get_data(){
    for x in 0hpa1 0hpa2 10dpa1 10dpa2 12hpa1 12hpa2 14dpa1 14dpa2 36hpa1 36hpa2 3dpa1 3dpa2 5dpa1 5dpa2 7dpa1 7dpa2 WT
    do
        sh run_pcgs_oneindv.sh $x  $1 &
    done
    wait
}

get_data 100
#get_data 75 
#get_data 50 
#get_data 25
get_data 10
