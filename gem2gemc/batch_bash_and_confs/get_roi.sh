#!/bin/bash

grep $1 /dellfsqd2/ST_OCEAN/USER/huangzhi/wt/small/check/result3.txt | \
                                                                            grep -w $2 | \
awk -F '\t' 'BEGIN{
                 print "[" 
             } 
             { 
                  if ( NR>1 ) printf(",\n");
                  printf("\t[\"%s\",[%s,%s,%s,%s], [%s,%s,%s,%s], %s]",$3,$4,$5,$6,$7, $8,$9,$10,$11, $12);
             } 
             END{ 
                 printf("\n]");
} ' 
