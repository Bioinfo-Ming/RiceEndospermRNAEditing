#!/bin/bash

python2 REDItoolDnaRna.py \
-i CRR826538_Alignment-sorted.bam \
-j sorted.bam \
-f rice_14.fa \
-o ./ \
 -t 20 -c 30,30 -m 25,25 -q 25,25 -e -E -d -D -u -U -l -L -V \
 -G rice_all_sorted.gff3
