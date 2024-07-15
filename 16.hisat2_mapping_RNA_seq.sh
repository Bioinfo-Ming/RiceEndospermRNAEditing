#!/bin/bash

hisat2 -p 20 -x rice_tran -1 CRR826539_f1_paired.fastq -2 CRR826539_r2_paired.fastq -S CRR826539_Alignment-unsorted.sam 2> CRR826539_hisat2_Mapping_Rate.txt

samtools view -bS -T rice_14.fa CRR826539_Alignment-unsorted.sam > CRR826539_Alignment-unsorted.bam

samtools sort -@ 20 -o CRR826539_Alignment-sorted.bam CRR826539_Alignment-unsorted.bam

samtools index CRR826539_Alignment-sorted.bam

