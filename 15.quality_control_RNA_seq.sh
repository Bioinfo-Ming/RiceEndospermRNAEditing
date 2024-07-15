#!/bin/bash

java -jar /xtdisk/apod-gen/chenming/00.software/Trimmomatic-0.35/trimmomatic-0.35.jar PE -threads 30 -phred33 \
/xtdisk/apod-gen/chenming/41.RNA_Editing_Leaf_Seedling/03.CRA011887/CRX743583/CRR826540_f1.fastq \
/xtdisk/apod-gen/chenming/41.RNA_Editing_Leaf_Seedling/03.CRA011887/CRX743583/CRR826540_r2.fastq \
/xtdisk/apod-gen/chenming/41.RNA_Editing_Leaf_Seedling/03.CRA011887/CRX743583/CRR826540_f1_paired.fastq \
/xtdisk/apod-gen/chenming/41.RNA_Editing_Leaf_Seedling/03.CRA011887/CRX743583/CRR826540_f1_unpaired.fastq \
/xtdisk/apod-gen/chenming/41.RNA_Editing_Leaf_Seedling/03.CRA011887/CRX743583/CRR826540_r2_paired.fastq \
/xtdisk/apod-gen/chenming/41.RNA_Editing_Leaf_Seedling/03.CRA011887/CRX743583/CRR826540_r2_unpaired.fastq \
ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
