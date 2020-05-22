#!/bin/bash

index="path to your mm10 index made by STAR"
save_dir="path to your directory where you want to save all files"
R1="path to your read 1 fastq file you want to map"
R2="path to your read 2 fastq file you want to map"
trimmed1="path to your trimmed fastq of read 1"
trimmed2="path to your trimmed fastq of read 2"
prefix="Please specify name prefix following the rule below."
#Name prefix rule
#If you are processing "RamDA_00000_A01_1.fastq.gz" and "RamDA_00000_A01_2.fastq.gz", name prefix is "RamDA_00000_A01_".
#This is required if you perform SCDE using our code.

#Adapter trimming with Trim Galore
trim_galore --trim1 --paired ${R1} ${R2}

#Mapping with STAR
STAR --genomeDir ${index} \
--runThreadN 4 --readFilesIn ${trimmed1} ${trimmed2} \
--outFileNamePrefix ${prefix} \
--readFilesCommand zcat \
--outSAMtype BAM SortedByCoordinate \
--outSAMmultNmax 1
