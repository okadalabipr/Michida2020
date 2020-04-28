#!/bin/bash

#Specify raw fastq file you want to process with this script at the first argument.
#Please change the file names you like.

index="path to your mm10 index made by bowtie"
adapter="path to your Trimmomatic/adapters/TruSeq2-SE.fa"
save_dir="path to your directory where you want to save all files"

#Adapter trimming with trimmomatic
trimmomatic SE -threads 4 -phred33 $1 ${save_dir}/trimmed.fastq ILLUMINACLIP:${adapter}:2:30:10 LEADING:30 TRAILING:30 SLIDINGWINDOW:4:15 MINLEN:30 HEADCROP:2

#Mapping with bowtie
bowtie -S -p 4 -m 1 ${index} -q ${save_dir}/trimmed.fastq ${save_dir}/mapped.sam

#Convert sam to bam with samtools
samtools view -@ 4 -b -o ${save_dir}/mapped.bam ${save_dir}/mapped.sam

#Sort bam with samtools
samtools sort -@ 4 ${save_dir}/mapped.bam -o ${save_dir}/mapped_sorted.bam

#Remove unmapped reads with samtools
samtools view -@ 4 -b -F 4 ${save_dir}/mapped_sorted.bam > ${save_dir}/mapped_sorted_rm_unmapped.bam

#Remove duplicates with picard
picard MarkDuplicates I=${save_dir}/mapped_sorted_rm_unmapped.bam O=${save_dir}/all_done.bam M=${save_dir}/report.txt REMOVE_DUPLICATES=true

#Remove intermediate files
rm -rf ${save_dir}/mapped.sam ${save_dir}/mapped.bam ${save_dir}/mapped_sorted.bam ${save_dir}/mapped_sorted_rm_unmapped.bam