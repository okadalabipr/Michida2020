#!/bin/bash

#Specify read 1 and read 2 fastq files you want to process with this script at the first argument and the second argument respectively.
#Please change the file names you like.
#Both bulk and single-cell ATAC-seq fastq files were processed by this script in the paper.

index="path to your mm10 index made by bowtie2"
adapter="path to your Trimmomatic/adapters/NexteraPE-PE.fa"
save_dir="path to your directory where you want to save all files"

#Adapter trimming with trimmomatic
trimmomatic PE -threads 4 -phred33 $1 $2 \
${save_dir}/trimmed_R1.fastq ${save_dir}/unpaired_R1.fastq \
${save_dir}/trimmed_R2.fastq ${save_dir}/unpaired_R2.fastq \
ILLUMINACLIP:${adapter}:2:30:10 TRAILING:30 LEADING:30 MINLEN:50

#Mapping with bowtie2
bowtie2 --maxins 2000 --minins 0 -p 4 --no-discordant --no-mixed --dovetail --very-sensitive -t -x ${index} \
-1 ${save_dir}/trimmed_R1.fastq -2 ${save_dir}/trimmed_R2.fastq -S ${save_dir}/mapped.sam

#Convert sam to bam with samtools
samtools view -@ 4 -bS ${save_dir}/mapped.sam -o ${save_dir}/mapped.bam

#Sort bam with samtools
samtools sort -@ 4 ${save_dir}/mapped.bam -o ${save_dir}/mapped_sorted.bam

#Remove duplicates with picard
picard MarkDuplicates I=${save_dir}/mapped_sorted.bam O=${save_dir}/mapped_sorted_rm_dups.bam M=${save_dir}/report.txt REMOVE_DUPLICATES=true

#Extract signle mapped reads with samtools
samtools view -@ 4 -b -q 30 ${save_dir}/mapped_sorted_rm_dups.bam > ${save_dir}/all_done.bam

#Remove intermediate files
rm -rf ${save_dir}/unpaired_R1.fastq ${save_dir}/unpaired_R2.fastq ${save_dir}/mapped.sam \
${save_dir}/mapped.bam ${save_dir}/mapped_sorted.bam ${save_dir}/mapped_sorted_rm_dups.bam