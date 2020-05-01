#!/bin/bash

#All bigwig files were made by bamCoverage of deepTools.

bam="path to your bam file you want to convert into a bigwig file"
bigwig="path to the bigwig file you are making by this script"

bamCoverage -b ${bam} -p 4 --ignoreDuplicates --normalizeUsing RPGC --effectiveGenomeSize 2150570000 --binSize 1 --ignoreForNormalization chrX -o ${bigwig}