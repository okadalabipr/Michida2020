#!/bin/bash

annotation="path to your annotation file (In this paper, GENCODE, release M17 (GRCm38.p6), Comprehensive gene annotation ALL)"
count="path to output count file"

featureCounts -p -B -t exon -g gene_name -a ${annotation} -o ${count} \
"path to bam of RamDA_00000_A01" \
"path to bam of RamDA_00000_A02" \
...
"path to bam of RamDA_00000_H12" \
"path to bam of RamDA_00010_A01" \
...
"path to bam of RamDA_00010_H12" \
"path to bam of RamDA_00100_A01" \
...
"path to bam of RamDA_10000_H12"