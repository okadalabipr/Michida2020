#!/bin/bash

#Map single-cell ATAC-seq fastqs using ATAC-seq_mapping.sh and put all bam files in the same directory which contains only the bams. 
#In this repository, we call the directory you prepared above "dir_scATAC".

cd dir_scATAC

#Peak call using MACS2
files=`find *.bam`

mkdir peakcall #Make dir_scATAC/peakcall

for f in ${files};do
    macs2 callpeak -t ${f} -g mm --outdir ./peakcall/${f} -n ${f} --nomodel --nolambda --keep-dup all --call-summits
done

#Convert *peaks.xls into simple tsv files used in downstream processing.
dirs=`find * -type d`

for d in ${dirs};do
    cd ${d}
    cat ${d}_peaks.xls | grep -v '^#' | sed '/^$/d' > ${d}_peaks.tsv
    cd ..
done