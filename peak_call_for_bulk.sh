#!/bin/bash

#RelA ChIP-seq of anti-IgM 060 min
RelA_IgM_bam="path to ypur bam file of RelA ChIP-seq of anti-IgM 060 min condition"
RelA_IgM_out="path to the output directory for RelA ChIP-seq of anti-IgM 060 min peak call"
macs2 callpeak -f BAM -t ${RelA_IgM_bam} -g mm --outdir ${RelA_IgM_out} -n RelA_IgM
cat ${RelA_IgM_out}/RelA_IgM_peaks.xls | grep -v '^#' | sed '/^$/d' | sed '1d' | awk -v 'OFS=\t' '{print $1, $2, $3}' > ${RelA_IgM_out}/RelA_IgM_peaks.bed #Convert peak.xls into BED format

#ATAC-seq of anti-IgM 000 min
ATAC_ctrl_bam="path to ypur bam file of ATAC-seq of anti-IgM 000 min condition"
ATAC_ctrl_out="path to the output directory for ATAC-seq of anti-IgM 000 min peak call"
macs2 callpeak -f BAM -t ${ATAC_ctrl_bam} -g mm --outdir ${ATAC_ctrl_out} -n ATAC_ctrl --keep-dup all
cat ${ATAC_ctrl_out}/ATAC_ctrl_peaks.xls | grep -v '^#' | sed '/^$/d' | sed '1d' | awk -v 'OFS=\t' '{print $1, $2, $3}' > ${ATAC_ctrl_out}/ATAC_ctrl.bed #Convert peak.xls into BED format

#ATAC-seq of anti-IgM 060 min
ATAC_IgM_bam="path to ypur bam file of ATAC-seq of anti-IgM 060 min condition"
ATAC_IgM_out="path to the output directory for ATAC-seq of anti-IgM 060 min peak call"
macs2 callpeak -f BAM -t ${ATAC_IgM_bam} -g mm --outdir ${ATAC_IgM_out} -n ATAC_IgM --keep-dup all
cat ${ATAC_IgM_out}/ATAC_IgM_peaks.xls | grep -v '^#' | sed '/^$/d' | sed '1d' | awk -v 'OFS=\t' '{print $1, $2, $3}' > ${ATAC_IgM_out}/ATAC_IgM.bed #Convert peak.xls into BED format