#!/bin/bash

#To get Fig1.C and Sup.Fig1. C, signal intensity was obtained using multiBigwigSummary of deepTools.

#In anti-IgM 000 min condition

H3K27Ac_ctrl="path to your bigwig file of H3K27Ac ChIP-seq data of anti-IgM 000 min"
RelA_ctrl="path to your bigwig file of RelA ChIP-seq data of anti-IgM 000 min"
ATAC_ctrl="path to your bigwig file of ATAC-seq data of anti-IgM 000 min"

SE_ctrl="path to your superEnhancers.txt of anti-IgM 000 min data obtained by findPeaks of HOMER"
TE_ctrl="path to your typicalEnhancers.txt of anti-IgM 000 min data obtained by findPeaks of HOMER"

SE_ctrl_BED="HOMER peak file, superEnhancers.txt is converting into BED format using following command. Please write path to the converted file here."
TE_ctrl_BED="As in SE."

TSS="path to your BED file of TSS. In this paper, mm10.tss in HOMER was used. (homer/data/genomes/mm10/mm10.tss Please install mm10 into HOMER fist.)"

SE_ctrl_npz="path to npz file of SEs in anti-IgM 000 min condition made by multiBigwigSummary"
TE_ctrl_npz="path to npz file of TEs in anti-IgM 000 min condition made by multiBigwigSummary"
TSS_ctrl_npz="path to npz file of TSSs in anti-IgM 000 min condition made by multiBigwigSummary"

SE_ctrl_tab="path to tab file of SEs in anti-IgM 000 min condition made by multiBigwigSummary"
TE_ctrl_tab="path to tab file of TEs in anti-IgM 000 min condition made by multiBigwigSummary"
TSS_ctrl_tab="path to tab file of TSSs in anti-IgM 000 min condition made by multiBigwigSummary"

#Convert HOMER peak file into BED format
sed '/^#/d' ${SE_ctrl} | awk -v 'OFS=\t' '{print $2, $3, $4}' > ${SE_ctrl_BED}
sed '/^#/d' ${TE_ctrl} | awk -v 'OFS=\t' '{print $2, $3, $4}' > ${TE_ctrl_BED}

#Intensities in SEs in anti-IgM 000 min conditions
multiBigwigSummary BED-file -p max -bs 1 -b ${H3K27Ac_ctrl} ${RelA_ctrl} ${ATAC_ctrl} -out ${SE_ctrl_npz} --BED ${SE_ctrl_BED} --outRawCounts ${SE_ctrl_tab}

#Intensities in TEs in anti-IgM 000 min conditions
multiBigwigSummary BED-file -p max -bs 1 -b ${H3K27Ac_ctrl} ${RelA_ctrl} ${ATAC_ctrl} -out ${TE_ctrl_npz} --BED ${TE_ctrl_BED} --outRawCounts ${TE_ctrl_tab}

#Intensities in TSSs in anti-IgM 000 min conditions
multiBigwigSummary BED-file -p max -bs 1 -b ${H3K27Ac_ctrl} ${RelA_ctrl} ${ATAC_ctrl} -out ${TSS_ctrl_npz} --BED ${TSS} --outRawCounts ${TSS_ctrl_tab}

#In anti-IgM 060 min condition

H3K27Ac_IgM="path to your bigwig file of H3K27Ac ChIP-seq data of anti-IgM 000 min"
RelA_IgM="path to your bigwig file of RelA ChIP-seq data of anti-IgM 000 min"
ATAC_IgM="path to your bigwig file of ATAC-seq data of anti-IgM 000 min"

SE_IgM="path to your superEnhancers.txt of anti-IgM 000 min data obtained by findPeaks of HOMER"
TE_IgM="path to your typicalEnhancers.txt of anti-IgM 000 min data obtained by findPeaks of HOMER"

SE_IgM_BED="HOMER peak file, superEnhancers.txt is converting into BED format using following command. Please write path to the converted file here."
TE_IgM_BED="As in SE."

SE_IgM_npz="path to npz file of SEs in anti-IgM 000 min condition made by multiBigwigSummary"
TE_IgM_npz="path to npz file of TEs in anti-IgM 000 min condition made by multiBigwigSummary"
TSS_IgM_npz="path to npz file of TSSs in anti-IgM 000 min condition made by multiBigwigSummary"

SE_IgM_tab="path to tab file of SEs in anti-IgM 000 min condition made by multiBigwigSummary"
TE_IgM_tab="path to tab file of TEs in anti-IgM 000 min condition made by multiBigwigSummary"
TSS_IgM_tab="path to tab file of TSSs in anti-IgM 000 min condition made by multiBigwigSummary"

#Convert HOMER peak file into BED format
sed '/^#/d' ${SE_IgM} | awk -v 'OFS=\t' '{print $2, $3, $4}' > ${SE_IgM_BED}
sed '/^#/d' ${TE_IgM} | awk -v 'OFS=\t' '{print $2, $3, $4}' > ${TE_IgM_BED}

#Intensities in SEs in anti-IgM 000 min conditions
multiBigwigSummary BED-file -p max -bs 1 -b ${H3K27Ac_IgM} ${RelA_IgM} ${ATAC_IgM} -out ${SE_IgM_npz} --BED ${SE_IgM_BED} --outRawCounts ${SE_IgM_tab}

#Intensities in TEs in anti-IgM 000 min conditions
multiBigwigSummary BED-file -p max -bs 1 -b ${H3K27Ac_IgM} ${RelA_IgM} ${ATAC_IgM} -out ${TE_IgM_npz} --BED ${TE_IgM_BED} --outRawCounts ${TE_IgM_tab}

#Intensities in TSSs in anti-IgM 000 min conditions
multiBigwigSummary BED-file -p max -bs 1 -b ${H3K27Ac_IgM} ${RelA_IgM} ${ATAC_IgM} -out ${TSS_IgM_npz} --BED ${TSS} --outRawCounts ${TSS_IgM_tab}