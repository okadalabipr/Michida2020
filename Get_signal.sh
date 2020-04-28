#!/bin/bash

#To get signal intensity used in Fig.1D, annotatePeaks.pl of HOMER was used.

SE_catalog="path to your SE catalog (BED format)"
TE_catalog="path to your TE catalog (BED format)"

H3K27Ac_TagDir_ctrl="path to your Tag Directory of H3K27Ac ChIP-seq data in anti-IgM 000 min condition"
H3K27Ac_TagDir_IgM="path to your Tag Directory of H3K27Ac ChIP-seq data in anti-IgM 060 min condition"

RelA_TagDir_ctrl="path to your Tag Directory of RelA ChIP-seq data in anti-IgM 000 min condition"
RelA_TagDir_IgM="path to your Tag Directory of RelA ChIP-seq data in anti-IgM 060 min condition"

ATAC_TagDir_ctrl="path to your Tag Directory of ATAC-seq data in anti-IgM 000 min condition"
ATAC_TagDir_IgM="path to your Tag Directory of ATAC-seq data in anti-IgM 060 min condition"

SE_signal="path to the output file for SE catalog"
TE_signal="path to the output file for TE catalog"
TSS_signal="path to the output file for TSS"

#Signal intensity in SE catalog
annotatePeaks.pl ${SE_catalog} mm10 -norm 1000000 -size given -d ${H3K27Ac_TagDir_ctrl} ${H3K27Ac_TagDir_IgM} ${RelA_TagDir_ctrl} ${RelA_TagDir_IgM} ${ATAC_TagDir_ctrl} ${ATAC_TagDir_IgM} > ${SE_signal}

#Signal intensity in TE catalog
annotatePeaks.pl ${TE_catalog} mm10 -norm 1000000 -size given -d ${H3K27Ac_TagDir_ctrl} ${H3K27Ac_TagDir_IgM} ${RelA_TagDir_ctrl} ${RelA_TagDir_IgM} ${ATAC_TagDir_ctrl} ${ATAC_TagDir_IgM} > ${TE_signal}

#Signal intensity in TSS
annotatePeaks.pl tss mm10 -norm 1000000 -size given -d ${H3K27Ac_TagDir_ctrl} ${H3K27Ac_TagDir_IgM} ${RelA_TagDir_ctrl} ${RelA_TagDir_IgM} ${ATAC_TagDir_ctrl} ${ATAC_TagDir_IgM} > ${TSS_signal}