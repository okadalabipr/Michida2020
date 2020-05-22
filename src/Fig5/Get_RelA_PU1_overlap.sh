#!/bin/bash

PU1_in_ATAC_in_SE_IgM="path to your PU.1 motifs positions in ATAC-seq peaks in gained SEs in anti-IgM 060 min condition obtained using Get_PU1_motifs_in_enhancers.sh and Get_motif_position.R"
PU1_in_ATAC_in_TE_IgM="path to your PU.1 motifs positions in ATAC-seq peaks in gained TEs in anti-IgM 060 min condition obtained using Get_PU1_motifs_in_enhancers.sh and Get_motif_position.R"
PU1_in_ATAC_in_SE_TE_IgM="path to output of below command"
cat ${PU1_in_ATAC_in_SE_IgM} ${PU1_in_ATAC_in_TE_IgM} > ${PU1_in_ATAC_in_SE_TE_IgM} #${PU1_in_ATAC_in_SE_TE_IgM} is required for for SFig6A.R.

RelA_in_ATAC_IgM="path to RelA ChIP-seq peaks in extended bulk ATAC-seq peaks in anti-IgM 060 min condition (BED format) obtained using Get_RelA_peaks_in_ATAC_peaks.sh"
PU1_in_RelA_in_ATAC_IgM="path to output file of below command"
intersectBed -a ${RelA_in_ATAC_IgM} -b ${PU1_in_ATAC_in_SE_TE_IgM} > ${PU1_in_RelA_in_ATAC_IgM}
