#!/bin/bash

filtered_Gained_SE_csv="path to you filtered_Gained_SE.csv made by Filtering_genes.R"
filtered_Gained_SE_bed="BED file of filtered gained SEs will be made in this script. Please specify its path here."
filtered_Gained_TE_csv="path to you filtered_Gained_TE.csv made by Filtering_genes.R"
filtered_Gained_TE_bed="BED file of filtered gained TEs will be made in this script. Please specify its path here."
extended_bulk_ATAC_peak_IgM="path to your extended bulk ATAC-seq peaks in anti-IgM 060 min made by Extend_bulk_ATAC_peaks.R"

#Make BED file of filtered gained SEs
sed -E "s/,/\t/g" ${filtered_Gained_SE_csv} | sed '1d' | sed -E "s/\"//g" | awk -v 'OFS=\t' '{print $2, $3, $4}' > ${filtered_Gained_SE_bed} #This command may be required to change depending on your data handling.

#Make BED file of filtered gained TEs
sed -E "s/,/\t/g" ${filtered_Gained_TE_csv} | sed '1d' | sed -E "s/\"//g" | awk -v 'OFS=\t' '{print $2, $3, $4}' > ${filtered_Gained_TE_bed} #This command may be required to change depending on your data handling.

#Get ATAC peaks in gained SEs and TEs in anti-IgM 060 min.
ATAC_in_Gained_SE="path to output file of below command"
intersectBed -a ${filtered_Gained_SE_bed} -b ${extended_bulk_ATAC_peak_IgM} > ${ATAC_in_Gained_SE}

ATAC_in_Gained_TE="path to output file of below command"
intersectBed -a ${filtered_Gained_TE_bed} -b ${extended_bulk_ATAC_peak_IgM} > ${ATAC_in_Gained_TE}

#Convert BED file into HOMER BED format.
ATAC_in_Gained_SE_hb="path to output file of below command"
bed2pos.pl ${ATAC_in_Gained_SE} > ${ATAC_in_Gained_SE_hb}

ATAC_in_Gained_TE_hb="path to output file of below command"
bed2pos.pl ${ATAC_in_Gained_TE} > ${ATAC_in_Gained_TE_hb}

#Search PU.1 motifs in ATAC-seq peak regions in gained SEs and TEs in anti-IgM 060 min condition.
PU1_in_ATAC_in_SE_IgM="path to output file of below command"
findMotifsGenome.pl ${ATAC_in_Gained_SE_hb} mm10 ./ -find PU1_homer.motif > ${PU1_in_ATAC_in_SE_IgM}

PU1_in_ATAC_in_TE_IgM="path to output file of below command"
findMotifsGenome.pl ${ATAC_in_Gained_TE_hb} mm10 ./ -find PU1_homer.motif > ${PU1_in_ATAC_in_TE_IgM}

#In this paper, PU.1 motifs in HOMER was used. PU1_homer.motif is provided in our repository.