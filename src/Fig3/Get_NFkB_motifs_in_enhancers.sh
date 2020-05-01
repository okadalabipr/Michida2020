#!/bin/bash

#This script is for SFig.4A.

ATAC_peak="path to your extended bulk ATAC-seq peaks of anti-IgM 060 min condition obtained using Peak_call_for_bulk.sh and Extend_bulk_ATAC_peaks.R"
ATAC_peak_HOMER_BED="BED format should be converted into HOMER BED (.hb) format when you input it in HOMER's program. Specify path of HOMER BED of ${ATAC_peak} here."
NFkB_motifs_in_ATAC_IgM="path to the output file"
motif="path to NFkB_homer.motif provided in our repository"

bed2pos.pl ${ATAC_peak} > ${ATAC_peak_HOMER_BED} #bed2pos.pl is HOMER's program which converts BED files into HOMER BED files. 

findMotifsGenome.pl ${ATAC_peak_HOMER_BED} mm10 ./ -find ${motif} > ${NFkB_motifs_in_ATAC_IgM}