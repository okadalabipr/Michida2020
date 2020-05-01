#!/bin/bash

SE_BED="path to your BED file of SEs made through running Make_enhancer_catalog.sh"
TE_BED="path to your BED file of TEs made through running Make_enhancer_catalog.sh"
ATAC_bigwig="path to your bigiwig file made by Make_bigwig.sh"

compute_matrix_object="path to output object of computeMatrix (deepTools) (.gz)"

computeMatrix scale-regions -p 4 -m 10000 \
-S ${ATAC_bigwig} \
-R ${SE_BED} \
   ${TE_BED} \
-a 5000 \
-b 5000 \
-bs 1 \
--averageTypeBins sum \
--skipZeros \
-out ${compute_matrix_object}

PDF="path to output figure pdf"
TXT="path to output text, in this paper, SFig.1D was plotted using this text file by SFig1D.sh"

plotProfile -m ${compute_matrix_object} \
--plotType=se --averageType mean \
--plotHeight 9 --plotWidth 9 \
--numPlotsPerRow 2 \
--colors "#e34a33" "gray" "#31a354" \
-z "SE" "TE" \
-y "Average RelA Signal (read)" \
-out ${PDF} --outFileNameData ${TXT}
