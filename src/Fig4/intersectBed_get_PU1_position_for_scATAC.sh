#!/bin/bash

SE="path to your SE catalog made by Make_enhancer_catalog.sh"
TE="path to your TE catalog made by Make_enhancer_catalog.sh"
motif="path to your PU.1 motif file. In this paper, PU1_homer.motif provided in /data of this repository was used."

cd dir_scATAC/peakcall

dirs=`find * -type d`

for d in ${dirs};do
    cd ${d}
    fromdos *.bed #For windows users.
    intersectBed -a ${SE} -b ${d}_fixed_peaks.bed > ${d}_ATAC_in_SE.bed
    intersectBed -a ${TE} -b ${d}_fixed_peaks.bed > ${d}_ATAC_in_TE.bed
    bed2pos.pl ${d}_ATAC_in_SE.bed > ${d}_ATAC_in_SE.hb
    bed2pos.pl ${d}_ATAC_in_TE.bed > ${d}_ATAC_in_TE.hb

    findMotifsGenome.pl ${d}_ATAC_in_SE.hb mm10 ./ -find ${motif} > ${d}_PU1_in_ATAC_in_SE.txt
    mv motifFindingParameters.txt SE_motifFindingParameters.txt
    findMotifsGenome.pl ${d}_ATAC_in_TE.hb mm10 ./ -find ${motif} > ${d}_PU1_in_ATAC_in_TE.txt
    mv motifFindingParameters.txt TE_motifFindingParameters.txt
    cd ..
done
