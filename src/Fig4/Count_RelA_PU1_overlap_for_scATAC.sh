#!/bin/bash

#Please redirect output of this script like bash Count_RelA_PU1_overlap_for_scATAC.sh > overlap_count_scATAC.txt

cd dir_scATAC/peakcall

RelA_ctrl="path to your RelA ChIP-seq peaks in anti-IgM 000 min condition (BED format) obtained by Peak_call_for_bulk.sh"

dirs=`find * -type d`

echo "SE"

for d in ${dirs};do
    cd ${d}
    fromdos *.bed #For windows users.
    intersectBed -a ${RelA_ctrl} -b ${d}_PU1_in_ATAC_in_SE.bed | wc -l
    cd ..
done

echo "TE"

for d in ${dirs};do
    cd ${d}
    fromdos *.bed #For windows users.
    intersectBed -a ${RelA_ctrl} -b ${d}_PU1_in_ATAC_in_TE.bed | wc -l
    cd ..
done