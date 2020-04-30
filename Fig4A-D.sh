#!/bin/bash

Enhancer="path to your BED file of SE catalog or TE catalog made by Make_enhancer_catalog.sh"

#Control
ATAC_ctrl="path to your extended ATAC-seq peaks (BED format) in anti-IgM 000 min condition made by peak_call_for_bulk.sh and Extend_bulk_ATAC_peaks.R"
RelA_ctrl="path to your RelA ChIP-seq peaks (BED format) in anti-IgM 000 min condition made by peak_call_bulk.sh"

ATAC_in_Enhancer_ctrl="path to output of below command"
intersectBed -a ${Enhancer} -b ${ATAC_ctrl} > ${ATAC_in_Enhancer_ctrl} #Get ATAC-seq peaks in enhancer region

ATAC_in_Enhancer_ctrl_hb="path to output of below command"
bed2pos.pl ${ATAC_in_Enhancer_ctrl} > ${ATAC_in_Enhancer_ctrl_hb} #Convert BED into HOMER BED

RelA_in_ATAC_in_Enhancer_ctrl="path to output of below command"
intersectBed -a ${ATAC_in_Enhancer_ctrl} -b ${RelA_ctrl} > ${RelA_in_ATAC_in_Enhancer_ctrl} #Get RelA ChIP-seq peaks in ATAC-seq peaks in enhancer region

RelA_ATAC_in_Enhancer_ctrl_hb="path to output of below command"
bed2pos.pl ${RelA_ATAC_in_Enhancer_ctrl} > ${RelA_ATAC_in_Enhancer_ctrl_hb} #Convert BED into HOMER BED

#anti-IgM
ATAC_IgM="path to your extended ATAC-seq peaks (BED format) in anti-IgM 060 min condition made by peak_call_for_bulk.sh and Extend_bulk_ATAC_peaks.R"
RelA_IgM="path to your RelA ChIP-seq peaks (BED format) in anti-IgM 060 min condition made by peak_call_bulk.sh"

ATAC_in_Enhancer_IgM="path to output of below command"
intersectBed -a ${Enhancer} -b ${ATAC_IgM} > ${ATAC_in_Enhancer_IgM} #Get ATAC-seq peaks in enhancer region

ATAC_in_Enhancer_IgM_hb="path to output of below command"
bed2pos.pl ${ATAC_in_Enhancer_IgM} > ${ATAC_in_Enhancer_IgM_hb} #Convert BED into HOMER BED

RelA_in_ATAC_in_Enhancer_IgM="path to output of below command"
intersectBed -a ${ATAC_in_Enhancer_IgM} -b ${RelA_IgM} > ${RelA_in_ATAC_in_Enhancer_IgM} #Get RelA ChIP-seq peaks in ATAC-seq peaks in enhancer region

RelA_ATAC_in_Enhancer_IgM_hb="path to output of below command"
bed2pos.pl ${RelA_ATAC_in_Enhancer_IgM} > ${RelA_ATAC_in_Enhancer_IgM_hb} #Convert BED into HOMER BED

#findMotifsGenome.pl for Fig.4A-D
RelA_in_ATAC_in_Enhancer_ctrl_dir="path to output directory of below command"
findMotifsgenome.pl ${RelA_ATAC_in_Enhancer_ctrl_hb} mm10 ${RelA_in_ATAC_in_Enhancer_ctrl_dir} -size given -bg ${RelA_ATAC_in_Enhancer_IgM_hb} -p 4 #Fig.4A,C

RelA_in_ATAC_in_Enhancer_IgM_dir="path to output directory of below command"
findMotifsgenome.pl ${RelA_ATAC_in_Enhancer_IgM_hb} mm10 ${RelA_in_ATAC_in_Enhancer_IgM_dir} -size given -bg ${RelA_ATAC_in_Enhancer_ctrl_hb} -p 4 #Fig.4B,D

#findMotifsGenome.pl for SFig.5C-F
ATAC_in_Enhancer_ctrl_dir="path to output directory of below command"
findMotifsgenome.pl ${ATAC_in_Enhancer_ctrl_hb} mm10 ${ATAC_in_Enhancer_ctrl_dir} -size given -bg ${ATAC_in_Enhancer_IgM_hb} -p 4 #SFig.5C,E

ATAC_in_Enhancer_IgM_dir="path to output directory of below command"
findMotifsgenome.pl ${ATAC_in_Enhancer_IgM_hb} mm10 ${ATAC_in_Enhancer_IgM_dir} -size given -bg ${ATAC_in_Enhancer_ctrl_hb} -p 4 #SFig.5D,F

#Figures were made from knownResults.html