#!/bin/bash

#Ctrl
#Get SE and TE
TagDir_ctrl="path to your Tag Directory of H3K27Ac anti-IgM 000 min data"
findPeaks ${TagDir_ctrl} -style super -typical ${TagDir_ctrl}/typicalEnhancers.txt -L 0 -fdr 0.00005 -o auto

#Annotate enhancers to their nearest TSS
#SE
annotatePeaks.pl ${TagDir_ctrl}/superEnhancers.txt mm10 > ${TagDir_ctrl}/ctrl_SE_annotated.txt
#TE
annotatePeaks.pl ${TagDir_ctrl}/typicalEnhancers.txt mm10 > ${TagDir_ctrl}/ctrl_TE_annotated.txt


#anti-IgM
#Get SE and TE
TagDir_IgM="path to your Tag Directory of H3K27Ac anti-IgM 060 min data"
findPeaks ${TagDir_IgM} -style super -typical ${TagDir_IgM}/typicalEnhancers.txt -L 0 -fdr 0.00005 -o auto

#Annotate enhancers to their nearest TSS
#SE
annotatePeaks.pl ${TagDir_IgM}/superEnhancers.txt mm10 > ${TagDir_IgM}/IgM_SE_annotated.txt
#TE
annotatePeaks.pl ${TagDir_IgM}/typicalEnhancers.txt mm10 > ${TagDir_IgM}/IgM_TE_annotated.txt