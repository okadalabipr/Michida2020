#!/bin/bash

SE_ctrl="path to your superEnhancers.txt of anti-IgM 000 min data obtained by findPeaks of HOMER"
TE_ctrl="path to your typicalEnhancers.txt of anti-IgM 000 min data obtained by findPeaks of HOMER"

SE_ctrl_BED="HOMER peak file, superEnhancers.txt (anti-IgM 000 min) is converting into BED format using following command. Please write path to the converted file here."
TE_ctrl_BED="As in SE."

SE_IgM="path to your superEnhancers.txt of anti-IgM 060 min data obtained by findPeaks of HOMER"
TE_IgM="path to your typicalEnhancers.txt of anti-IgM 060 min data obtained by findPeaks of HOMER"

SE_IgM_BED="HOMER peak file, superEnhancers.txt (anti-IgM 060 min) is converting into BED format using following command. Please write path to the converted file here."
TE_IgM_BED="As in SE."

SE_catalog="path to SE catalog you are making with this script"
TE_catalog="path to TE catalog you are making with this script"

#2 tmp BED files are required. Please specify.

tmp1="path to tmp file 1"
tmp2="path to tmp file 2"

#Convert HOMER peak file into BED format
sed '/^#/d' ${SE_ctrl} | awk -v 'OFS=\t' '{print $2, $3, $4}' > ${SE_ctrl_BED}
sed '/^#/d' ${TE_ctrl} | awk -v 'OFS=\t' '{print $2, $3, $4}' > ${TE_ctrl_BED}

sed '/^#/d' ${SE_IgM} | awk -v 'OFS=\t' '{print $2, $3, $4}' > ${SE_IgM_BED}
sed '/^#/d' ${TE_IgM} | awk -v 'OFS=\t' '{print $2, $3, $4}' > ${TE_IgM_BED}

#Make SE catalog
cat ${SE_ctrl_BED} ${SE_IgM_BED} > ${tmp1}
sortBed -i ${tmp1} > ${tmp2}
mergeBed -i ${tmp2} > ${SE_catalog}
rm -rf ${tmp1} ${tmp2}

#Make TE catalog
cat ${TE_ctrl_BED} ${TE_IgM_BED} > ${tmp1}
sortBed -i ${tmp1} > ${tmp2}
mergeBed -i ${tmp2} > ${TE_catalog}
rm -rf ${tmp1} ${tmp2}



