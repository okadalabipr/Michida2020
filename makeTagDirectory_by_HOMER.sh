#!/bin/bash

#H3K27Ac anti-IgM 000 min
H3K27Ac_IgM_000min_bam="path to your bam file of H3K27Ac anti-IgM 000 min data"
H3K27Ac_IgM_000min_TagDir="path to tag directory of H3K27Ac anti-IgM 000 min data you are making by this script"

makeTagDirectory ${H3K27Ac_IgM_000min_TagDir} ${H3K27Ac_No1_IgM_000min_bam} -single -tbp 1 -fragLength 200

#H3K27Ac anti-IgM 060 min
H3K27Ac_IgM_060min_bam="path to your bam file of H3K27Ac anti-IgM 060 min data"
H3K27Ac_IgM_060min_TagDir="path to tag directory of H3K27Ac anti-IgM 060 min data you are making by this script"

makeTagDirectory ${H3K27Ac_IgM_060min_TagDir} ${H3K27Ac_No1_IgM_060min_bam} -single -tbp 1 -fragLength 200

#RelA anti-IgM 000 min
RelA_IgM_000min_bam="path to your bam file of RelA anti-IgM 000 min data"
RelA_IgM_000min_TagDir="path to tag directory of RelA anti-IgM 000 min data you are making by this script"

makeTagDirectory ${RelA_IgM_000min_TagDir} ${RelA_IgM_000min_bam} -single -tbp 1 -fragLength 200

#RelA anti-IgM 060 min
RelA_IgM_060min_bam="path to your bam file of RelA anti-IgM 000 min data"
RelA_IgM_060min_TagDir="path to tag directory of RelA anti-IgM 000 min data you are making by this script"

makeTagDirectory ${RelA_IgM_000min_TagDir} ${RelA_IgM_000min_bam} -single -tbp 1 -fragLength 200

#ATAC anti-IgM 000 min
ATAC_IgM_000min_bam="path to your bam file of ATAC anti-IgM 000 min data"
ATAC_IgM_000min_TagDir="path to tag directory of ATAC anti-IgM 000 min data you are making by this script"

makeTagDirectory ${ATAC_IgM_000min_TagDir} ${ATAC_IgM_000min_bam} -single -tbp 1 -fragLength 200

#ATAC anti-IgM 060 min
ATAC_IgM_060min_bam="path to your bam file of ATAC anti-IgM 060 min data"
ATAC_IgM_060min_TagDir="path to tag directory of ATAC anti-IgM 060 min data you are making by this script"

makeTagDirectory ${ATAC_IgM_060min_TagDir} ${ATAC_IgM_060min_bam} -single -tbp 1 -fragLength 200