#!/bin/bash

extended_ATAC_IgM_peak="path to your extended bulk ATAC-seq peaks (BED format) in anti-IgM 060 min condition obtained using peak_call_for_bulk.sh and Extend_bulk_ATAC_peaks.R"
RelA_IgM_peak="path to your RelA ChIP-seq peaks (BED format) in anti-IgM 060 min condition obtained using peak_call_for_bulk.sh (RelA_IgM_peaks.bed)"
RelA_in_ATAC_IgM="path to the output file"

intersectBed -a ${extended_ATAC_IgM_peak} -b ${RelA_IgM_peak} > ${RelA_in_ATAC_IgM}