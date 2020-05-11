# Michida et al., 2020 <br>The number of transcription factors at an enhancer determines switch-like gene expression (To appear)


## Introduction

This repository contains the source code for the sequence analysis used in the above paper. 

Relative paths under “src” are shown.

### ATAC sequencing and genomic alignment

- Primary_analysis/ATAC-seq_mapping.sh

### ChIP sequencing and genomic alignment

- Primary_analysis/ChIP-seq_mapping.sh

### Single cell RNA sequencing and genomic alignment

- Primary_analysis/scRNA-seq_mapping.sh
- Primary_analysis/featureCounts.sh
- Fig2/SCDE.R

### Identification of SEs and TEs 

- Primary_analysis/Call_SE_TE.sh
- Fig1/Fig1A.R
- Fig1/Get_venn_num.sh
- Fig1/Fig1B.R
- Primary_analysis/Make_enhancer_catalog.sh
- Fig1/Fig1D.R

### ATAC-seq coverage in SEs and TEs

- Sup_figs/SFig1D.sh
- Sup_figs/SFig1D.R

### Quantification of ChIP and ATAC signals in enhancer and proximal TSS regions 

- Primary_analysis/Make_bigwig.sh
- Primary_analysis/multiBigwigSummary.sh
- Fig1/Fig1C_E-H.R

### Gene assignment to enhancers 

- Fig2/Gene_assignment.R

### Identification of RelA binding numbers from sequence data

- Primary_analysis/Peak_call_for_bulk.sh
- Fig3/Extend_bulk_ATAC_peaks.R
- Fig3/Get_RelA_peaks_in_ATAC_peaks.sh
- Fig3/RelA_peak_motif_count.R

### Motif analysis 

- Fig4/Fig4A-D.sh
- Fig4/Fig4E.R

### Identification of NF-κB and PU.1 motifs

- Fig4/Peak_call_for_scATAC.sh
- Fig4/Extend_sc-ATAC_peaks.R
- Fig4/intersectBed_get_PU1_position_for_scATAC.sh
- Fig4/Count_RelA_PU1_overlap_for_scATAC.sh
- Fig4/Fig4E.R
- Fig5/Get_PU1_motifs_in_enhancers.sh
- Fig5/Get_RelA_PU1_overlap.sh
- Fig3/Get_NFkB_motifs_in_enhancers.sh
- Fig3/RelA_peak_motif_count.R
- Sup_figs/SFig6A.R

### Gene Ontology (GO) analysis 

- Sup_figs/SFig1A.R

### Mathematical modeling 

- Fig3/Filtering_genes.R
- Fig3/Fig3B_D.R
- Sup_figs/SFig5B.R
- Fig5/Fig5A-E.R

### Prediction of single cell mRNA distribution using a mathematical model 

- Fig3/Parameter_estimation.py

## Authors

Hiroki Michida and Hiroaki Imoto

## License

[MIT](LICENSE)
