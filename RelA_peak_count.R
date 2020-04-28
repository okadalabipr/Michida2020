library(tidyverse)

RelA_peaks <- as.data.frame(read.table("path to your RelA peaks in ATAC peaks obtained using Get_RelA_peaks_in_ATAC_peaks.sh")[,1:3])
colnames(RelA_peaks) <- c("chr","start","end")

Enhancers <- read.csv("path to your SE_TE_3classified_assigned.csv obtained with gene_assignment.R",stringsAsFactors = F)

chr <- c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10",
        "chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19",
        "chrX","chrY")

#Count peaks in signal-gained SE--------------------------------------------------------------------------------------------------
SE_sg <- NULL
for (c in chr){
  En_chr = dplyr::filter(Enhancers, group == "SE_sg" & chr == c)
  if (nrow(En_chr) != 0){
    peak_chr = dplyr::filter(RelA_peaks, chr == c)
    v = c()
    for (i in 1:nrow(En_chr)){
      peak_num = length(which(peak_chr$end >= En_chr$start[i] & En_chr$end[i] >= peak_chr$start))
      v = c(v,peak_num)
    }
    En_chr$peak_num = v
    SE_sg <- rbind(SE_sg,En_chr)
  }
}

SE_sg <- dplyr::filter(SE_sg, peak_num > 0) #Discard peaks which doesn't have RelA peaks.

#Count peaks in signal-gained TE-------------------------------------------------------------------------------------
TE_sg <- NULL
for (c in chr){
  En_chr = dplyr::filter(Enhancers, group == "TE_sg" & chr == c)
  if (nrow(En_chr) != 0){
    peak_chr = dplyr::filter(RelA_peaks, chr == c)
    v = c()
    for (i in 1:nrow(En_chr)){
      peak_num = length(which(peak_chr$end >= En_chr$start[i] & En_chr$end[i] >= peak_chr$start))
      v = c(v,peak_num)
    }
    En_chr$peak_num = v
    TE_sg <- rbind(TE_sg,En_chr)
  }
}

TE_sg <- dplyr::filter(TE_sg, peak_num > 0) #Discard peaks which doesn't have RelA peaks.

#Extract enhancers whose target genes' expression data are exist in RNA-seq data-------------------------------------
#signal-gained SE
RamDA_SE_sg = read.csv("path to your SE_sg_genes.csv obtaind using Fig2C.R",stringsAsFactors = F)

n = c()

for (i in 1:nrow(RamDA_SE_sg)){
  p = which(as.vector(SE_sg$Target) == as.vector(RamDA_SE_sg$gene[i]))
  if (length(p) != 0) n = c(n,p)
}


SE_sg = SE_sg[n,]

colnames(SE_sg)[17] <- "gene"

df_SE_sg <- dplyr::left_join(SE_sg,RamDA_SE_sg,by = "gene")

write.csv(df_SE_sg,"SE_sg_peak_count.csv",row.names = F) #This file is required in downstream analysis.

#signal-gained TE
RamDA_TE_sg = read.csv("path to your TE_sg_genes.csv obtaind using Fig2C.R",stringsAsFactors = F)

n = c()

for (i in 1:nrow(RamDA_TE_sg)){
  p = which(as.vector(TE_sg$Target) == as.vector(RamDA_TE_sg$gene[i]))
  if (length(p) != 0) n = c(n,p)
}

TE_sg = TE_sg[n,]

colnames(TE_sg)[17] <- "gene"

df_TE_sg <- dplyr::left_join(TE_sg,RamDA_TE_sg,by = "gene")

write.csv(df_TE_sg,"TE_sg_peak_count.csv",row.names = F) #This file is required in downstream analysis.
