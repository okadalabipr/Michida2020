library(tidyverse)
library(ggplot2)

motif_in_peak <- read.table("path to ${PU1_in_ATAC_in_SE_TE_IgM} made by Get_RelA_PU1_overlap.sh")
colnames(motif_in_peak) <- c("chr","start","end")

#These two files are obtained using Filtering_genes.R. Please specify their path depending on your environment.
SE <- read.csv("filtered_Gained_SE.csv",stringsAsFactors = F)
TE <- read.csv("filtered_Gained_TE.csv",stringsAsFactors = F)

Enhancers <- rbind(SE,TE)

chr = c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10",
        "chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19",
        "chrX","chrY")

#Count RelA-PU1 overlaps in gained SEs and gained TEs------------------------------------------------------
#Gained SE
SE_sg <- NULL
for (c in chr){
  En_chr = dplyr::filter(Enhancers, group == "SE_sg" & chr == c)
  if (nrow(En_chr) != 0){
    peak_chr = dplyr::filter(motif_in_peak, chr == c)
    v = c()
    for (i in 1:nrow(En_chr)){
      peak_num = length(which(peak_chr$end >= En_chr$start[i] & En_chr$end[i] >= peak_chr$start))
      v = c(v,peak_num)
    }
    En_chr$PU1_num = v
    SE_sg <- rbind(SE_sg,En_chr)
  }
}

#Gained TE
TE_sg <- NULL
for (c in chr){
  En_chr = dplyr::filter(Enhancers, group == "TE_sg" & chr == c)
  if (nrow(En_chr) != 0){
    peak_chr = dplyr::filter(motif_in_peak, chr == c)
    v = c()
    for (i in 1:nrow(En_chr)){
      peak_num = length(which(peak_chr$end >= En_chr$start[i] & En_chr$end[i] >= peak_chr$start))
      v = c(v,peak_num)
    }
    En_chr$PU1_num = v
    TE_sg <- rbind(TE_sg,En_chr)
  }
}

stats_se <- quantile(SE_sg$PU1_num,c(1/10,1/4,1/2,3/4,9/10)) %>% as.numeric()
stats_te <- quantile(TE_sg$PU1_num,c(1/10,1/4,1/2,3/4,9/10)) %>% as.numeric()


df <- data.frame(rbind(stats_se,stats_te))
df$L <- c("Gained SE","Gained TE")

df$mean <- c(mean(SE_sg$PU1_num),mean(TE_sg$PU1_num))

#SFig.6A
p <- ggplot(df, aes(x=L, lower=X2, upper=X4, middle=X3, ymin=X1, ymax=X5, fill=L))+
  theme_classic(base_size = 20)+
  theme(legend.position = "none")+
  scale_fill_manual(values = c("firebrick3","royalblue") ) +
  #scale_y_continuous(limits=c(0,15), breaks=seq(0,15,5))+
  labs(x="", y="Number of PU.1 motifs")+
  geom_errorbar(width=0.5)+
  theme(axis.text.x = element_text(size = 20))+
  theme(axis.text.y = element_text(size = 20))+
  geom_boxplot(stat="identity")

p <- p + geom_point(data = df,aes(x = L,y = mean),fill = "white", size = 3, shape = 23)

plot(p)

t.test(SE_sg$PU1_num,TE_sg$PU1_num,var.equal = F)$p.value #t-test

df$mean #Mean values

#non-zero percentage
100*(which(SE_sg$PU1_num > 0) %>% length())/nrow(SE_sg)
100*(which(TE_sg$PU1_num > 0) %>% length())/nrow(TE_sg)