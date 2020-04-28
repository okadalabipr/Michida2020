library(ggplot2)
library(tidyverse)

#This code is assuming you specified Tag Directories in H3K27Ac_ctrl, H3K27Ac_IgM, 
#RelA_ctrl, RelA_IgM, ATAC_ctrl, ATAC_IgM order at annotatePeaks.pl in Get_signal_for_Fig1D.sh.


SE <- read.csv("path to your signal intensity in SE catalog obtained with Get_signal_for_Fig1D.sh", 
               sep = "\t", stringsAsFactors = F)[,c(1:4,20:25,16)]

colnames(SE) <- c("ID","chr","start","end","H3K27Ac_0H","H3K27Ac_1H",
                  "RelA_0H","RelA_1H","ATAC_0H","ATAC_1H","gene_name")

TE <- read.csv("path to your signal intensity in TE catalog obtained with Get_signal_for_Fig1D.sh", 
               sep = "\t", stringsAsFactors = F)[,c(1:4,20:25,16)]

colnames(TE) <- c("ID","chr","start","end","H3K27Ac_0H","H3K27Ac_1H",
                  "RelA_0H","RelA_1H","ATAC_0H","ATAC_1H","gene_name")

SE$H3K27Ac_FC = log2(SE$H3K27Ac_1H/SE$H3K27Ac_0H)
SE$RelA_FC = log2(SE$RelA_1H/SE$RelA_0H)
SE$ATAC_FC = log2(SE$ATAC_1H/SE$ATAC_0H)
SE = dplyr::filter(SE, !is.na(H3K27Ac_FC) & is.finite(H3K27Ac_FC) & !is.na(RelA_FC) 
                   & is.finite(RelA_FC) & !is.na(ATAC_FC) & is.finite(ATAC_FC))
SE$mean_RelA_signal = (SE$RelA_0H + SE$RelA_1H)/2

TE$H3K27Ac_FC = log2(TE$H3K27Ac_1H/TE$H3K27Ac_0H)
TE$RelA_FC = log2(TE$RelA_1H/TE$RelA_0H)
TE$ATAC_FC = log2(TE$ATAC_1H/TE$ATAC_0H)
TE = dplyr::filter(TE, !is.na(H3K27Ac_FC) & is.finite(H3K27Ac_FC) & !is.na(RelA_FC) 
                   & is.finite(RelA_FC) & !is.na(ATAC_FC) & is.finite(ATAC_FC))
TE$mean_RelA_signal = (TE$RelA_0H + TE$RelA_1H)/2

RelA_enriched_TE <- dplyr::filter(TE, mean_RelA_signal >= 1) #Obtain RelA enriched TEs

#H3K27Ac vs RelA (SE)
SE_H3_vs_RelA = ggplot(SE, aes(x = RelA_FC, y = H3K27Ac_FC))+
  geom_point(color = "firebrick2",alpha = 0.5)+
  xlab("Fold Change in RelA Signal (log2)")+
  ylab("Fold Change in H3K27Ac Signal (log2)")+
  theme_bw()+
  theme(axis.text.x = element_text(size = 20),axis.text.y = element_text(size = 20))+
  theme(axis.title.x = element_text(size = 20),axis.title.y = element_text(size = 20))+
  theme(panel.grid = element_blank(), strip.background = element_blank(), strip.text = element_blank())+
  theme(panel.background = element_rect(fill = "white", colour = "black", size = 1.2))

plot(SE_H3_vs_RelA)

cor(SE$H3K27Ac_FC,SE$RelA_FC,method="p") #Peason correlation

#H3K27Ac vs RelA (TE)
TE_H3_vs_RelA = ggplot(TE, aes(x = RelA_FC, y = H3K27Ac_FC))+
  geom_point(color = "royalblue2",alpha = 0.5)+
  xlab("Fold Change in RelA Signal (log2)")+
  ylab("Fold Change in H3K27Ac Signal (log2)")+
  theme_bw()+
  theme(axis.text.x = element_text(size = 20),axis.text.y = element_text(size = 20))+
  theme(axis.title.x = element_text(size = 20),axis.title.y = element_text(size = 20))+
  theme(panel.grid = element_blank(), strip.background = element_blank(), strip.text = element_blank())+
  theme(panel.background = element_rect(fill = "white", colour = "black", size = 1.2))

plot(TE_H3_vs_RelA)

cor(TE$H3K27Ac_FC,TE$RelA_FC,method="p") #Peason correlation

#H3K27Ac vs RelA (RelA enriched TE)
RelA_enriched_TE_H3_vs_RelA = ggplot(RelA_enriched_TE, aes(x = RelA_FC, y = H3K27Ac_FC))+
  geom_point(color = "seagreen3",alpha = 0.5)+
  xlab("Fold Change in RelA Signal (log2)")+
  ylab("Fold Change in H3K27Ac Signal (log2)")+
  theme_bw()+
  theme(axis.text.x = element_text(size = 20),axis.text.y = element_text(size = 20))+
  theme(axis.title.x = element_text(size = 20),axis.title.y = element_text(size = 20))+
  theme(panel.grid = element_blank(), strip.background = element_blank(), strip.text = element_blank())+
  theme(panel.background = element_rect(fill = "white", colour = "black", size = 1.2))

plot(RelA_enriched_TE_H3_vs_RelA)

cor(RelA_enriched_TE$H3K27Ac_FC,RelA_enriched_TE$RelA_FC,method="p") #Peason correlation

#H3K27Ac vs ATAC (SE)
SE_H3_vs_ATAC = ggplot(SE, aes(x = ATAC_FC, y = H3K27Ac_FC))+
  geom_point(color = "firebrick2",alpha = 0.5)+
  xlab("Fold Change in ATAC Signal (log2)")+
  ylab("Fold Change in H3K27Ac Signal (log2)")+
  theme_bw()+
  theme(axis.text.x = element_text(size = 20),axis.text.y = element_text(size = 20))+
  theme(axis.title.x = element_text(size = 20),axis.title.y = element_text(size = 20))+
  theme(panel.grid = element_blank(), strip.background = element_blank(), strip.text = element_blank())+
  theme(panel.background = element_rect(fill = "white", colour = "black", size = 1.2))

plot(SE_H3_vs_ATAC)

cor(SE$H3K27Ac_FC,SE$ATAC_FC,method="p") #Peason correlation

#H3K27Ac vs ATAC (TE)
TE_H3_vs_ATAC = ggplot(TE, aes(x = ATAC_FC, y = H3K27Ac_FC))+
  geom_point(color = "royalblue2",alpha = 0.5)+
  xlab("Fold Change in ATAC Signal (log2)")+
  ylab("Fold Change in H3K27Ac Signal (log2)")+
  theme_bw()+
  theme(axis.text.x = element_text(size = 20),axis.text.y = element_text(size = 20))+
  theme(axis.title.x = element_text(size = 20),axis.title.y = element_text(size = 20))+
  theme(panel.grid = element_blank(), strip.background = element_blank(), strip.text = element_blank())+
  theme(panel.background = element_rect(fill = "white", colour = "black", size = 1.2))

plot(TE_H3_vs_ATAC)

cor(TE$H3K27Ac_FC,TE$ATAC_FC,method="p") #Peason correlation

#H3K27Ac vs ATAC (ATAC enriched TE)
RelA_enriched_TE_H3_vs_ATAC = ggplot(RelA_enriched_TE, aes(x = ATAC_FC, y = H3K27Ac_FC))+
  geom_point(color = "seagreen3",alpha = 0.5)+
  xlab("Fold Change in ATAC Signal (log2)")+
  ylab("Fold Change in H3K27Ac Signal (log2)")+
  theme_bw()+
  theme(axis.text.x = element_text(size = 20),axis.text.y = element_text(size = 20))+
  theme(axis.title.x = element_text(size = 20),axis.title.y = element_text(size = 20))+
  theme(panel.grid = element_blank(), strip.background = element_blank(), strip.text = element_blank())+
  theme(panel.background = element_rect(fill = "white", colour = "black", size = 1.2))

plot(RelA_enriched_TE_H3_vs_ATAC)

cor(RelA_enriched_TE$H3K27Ac_FC,RelA_enriched_TE$ATAC_FC,method="p") #Peason correlation
