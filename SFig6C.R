library(tidyverse)
library(ggplot2)
library(gridExtra)

TPM <- read.table("path to your cleaned_TPM.txt made by Make_data_table.R", header=T, stringsAsFactors=F,sep="\t")

Irf4 <- as.numeric(TPM[which(TPM$Geneid == "Irf4"),2:477])
Smad3 <- as.numeric(TPM[which(TPM$Geneid == "Smad3"),2:477])
Fam43a <- as.numeric(TPM[which(TPM$Geneid == "Fam43a"),2:477])

Rela <- as.numeric(TPM[which(TPM$Geneid == "Rela"),2:477])
Ran <- as.numeric(TPM[which(TPM$Geneid == "Ran"),2:477])
Eif3d <- as.numeric(TPM[which(TPM$Geneid == "Eif3d"),2:477])

df <- rbind(Irf4,Smad3,Fam43a,Rela,Ran,Eif3d)

plot_df <- NULL

gene <- c("Irf4", "Smad3", "Fam43a", "Rela", "Ran", "Eif3d")

for(i in 1:nrow(df)){
  dose0_mean <- df[i,1:95] %>% as.numeric() %>% mean()
  dose0_sd <- df[i,1:95] %>% as.numeric() %>% sd()
  
  dose0.01_mean <- df[i,96:191] %>% as.numeric() %>% mean()
  dose0.01_sd <- df[i,96:191] %>% as.numeric() %>% sd()
  
  dose0.1_mean <- df[i,192:286] %>% as.numeric() %>% mean()
  dose0.1_sd <- df[i,192:286] %>% as.numeric() %>% sd()
  
  dose1_mean <- df[i,287:381] %>% as.numeric() %>% mean()
  dose1_sd <- df[i,287:381] %>% as.numeric() %>% sd()
  
  dose10_mean <- df[i,382:476] %>% as.numeric() %>% mean()
  dose10_sd <- df[i,382:476] %>% as.numeric() %>% sd()
  
  mean <- c(dose0_mean, dose0.01_mean, dose0.1_mean, dose1_mean, dose10_mean)
  sd <- c(dose0_sd, dose0.01_sd, dose0.1_sd, dose1_sd, dose10_sd)
  
  tmp <- cbind(mean,sd) %>% as.data.frame()
  
  tmp$gene <- gene[i]
  
  plot_df <- rbind(plot_df,tmp)
}

colnames(plot_df)[1:2] <- c("MEAN", "SD")

plot_df$gene <- factor(plot_df$gene, levels=c("Irf4", "Smad3", "Fam43a", "Rela", "Ran", "Eif3d"))

plot_df$CV <- plot_df$SD/plot_df$MEAN

plot_df$tmp <- rep(1:5,6)

plot_df$type <- c(rep("SE",15),rep("TE",15))

#Plot the graph
g2 <- ggplot(plot_df,aes(x = tmp, y = CV, fill = type))+
  geom_bar(stat = "identity")+
  facet_wrap(~gene, ncol = 3)+
  theme_bw()+
  ylab("CV of TPM")+
  xlab("anti-IgM concentration (ug/ml)")+
  scale_x_continuous(breaks = seq(1,5,1),labels=c(0, 0.01, 0.1, 1 ,10))+
  scale_fill_manual(values =c("firebrick3","royalblue"))+
  guides(fill = F)+
  theme(panel.grid = element_blank(), strip.background = element_blank())+
  theme(axis.title.y = element_blank())+
  theme(axis.title.x = element_blank())+
  theme(axis.text.x = element_text(size = 15))+
  theme(axis.text.y = element_text(size = 15))+
  theme(strip.text = element_text(face = "italic", size = 20))+
  theme(panel.background = element_rect(fill = "white", colour = "black", size = 1.2))

plot(g2)