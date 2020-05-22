library(ggplot2)

tss <- read.table("path to your table of signal intensity in TSSs obtained with Get_signal.sh", 
                  header=T, stringsAsFactors=F, sep="\t", quote="", comment.char = "")
colnames(tss)[20:25] = c("H3K27Ac_0h","H3K27Ac_1h","RelA_0h","RelA_1h","ATAC_0h","ATAC_1h")

tss$H3K27Ac_FC = log2(tss$H3K27Ac_1h/tss$H3K27Ac_0h)
tss$RelA_FC = log2(tss$RelA_1h/tss$RelA_0h)
tss$ATAC_FC = log2(tss$ATAC_1h/tss$ATAC_0h)
tss$Mean_H3K27Ac = (tss$H3K27Ac_0h + tss$H3K27Ac_1h)/2
tss$Mean_ATAC = (tss$ATAC_0h + tss$ATAC_1h)/2

tss = tss[which(!is.na(tss$H3K27Ac_FC) & is.finite(tss$H3K27Ac_FC) & !is.na(tss$ATAC_FC) & is.finite(tss$ATAC_FC)),]

#SFig.3D-E-------------------------------------------------------------------------------------------------
#SFig.3D
g <- ggplot(tss,aes(x = log10(Mean_H3K27Ac)))+
  geom_histogram(bins = 100)+
  theme_bw()+
  theme(panel.grid = element_blank(), strip.background = element_blank(),strip.text = element_blank())+
  labs(x="Mean H3K27Ac signal (log10rpm)", y="Count")+
  theme(axis.title.y = element_text(size = 25))+
  theme(axis.title.x = element_text(size = 25))+
  theme(axis.text.x = element_text(size = 20))+
  theme(axis.text.y = element_text(size = 20))+
  theme(panel.background = element_rect(fill = "white", colour = "black", size = 1.2))+
  geom_vline(xintercept = log10(2.00), linetype="solid")

plot(g)

#SFig.3E
g <- ggplot(tss,aes(x = log10(Mean_ATAC)))+
  geom_histogram(bins = 100)+
  theme_bw()+
  theme(panel.grid = element_blank(), strip.background = element_blank(),strip.text = element_blank())+
  labs(x="Mean ATAC signal (log10rpm)", y="Count")+
  theme(axis.title.y = element_text(size = 25))+
  theme(axis.title.x = element_text(size = 25))+
  theme(axis.text.x = element_text(size = 20))+
  theme(axis.text.y = element_text(size = 20))+
  theme(panel.background = element_rect(fill = "white", colour = "black", size = 1.2))+
  geom_vline(xintercept = log10(3.16), linetype="solid")

plot(g)
#----------------------------------------------------------------------------------------------------------

#Get expressing genes
expressed = tss[which((tss$Mean_H3K27Ac>2.00)&(tss$Mean_ATAC>3.16)),] #threshold: 10^0.3=2.00, 10^0.5=3.16

up = 0.1
dw = -0.1

expressed$type = "XXX"
expressed$type[which((expressed$H3K27Ac_FC > up) & (expressed$ATAC_FC > up))] = "gene_up" #Gained genes
expressed$type[which((expressed$H3K27Ac_FC < dw) & (expressed$ATAC_FC < dw))] = "gene_dw" #Lost genes
expressed$type[which((abs(expressed$H3K27Ac_FC) <= up) & (abs(expressed$ATAC_FC) <= up))] = "gene_cons" #Conserved genes
expressed = expressed[which(expressed$type != "XXX"),]

expressed$TSS = (expressed$Start + expressed$End)/2

#SE assign
SE <- read.csv("path to your SE_3classified.csv obtained using Fig2A-B.R",stringsAsFactors = F)
SE$Center <- (SE$start + SE$end)/2

SE$Target <- "XXX"

expressed_up = expressed[which(expressed$type == "gene_up"),]

expressed_cons = expressed[which(expressed$type == "gene_cons"),]

expressed_dw = expressed[which(expressed$type == "gene_dw"),]

for(i in 1:nrow(SE)){
  if(SE$group[i] == "SE_sg"){
    tmp <- expressed_up[which(expressed_up$Chr == SE$chr[i]),]
    if(nrow(tmp) == 0) next
    SE$Target[i] <- tmp$Gene.Name[which.min(abs(tmp$TSS - SE$Center[i]))]
    
  }else if(SE$group[i] == "SE_sl"){
    tmp <- expressed_dw[which(expressed_dw$Chr == SE$chr[i]),]
    if(nrow(tmp) == 0) next
    SE$Target[i] <- tmp$Gene.Name[which.min(abs(tmp$TSS - SE$Center[i]))]
    
  }else{
    tmp <- expressed_cons[which(expressed_cons$Chr == SE$chr[i]),]
    if(nrow(tmp) == 0) next
    SE$Target[i] <- tmp$Gene.Name[which.min(abs(tmp$TSS - SE$Center[i]))]
  }
}

SE <- SE[which(SE$Target != "XXX"),]

#TE assign
TE <- read.csv("path to your TE_3classified.csv obtained using Fig2A-B.R", stringsAsFactors = F)
TE$Center <- (TE$start + TE$end)/2

TE$Target <- "XXX"

for(i in 1:nrow(TE)){
  if(TE$group[i] == "TE_sg"){
    tmp <- expressed_up[which(expressed_up$Chr == TE$chr[i]),]
    if(nrow(tmp) == 0) next
    TE$Target[i] <- tmp$Gene.Name[which.min(abs(tmp$TSS - TE$Center[i]))]
    
  }else if(TE$group[i] == "TE_sl"){
    tmp <- expressed_dw[which(expressed_dw$Chr == TE$chr[i]),]
    if(nrow(tmp) == 0) next
    TE$Target[i] <- tmp$Gene.Name[which.min(abs(tmp$TSS - TE$Center[i]))]
    
  }else{
    tmp <- expressed_cons[which(expressed_cons$Chr == TE$chr[i]),]
    if(nrow(tmp) == 0) next
    TE$Target[i] <- tmp$Gene.Name[which.min(abs(tmp$TSS - TE$Center[i]))]

  }
}

TE <- TE[which(TE$Target != "XXX"),]

df <- rbind(SE,TE)
write.csv(df,"SE_TE_3classified_assigned.csv",row.names = F) #This file is required for downstream analysis.
