library(ggplot2)
library(ggrepel)

header <- c("PeakID", "chr", "start", "end", "strand", "Normalized_Tag_Count", "superEnhancer_slope", 
            "findPeaks_Score", "Clonal_Fold_Change")

SE <- read.table("superEnhancers.txt", header = F) #Output of HOMER findPeaks
colnames(SE) <- header

SE_anno <- read.csv("ctrl_SE_annotated.txt", header = T, stringsAsFactors = F, sep = "\t") #Output of HOMER annotatePeaks.pl
colnames(SE_anno)[1] <- "PeakID"

TE <- read.table("typicalEnhancer.txt", header = F) #Output of HOMER findPeaks
colnames(TE) <- header

TE_anno <- read.csv("ctrl_TE_annotated.txt", header = T, stringsAsFactors = F, sep = "\t") #Output of HOMER annotatePeaks.pl
colnames(TE_anno)[1] <- "PeakID"

df <- rbind(SE,TE)
df$rank <- seq(1,nrow(df),1)
df$gene_name <- rep("",nrow(df))
df$type <- c(rep("SE",nrow(SE)),rep("TE",nrow(TE)))

df_anno <- rbind(SE_anno,TE_anno)


#Assign labels to the peak whose rank is the highest among the peaks assinged to the same gene.
selected_genes = c("Irf4", "Tnfaip3", "Cd83", "Cr2", "Nfkbia","Cd44", "Rel")

selected_row <- c()
for(i in 1:length(selected_genes)){
  tmp <- which(selected_genes[i] == df_anno$Gene.Name)[1]
  tmp <- which(df_anno$PeakID[tmp] == df$PeakID)
  selected_row <- append(selected_row,tmp)
}

df$gene_name[selected_row] <- selected_genes


#Plot the graph
g <- ggplot(df, aes(x = rank, y = Normalized_Tag_Count/10, color = type, label = gene_name))+
  geom_point()+
  theme_bw()+
  xlab("Enhancers ranked by H3K27Ac signal")+
  ylab("H3K27Ac ChIP signal (rpm)")+
  scale_x_reverse(limits = c(4500,1), breaks = c(seq(4000,1000,-1000),1))+
  scale_color_manual(values =c("firebrick2","royalblue2"))+
  labs(color = "Enhancer type")+
  geom_text_repel(na.rm = TRUE, size = 7.0, nudge_x = -2000, nudge_y = 200, max.iter = 2000,segment.alpha = 0.5)+
  guides(color = F)+
  theme(axis.text.x = element_text(size=20),axis.text.y = element_text(size=15))+
  theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20))+
  theme(panel.grid = element_blank(), strip.background = element_blank(), strip.text = element_blank())+
  theme(panel.background = element_rect(fill = "white", colour = "black", size = 1.2))

plot(g)