library(ggplot2)
library(tidyverse)

tss <- read.table("path to your table of signal intensity in TSSs obtained with Get_signal.sh", 
                  header=T, stringsAsFactors=F, sep="\t", quote="", comment.char = "") 
colnames(tss)[1] = "PeakID"
colnames(tss)[20:25] = c("H3K27Ac_0h","H3K27Ac_1h","RelA_0h","RelA_1h","ATAC_0h","ATAC_1h")

tss$H3K27Ac_FC = log2(tss$H3K27Ac_1h/tss$H3K27Ac_0h)
tss$RelA_FC = log2(tss$RelA_1h/tss$RelA_0h)
tss$ATAC_FC = log2(tss$ATAC_1h/tss$ATAC_0h)
tss$Mean_H3K27Ac = (tss$H3K27Ac_0h + tss$H3K27Ac_1h)/2
tss$Mean_ATAC = (tss$ATAC_0h + tss$ATAC_1h)/2

tss = tss[which(!is.na(tss$H3K27Ac_FC) & is.finite(tss$H3K27Ac_FC) & !is.na(tss$ATAC_FC) & is.finite(tss$ATAC_FC)),] 


expressed = tss[which( (tss$Mean_H3K27Ac>2.00)&(tss$Mean_ATAC>3.16) ),]

expressed$TSS = (expressed$Start + expressed$End)/2

#assign----------------------------------------------------------------------------------
assign <- function(peak){
  peak <- read.table(peak,
                     header=T, stringsAsFactors=F, sep="\t", quote="", comment.char = "")
  
  colnames(peak)[1] <- "ID"
  peak$Center <- (peak$start + peak$end)/2
  
  peak$Target <- "XXX"
  
  for(i in 1:nrow(peak)){
    tmp <- expressed[which(expressed$Chr == peak$chr[i]),]
    if(nrow(tmp) == 0) next
    peak$Target[i] <- tmp$Gene.Name[which.min(abs(tmp$TSS - peak$Center[i]))]
  }
  
  peak <- peak[which(peak$Target != "XXX"),]
  
  genes <- peak$Target
  
  return(genes)
}

#When you ran Get_venn_num.sh, you got peak files of each venn region made by mergepeaks. 
#Please specify the path of them here in this order.
files <- c("anti-IgM 000 min unique superEnhancers",
           "Sheared superEnhancers",
           "anti-IgM 060 min unique superEnhancers",
           "anti-IgM 000 min unique typicalEnhancers",
           "Sheared typicalEnhancers",
           "anti-IgM 060 min unique typicalEnhancers")

Ns <- c()
genes <- c()

for(f in files){
  tmp <- assign(f)
  genes <- append(genes,tmp)
  Ns <- append(Ns,length(tmp))
}

#Get expression data in anti-IgM 10 ug/ml condition--------------------------------------------------------
RamDA <- read.csv("path to your RamDA.csv obtained with SCDE.R", stringsAsFactors = F)
RamDA <- RamDA[,c(1,15)]
colnames(RamDA) <- c("gene","FC_10")
RamDA$FC_10 <- 2^RamDA$FC_10

get_FC <- function(gene){
  tmp <- which(gene == RamDA$gene)
  if(length(tmp) == 1) return(RamDA$FC_10[tmp])
  else return(NA)
}

FCs <- c()
for(i in 1:length(genes)) FCs <- append(FCs, get_FC(genes[i]))

df <- data.frame(FC_10 = FCs,
                 group = c(rep("Ctrl unique SEs", Ns[1]),
                            rep("Sheared SEs", Ns[2]),
                            rep("anti-IgM unique SEs", Ns[3]),
                            rep("Ctrl unique TEs", Ns[4]),
                            rep("Sheared TEs", Ns[5]),
                            rep("anti-IgM unique TEs", Ns[6])))

df <- df[which(!is.na(df$FC_10)),]

#Plot the graph-------------------------------------------------------------------------------------------------------
stats_se_only_ctrl <- quantile(df$FC_10[grep("Ctrl unique SEs", df$group)],c(1/10,1/4,1/2,3/4,9/10)) %>% as.numeric()
stats_se_overlap <- quantile(df$FC_10[grep("Sheared SEs", df$group)],c(1/10,1/4,1/2,3/4,9/10)) %>% as.numeric()
stats_se_only_IgM <- quantile(df$FC_10[grep("anti-IgM unique SEs", df$group)],c(1/10,1/4,1/2,3/4,9/10)) %>% as.numeric()

stats_te_only_ctrl <- quantile(df$FC_10[grep("Ctrl unique TEs", df$group)],c(1/10,1/4,1/2,3/4,9/10)) %>% as.numeric()
stats_te_overlap <- quantile(df$FC_10[grep("Ctrl unique TEs", df$group)],c(1/10,1/4,1/2,3/4,9/10)) %>% as.numeric()
stats_te_only_IgM <- quantile(df$FC_10[grep("anti-IgM unique TEs", df$group)],c(1/10,1/4,1/2,3/4,9/10)) %>% as.numeric()


plot_df <- data.frame(rbind(stats_se_only_ctrl, stats_se_overlap, stats_se_only_IgM,
                       stats_te_only_ctrl, stats_te_overlap, stats_te_only_IgM))

plot_df$L <- c("SE only in ctrl","Sheared SE", "SE only in anti-IgM", 
          "TE only in ctrl","Sheared TE", "TE only in anti-IgM")

plot_df$L <- factor(plot_df$L, levels = plot_df$L)

plot_df$mean <- c(mean(df$FC_10[grep("Ctrl unique SEs", df$group)]),
             mean(df$FC_10[grep("Sheared SEs", df$group)]), 
             mean(df$FC_10[grep("anti-IgM unique SEs", df$group)]),
             mean(df$FC_10[grep("Ctrl unique TEs", df$group)]),
             mean(df$FC_10[grep("Ctrl unique TEs", df$group)]), 
             mean(df$FC_10[grep("anti-IgM unique TEs", df$group)]))


p <- ggplot(plot_df, aes(x=L, lower=X2, upper=X4, middle=X3, ymin=X1, ymax=X5, fill=L))+
  theme_classic(base_size = 20)+
  theme(legend.position = "none")+
  #scale_fill_manual(values = c("firebrick3","royalblue") ) +
  #scale_y_continuous(limits=c(0,15), breaks=seq(0,15,5))+
  labs(x="", y="FC anti-IgM 0 ƒÊg/ml to 10 ƒÊg/ml")+
  geom_errorbar(width=0.5)+
  theme(axis.text.x = element_text(size = 20))+
  theme(axis.text.y = element_text(size = 20))+
  geom_boxplot(stat="identity")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

p <- p + geom_point(data = plot_df,aes(x = L,y = mean),fill = "white", size = 3, shape = 23)

plot(p)