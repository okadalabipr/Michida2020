library(ggplot2)
library(tidyverse)

RNA_seq <- read.csv("path to your RamDA.csv obtained with SCDE.R",stringsAsFactors = F)

colnames(RNA_seq)[1] <- "gene"

assigned <- read.csv("path to your SE_TE_3classified_assigned.csv obtained with gene_assignment.R", 
                     stringsAsFactors = F)

assigned$RelA_mean <- (assigned$RelA_0H + assigned$RelA_1H)/2

#Replace label, "group" of RelA enriched TEs
#"RE" means "RelA enriched"
for(i in (length(which(grepl("SE",assigned$group))) + 1): nrow(assigned)){
  if(assigned$RelA_mean[i] >= 1){                                       
    if(assigned$group[i] == "TE_sg") assigned$group[i] = "RE_TE_sg"
    else if(assigned$group[i] == "TE_su") assigned$group[i] = "RE_TE_su"
    else assigned$group[i] = "RE_TE_sl"
  }
}

#Assign RNA-seq results of target genes of each enhancer.
v = c()

for (i in 1:nrow(RNA_seq)){
  gyou = match(RNA_seq$gene[i],assigned$Target)
  v <- c(v, as.character(assigned$group[gyou]))
}

RNA_seq$type <- v

RNA_seq <- dplyr::filter(RNA_seq, !is.na(type))

RNA_seq <- RNA_seq[,c(1,3,7,11,15,18)]

#Save genes which were assigned to signal gained enhancers (Used later analysis)
a <- RNA_seq[which(grepl("SE_sg",RNA_seq$type)),]
write.csv(a,"SE_sg_genes.csv", row.names = F)

b <- RNA_seq[which(grepl("TE_sg",RNA_seq$type)),]
write.csv(b,"TE_sg_genes.csv", row.names = F)

#Calculate mean values
#sg: signal gained, su: signal unchanged, sl: siganl lost
log2FC_mean_SE_sg <- c(0,apply(RNA_seq[which(grepl("SE_sg",RNA_seq$type)),2:5],2,mean))

log2FC_mean_SE_su <- c(0,apply(RNA_seq[which(grepl("SE_su",RNA_seq$type)),2:5],2,mean))

log2FC_mean_SE_sl <- c(0,apply(RNA_seq[which(grepl("SE_sl",RNA_seq$type)),2:5],2,mean))

log2FC_mean_TE_sg <- c(0,apply(RNA_seq[which(grepl("TE_sg",RNA_seq$type)),2:5],2,mean))

log2FC_mean_TE_su <- c(0,apply(RNA_seq[which(grepl("TE_su",RNA_seq$type)),2:5],2,mean))

log2FC_mean_TE_sl <- c(0,apply(RNA_seq[which(grepl("TE_sl",RNA_seq$type)),2:5],2,mean))

log2FC_mean_RE_TE_sg <- c(0,apply(RNA_seq[which(grepl("RE_TE_sg",RNA_seq$type)),2:5],2,mean))

log2FC_mean_RE_TE_su <- c(0,apply(RNA_seq[which(grepl("RE_TE_su",RNA_seq$type)),2:5],2,mean))

log2FC_mean_RE_TE_sl <- c(0,apply(RNA_seq[which(grepl("RE_TE_sl",RNA_seq$type)),2:5],2,mean))

#Calculate SD values
log2FC_sd_SE_sg <- c(0,apply(RNA_seq[which(grepl("SE_sg",RNA_seq$type)),2:5],2,sd))

log2FC_sd_SE_su <- c(0,apply(RNA_seq[which(grepl("SE_su",RNA_seq$type)),2:5],2,sd))

log2FC_sd_SE_sl <- c(0,apply(RNA_seq[which(grepl("SE_sl",RNA_seq$type)),2:5],2,sd))

log2FC_sd_TE_sg <- c(0,apply(RNA_seq[which(grepl("TE_sg",RNA_seq$type)),2:5],2,sd))

log2FC_sd_TE_su <- c(0,apply(RNA_seq[which(grepl("TE_su",RNA_seq$type)),2:5],2,sd))

log2FC_sd_TE_sl <- c(0,apply(RNA_seq[which(grepl("TE_sl",RNA_seq$type)),2:5],2,sd))

log2FC_sd_RE_TE_sg <- c(0,apply(RNA_seq[which(grepl("RE_TE_sg",RNA_seq$type)),2:5],2,sd))

log2FC_sd_RE_TE_su <- c(0,apply(RNA_seq[which(grepl("RE_TE_su",RNA_seq$type)),2:5],2,sd))

log2FC_sd_RE_TE_sl <- c(0,apply(RNA_seq[which(grepl("RE_TE_sl",RNA_seq$type)),2:5],2,sd))

plot_df <- data.frame(dose = rep(c(1,2,3,4,5),9), FC = c(log2FC_mean_SE_sg,log2FC_mean_SE_su,log2FC_mean_SE_sl,
                                                         log2FC_mean_TE_sg,log2FC_mean_TE_su,log2FC_mean_TE_sl,
                                                         log2FC_mean_RE_TE_sg,log2FC_mean_RE_TE_su,log2FC_mean_RE_TE_sl),
                      sd = c(log2FC_sd_SE_sg,log2FC_sd_SE_su,log2FC_sd_SE_sl,
                             log2FC_sd_TE_sg,log2FC_sd_TE_su,log2FC_sd_TE_sl,
                             log2FC_sd_RE_TE_sg,log2FC_sd_RE_TE_su,log2FC_sd_RE_TE_sl
                      ),
                      class <- c(rep("sg",5),rep("su",5),rep("sl",5),
                                 rep("sg",5),rep("su",5),rep("sl",5),
                                 rep("sg",5),rep("su",5),rep("sl",5)),
                      type <- c(rep("SE",15),rep("TE",15),rep("RelA enriched TE",15)))

colnames(plot_df) <- c("dose","FC","sd","class","type")

plot_df$type <- factor(plot_df$type, levels = c("SE","TE","RelA enriched TE"))

g <- ggplot(plot_df,aes(x = dose,y = FC, color = class))+
  geom_line(size = 1)+
  geom_point(size = 2)+
  geom_errorbar(aes(ymin = FC - sd, ymax = FC + sd),width = 0.3,alpha = 0.5)+
  facet_wrap(~type,ncol = 3)+
  theme_bw()+
  xlab("anti-IgM concentration (ug/ml)")+
  ylab("Fold change in \nmRNA level (log2)")+
  scale_color_manual(values =c("firebrick3","royalblue","gray"))+
  guides(color = F)+
  theme(axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size=20))+
  theme(axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20))+
  theme(panel.grid = element_blank(), strip.background = element_blank(), 
        strip.text = element_text(size = 20))+
  theme(panel.background = element_rect(fill = "white", colour = "black", size = 1.2))+
  scale_x_continuous(breaks = seq(1,5,1),labels=c(0,0.01,0.1,1,10))+
  geom_hline(yintercept = 0, linetype="solid")

plot(g)
