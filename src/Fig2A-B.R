library(ggplot2)
library(tidyverse)
library(pals)
#This code is assuming you specified Tag Directories in H3K27Ac_ctrl, H3K27Ac_IgM, 
#RelA_ctrl, RelA_IgM, ATAC_ctrl, ATAC_IgM order at annotatePeaks.pl in Get_signal.sh.

SE <- read.csv("path to your signal intensity in SE catalog obtained with Get_signal.sh",
               sep = "\t", stringsAsFactors = F)[,c(1:4,20:25)]

colnames(SE) <- c("ID","chr","start","end","H3K27Ac_0H","H3K27Ac_1H","RelA_0H","RelA_1H","ATAC_0H","ATAC_1H")

TE <- read.csv("path to your signal intensity in TE catalog obtained with Get_signal.sh",
               sep = "\t", stringsAsFactors = F)[,c(1:4,20:25)]

colnames(TE) <- c("ID","chr","start","end","H3K27Ac_0H","H3K27Ac_1H","RelA_0H","RelA_1H","ATAC_0H","ATAC_1H")

SE$H3K27Ac_FC = log2(SE$H3K27Ac_1H/SE$H3K27Ac_0H)
SE$RelA_FC = log2(SE$RelA_1H/SE$RelA_0H)
SE$ATAC_FC = log2(SE$ATAC_1H/SE$ATAC_0H)
SE = dplyr::filter(SE, !is.na(H3K27Ac_FC) & is.finite(H3K27Ac_FC) & !is.na(RelA_FC)
                   & is.finite(RelA_FC) & !is.na(ATAC_FC) & is.finite(ATAC_FC))

TE$H3K27Ac_FC = log2(TE$H3K27Ac_1H/TE$H3K27Ac_0H)
TE$RelA_FC = log2(TE$RelA_1H/TE$RelA_0H)
TE$ATAC_FC = log2(TE$ATAC_1H/TE$ATAC_0H)
TE = dplyr::filter(TE, !is.na(H3K27Ac_FC) & is.finite(H3K27Ac_FC) & !is.na(RelA_FC)
                   & is.finite(RelA_FC) & !is.na(ATAC_FC) & is.finite(ATAC_FC))

SE_sg = dplyr::filter(SE, H3K27Ac_FC > 0.2) #sg: signal gained, su: signal unchanged, sl: siganl lost
SE_su = dplyr::filter(SE, abs(H3K27Ac_FC) <= 0.2)
SE_sl = dplyr::filter(SE, H3K27Ac_FC < -0.2)

TE_sg = dplyr::filter(TE, H3K27Ac_FC > 0.2)
TE_su = dplyr::filter(TE, abs(H3K27Ac_FC) <= 0.2)
TE_sl = dplyr::filter(TE, H3K27Ac_FC < -0.2)

SE_sg$group = rep("SE_sg", nrow(SE_sg))
SE_su$group = rep("SE_su", nrow(SE_su))
SE_sl$group = rep("SE_sl", nrow(SE_sl))

TE_sg$group = rep("TE_sg", nrow(TE_sg))
TE_su$group = rep("TE_su", nrow(TE_su))
TE_sl$group = rep("TE_sl", nrow(TE_sl))

SE = dplyr::bind_rows(SE_sg,SE_su)
SE = dplyr::bind_rows(SE,SE_sl)

TE = dplyr::bind_rows(TE_sg,TE_su)
TE = dplyr::bind_rows(TE,TE_sl)

SE = dplyr::arrange(SE, desc(H3K27Ac_FC)) #Sort in descending order of FC of H3K27Ac
TE = dplyr::arrange(TE, desc(H3K27Ac_FC))

SE$N = seq(1,nrow(SE),1)
TE$N = seq(1,nrow(TE),1)

#These files are needed for gene_assingment.R. Please specify path you like.
write.csv(SE,"SE_3classified.csv", row.names = F)
write.csv(TE,"TE_3classified.csv", row.names = F)

#Fig.2A--------------------------------------------------------------------------------------
#Fig.2A (SE)
SE_g <- ggplot(SE,aes(x = N , y = H3K27Ac_FC, fill = -N))+
  geom_bar(stat = "identity")+
  theme_bw()+
  guides(fill = F)+
  coord_flip()+
  scale_fill_gradientn(colours=coolwarm(10))+
  scale_x_reverse(breaks=c(1,100,200))+
  xlab("Enhancers ranked by H3K27Ac Signal \n Fold Change (log2)")+
  ylab("H3K27Ac Signal Fold Change (log2)")+
  theme(axis.text.x = element_text(size=15),axis.text.y = element_text(size=15))+
  theme(axis.title.x = element_text(size=15),axis.title.y = element_text(size=15))+
  theme(panel.grid = element_blank(), strip.background = element_blank())+
  scale_y_continuous(limits=c(-1,1), breaks=seq(-1,1,0.5))+
  geom_hline(yintercept = c(-0.2,0.2), linetype="dashed", size = 0.8)+
  theme(panel.background = element_rect(fill = "white", colour = "black", size = 2.0))

plot(SE_g)

#Fig.2A (TE)
TE_g <- ggplot(TE,aes(x = N , y = H3K27Ac_FC, fill = -N))+
  geom_bar(stat = "identity")+
  theme_bw()+
  guides(fill = F)+
  coord_flip()+
  scale_fill_gradientn(colours=coolwarm(10))+
  scale_x_reverse(breaks=c(1,2000,4000,6000,8000))+
  xlab("Enhancers ranked by H3K27Ac Signal \n Fold Change (log2)")+
  ylab("H3K27Ac Signal Fold Change (log2)")+
  theme(axis.text.x = element_text(size=15),axis.text.y = element_text(size=15))+
  theme(axis.title.x = element_text(size=15),axis.title.y = element_text(size=15))+
  theme(panel.grid = element_blank(), strip.background = element_blank())+
  scale_y_continuous(limits=c(-6,6), breaks=seq(-6,6,2))+
  geom_hline(yintercept = c(-0.2,0.2), linetype="dashed", size = 0.8)+
  theme(panel.background = element_rect(fill = "white", colour = "black", size = 2.0))

plot(TE_g)

#Fig.2B-------------------------------------------------------------------------------------
#SE
stats_SE_sg <- quantile(SE_sg$RelA_FC,c(1/10,1/4,1/2,3/4,9/10)) %>% as.numeric()
stats_SE_sg <- rbind(stats_SE_sg,quantile(SE_sg$ATAC_FC,c(1/10,1/4,1/2,3/4,9/10)) %>% as.numeric())

stats_SE_su <- quantile(SE_su$RelA_FC,c(1/10,1/4,1/2,3/4,9/10)) %>% as.numeric()
stats_SE_su <- rbind(stats_SE_su,quantile(SE_su$ATAC_FC,c(1/10,1/4,1/2,3/4,9/10)) %>% as.numeric())

stats_SE_sl <- quantile(SE_sl$RelA_FC,c(1/10,1/4,1/2,3/4,9/10)) %>% as.numeric()
stats_SE_sl <- rbind(stats_SE_sl,quantile(SE_sl$ATAC_FC,c(1/10,1/4,1/2,3/4,9/10)) %>% as.numeric())

#Fig.2B RelA SE
df <- data.frame((rbind(stats_SE_sg[1,],stats_SE_su[1,],stats_SE_sl[1,])))
df$L <- c("Gained","Unchanged", "Lost")
df$L <- factor(df$L, levels=c("Gained","Unchanged", "Lost"))

p <- ggplot(df, aes(x=L, lower=X2, upper=X4, middle=X3, ymin=X1, ymax=X5, fill=L))+
  theme_classic(base_size = 20)+
  theme(legend.position = "none")+
  scale_fill_manual(values = c("firebrick3","gray","royalblue")) +
  labs(x="", y="log2 FC in RelA Signal")+
  theme(axis.title.x = element_blank())+
  theme(axis.text.x = element_text(angle = 30, vjust = 0.7,size = 20))+
  geom_errorbar(width=0.5)+
  geom_boxplot(stat = "identity")

plot(p)

#Fig.2B ATAC SE
df <- data.frame((rbind(stats_SE_sg[2,],stats_SE_su[2,],stats_SE_sl[2,])))
df$L <- c("Gained","Unchanged", "Lost")

df$L <- factor(df$L, levels=c("Gained","Unchanged", "Lost"))


p <- ggplot(df, aes(x=L, lower=X2, upper=X4, middle=X3, ymin=X1, ymax=X5, fill=L))+
  theme_classic(base_size = 20)+
  theme(legend.position = "none")+
  scale_fill_manual(values = c("firebrick3","gray","royalblue")) +
  labs(x="", y="log2 FC in ATAC Signal")+
  theme(axis.title.x = element_blank())+
  theme(axis.text.x = element_text(angle = 30, vjust = 0.7,size = 20))+
  geom_errorbar(width=0.5)+
  geom_boxplot(stat = "identity")

plot(p)

#TE
stats_TE_sg <- quantile(TE_sg$RelA_FC,c(1/10,1/4,1/2,3/4,9/10)) %>% as.numeric()
stats_TE_sg <- rbind(stats_TE_sg,quantile(TE_sg$ATAC_FC,c(1/10,1/4,1/2,3/4,9/10)) %>% as.numeric())

stats_TE_su <- quantile(TE_su$RelA_FC,c(1/10,1/4,1/2,3/4,9/10)) %>% as.numeric()
stats_TE_su <- rbind(stats_TE_su,quantile(TE_su$ATAC_FC,c(1/10,1/4,1/2,3/4,9/10)) %>% as.numeric())

stats_TE_sl <- quantile(TE_sl$RelA_FC,c(1/10,1/4,1/2,3/4,9/10)) %>% as.numeric()
stats_TE_sl <- rbind(stats_TE_sl,quantile(TE_sl$ATAC_FC,c(1/10,1/4,1/2,3/4,9/10)) %>% as.numeric())

#Fig.2B RelA TE
df <- data.frame((rbind(stats_TE_sg[1,],stats_TE_su[1,],stats_TE_sl[1,])))
df$L <- c("Gained","Unchanged", "Lost")

df$L <- factor(df$L, levels=c("Gained","Unchanged", "Lost"))

p <- ggplot(df, aes(x=L, lower=X2, upper=X4, middle=X3, ymin=X1, ymax=X5, fill=L))+
  theme_classic(base_size = 20)+
  theme(legend.position = "none")+
  scale_fill_manual(values = c("firebrick3","gray","royalblue")) +
  labs(x="", y="log2 FC in RelA Signal")+
  theme(axis.title.x = element_blank())+
  theme(axis.text.x = element_text(angle = 30, vjust = 0.7,size = 20))+
  geom_errorbar(width=0.5)+
  geom_boxplot(stat="identity")

plot(p)

#Fig.2B ATAC TE
df <- data.frame((rbind(stats_TE_sg[2,],stats_TE_su[2,],stats_TE_sl[2,])))
df$L <- c("Gained","Unchanged", "Lost")

df$L <- factor(df$L, levels=c("Gained","Unchanged", "Lost"))

p <- ggplot(df, aes(x=L, lower=X2, upper=X4, middle=X3, ymin=X1, ymax=X5, fill=L))+
  theme_classic(base_size = 20)+
  theme(legend.position = "none")+
  scale_fill_manual(values = c("firebrick3","gray","royalblue")) +
  labs(x="", y="log2 FC in ATAC Signal")+
  theme(axis.title.x = element_blank())+
  theme(axis.text.x = element_text(angle = 30, vjust = 0.7,size = 20))+
  geom_errorbar(width=0.5)+
  geom_boxplot(stat = "identity")

plot(p)

#SFig.3B-C--------------------------------------------------------------------------------------
make_rand <- function(x){       #This function non-restoring extracts random indexies.
  rand_list <- list(
    sample(1:nrow(SE_sg),x),
    sample(1:nrow(SE_su),x),
    sample(1:nrow(SE_sl),x),
    sample(1:nrow(TE_sg),x),
    sample(1:nrow(TE_su),x),
    sample(1:nrow(TE_sl),x)
  )
  return(rand_list)
}


#SFig.3B
SE_mean = c()
SE_sd = c()
TE_mean = c()
TE_sd = c()

for(i in 5:22){
  SE_tmp = c()
  TE_tmp = c()
  for(j in 1:100){
    rand <- make_rand(i) #Get i random indexies
    SE_vx <- c(SE_sg$RelA_FC[rand[[1]]],SE_su$RelA_FC[rand[[2]]],SE_sl$RelA_FC[rand[[3]]]) #Extract data
    SE_fx <- factor(c(rep("SE_sg",i),rep("SE_su",i),rep("SE_sl",i)))
    SE_anova <- anova(aov(SE_vx~SE_fx)) #Perform One-way ANOVA
    SE_tmp <- c(SE_tmp,SE_anova$`Pr(>F)`[1]) #Get p-value
    
    TE_vx <- c(TE_sg$RelA_FC[rand[[4]]],TE_su$RelA_FC[rand[[5]]],TE_sl$RelA_FC[rand[[6]]])
    TE_fx <- factor(c(rep("TE_sg",i),rep("TE_su",i),rep("TE_sl",i)))
    TE_anova <- anova(aov(TE_vx~TE_fx))
    TE_tmp <- c(TE_tmp,TE_anova$`Pr(>F)`[1])
  }
  
  SE_mean = c(SE_mean,mean(SE_tmp)) #Calculate average of p-values
  SE_sd = c(SE_sd,sd(SE_tmp)) #Calculate SD of p-values
  
  TE_mean = c(TE_mean,mean(TE_tmp))
  TE_sd = c(TE_sd,sd(TE_tmp))
}

pvalue_df <- data.frame(p = c(SE_mean,TE_mean),sn = c(5:22,5:22),
                        label = c(rep("SE",18),rep("TE",18)),error = c(SE_sd,TE_sd))

pvalue_df$p <- pvalue_df$p


p <- ggplot(pvalue_df,aes(x = sn,y = p,color = label))+                        
  geom_point()+
  geom_line()+
  scale_y_continuous(breaks = seq(0,0.6,0.1))+
  #geom_errorbar(aes(ymin = p - error, ymax = p + error), width = 0.5)+
  guides(color = F)+
  theme_bw()+
  xlab("Sampling size")+
  ylab("Mean p-value (one-way anova)")+
  theme(axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size=20))+
  theme(axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20))+
  scale_color_manual(values =c("firebrick3","royalblue"))+
  theme(panel.grid = element_blank())+
  theme(panel.background = element_rect(fill = "white", colour = "black", size = 1.2))+
  geom_hline(yintercept = 0.01, linetype="solid")


plot(p)


#SFig.3C
SE_mean = c()
SE_sd = c()
TE_mean = c()
TE_sd = c()

for(i in 5:22){
  SE_tmp = c()
  TE_tmp = c()
  for(j in 1:100){
    rand <- make_rand(i) #Get i random indexies
    SE_vx <- c(SE_sg$ATAC_FC[rand[[1]]],SE_su$ATAC_FC[rand[[2]]],SE_sl$ATAC_FC[rand[[3]]]) #Extract data
    SE_fx <- factor(c(rep("SE_sg",i),rep("SE_su",i),rep("SE_sl",i)))
    SE_anova <- anova(aov(SE_vx~SE_fx)) #Perform One-way ANOVA
    SE_tmp <- c(SE_tmp,SE_anova$`Pr(>F)`[1]) #Get p-value
    
    TE_vx <- c(TE_sg$ATAC_FC[rand[[4]]],TE_su$ATAC_FC[rand[[5]]],TE_sl$ATAC_FC[rand[[6]]])
    TE_fx <- factor(c(rep("TE_sg",i),rep("TE_su",i),rep("TE_sl",i)))
    TE_anova <- anova(aov(TE_vx~TE_fx))
    TE_tmp <- c(TE_tmp,TE_anova$`Pr(>F)`[1])
  }
  
  SE_mean = c(SE_mean,mean(SE_tmp)) #Calculate average of p-values
  SE_sd = c(SE_sd,sd(SE_tmp)) #Calculate SD of p-values
  
  TE_mean = c(TE_mean,mean(TE_tmp))
  TE_sd = c(TE_sd,sd(TE_tmp))
}

pvalue_df <- data.frame(p = c(SE_mean,TE_mean),sn = c(5:22,5:22),
                        label = c(rep("SE",18),rep("TE",18)),error = c(SE_sd,TE_sd))

pvalue_df$p <- pvalue_df$p


p <- ggplot(pvalue_df,aes(x = sn,y = p,color = label))+                        
  geom_point()+
  geom_line()+
  scale_y_continuous(breaks = seq(0,0.6,0.1))+
  #geom_errorbar(aes(ymin = p - error, ymax = p + error), width = 0.5)+
  guides(color = F)+
  theme_bw()+
  xlab("Sampling size")+
  ylab("Mean p-value (one-way anova)")+
  theme(axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size=20))+
  theme(axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20))+
  scale_color_manual(values =c("firebrick3","royalblue"))+
  theme(panel.grid = element_blank())+
  theme(panel.background = element_rect(fill = "white", colour = "black", size = 1.2))+
  geom_hline(yintercept = 0.01, linetype="solid")


plot(p)

#I'm sorry for a bit messy source code...




