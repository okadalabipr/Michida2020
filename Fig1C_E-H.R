library(ggplot2)

#This function is assuming you specified bigwig files in "H3K27Ac", "RelA", "ATAC" order at multiBigwigSummary.
getdata <- function(n) {
  labels <- c("Chr","Start","End","H3K27Ac","RelA","ATAC")
  d <- read.table(n,header = F, sep="\t")
  colnames(d) <- labels
  d$W <- abs(d$Start-d$End)
  return(d)
}

SE <- getdata("path to your tab file of SEs obtained using multiBigwigSummary")
TE <- getdata("path to your tab file of TEs obtained using multiBigwigSummary")
TSS <- getdata("path to your tab file of TSSs obtained using multiBigwigSummary")

#Fig.1c-------------------------------------------------------------------------------------------
stats_SE <- quantile(SE$RelA*SE$W,c(1/10,1/4,1/2,3/4,9/10)) %>% as.numeric()
stats_SE <- rbind(stats_SE,quantile(SE$ATAC*SE$W,c(1/10,1/4,1/2,3/4,9/10)) %>% as.numeric())

stats_TE <- quantile(TE$RelA*TE$W,c(1/10,1/4,1/2,3/4,9/10)) %>% as.numeric()
stats_TE <- rbind(stats_TE,quantile(TE$ATAC*TE$W,c(1/10,1/4,1/2,3/4,9/10)) %>% as.numeric())

stats_TSS <- quantile(TSS$RelA*TSS$W,c(1/10,1/4,1/2,3/4,9/10)) %>% as.numeric()
stats_TSS <- rbind(stats_TSS,quantile(TSS$ATAC*TSS$W,c(1/10,1/4,1/2,3/4,9/10)) %>% as.numeric())

#Fig.1C RelA
df <- data.frame((rbind(stats_SE[1,],stats_TE[1,],stats_TSS[1,]))/1000)
df$L <- c("SE","TE", "TSS")
p <- ggplot(df, aes(x=L, lower=X2, upper=X4, middle=X3, ymin=X1, ymax=X5, fill=L))+
  theme_classic(base_size = 20)+
  theme(legend.position = "none")+
  scale_fill_manual(values = c("firebrick2","royalblue2","seagreen3") ) +
  labs(x="", y="Norm. RelA signal")+
  geom_errorbar(width=0.5 )+
  geom_boxplot(stat="identity")

plot(p)

#One-way ANOVA test with under random sampling
index_TE <- sample(1:nrow(TE),nrow(SE))
index_TSS <- sample(1:nrow(TSS),nrow(SE))
vx <- c(SE$RelA*SE$W,TE$RelA[index_TE]*TE$W[index_TE],TSS$RelA[index_TSS]*TSS$W[index_TSS])
fx <- factor(c(rep("SE",nrow(SE)),rep("TE",nrow(SE)),rep("TSS",nrow(SE))))
anova(aov(vx~fx))


#Fig.1C ATAC
df <- data.frame((rbind(stats_SE[2,],stats_TE[2,],stats_TSS[2,]))/1000)
df$L <- c("SE","TE", "TSS")
p <- ggplot(df, aes(x=L, lower=X2, upper=X4, middle=X3, ymin=X1, ymax=X5, fill=L))+
  theme_classic(base_size = 20)+
  theme(legend.position = "none")+
  scale_fill_manual(values = c("firebrick2","royalblue2","seagreen3") ) +
  labs(x="", y="Norm. ATAC signal")+
  geom_errorbar(width=0.5 )+
  geom_boxplot(stat="identity")

plot(p)

#One-way ANOVA test with under random sampling
index_TE <- sample(1:nrow(TE),nrow(SE))
index_TSS <- sample(1:nrow(TSS),nrow(SE))
vx <- c(SE$ATAC*SE$W,TE$ATAC[index_TE]*TE$W[index_TE],TSS$ATAC[index_TSS]*TSS$W[index_TSS])
fx <- factor(c(rep("SE",nrow(SE)),rep("TE",nrow(SE)),rep("TSS",nrow(SE))))
anova(aov(vx~fx))

#Fig.1E-H----------------------------------------------------------------------------------------
H <- c(SE$H3K27Ac*SE$W,TE$H3K27Ac*TE$W)
R <- c(SE$RelA*SE$W,TE$RelA*TE$W)
A <- c(SE$ATAC*SE$W,TE$ATAC*TE$W)
w  <- c(SE$W, TE$W)
l  <- c(rep("SE",nrow(SE)), rep("TE",nrow(TE)) )
df <- data.frame(L=l, H=H, R=R, A=A, W=w)

#Fig.1E (Axis labels were added using PowerPoint after this plot)
g <- ggplot(df,aes(x = W/100000, fill = L))+
  geom_histogram(bins = 100)+
  theme_bw(base_size = 20)+
  theme(panel.grid = element_blank(), strip.background = element_blank(),strip.text = element_blank())+
  facet_wrap(~L,ncol=1,scales = "free_y")+
  scale_fill_manual(values = c("firebrick2","royalblue2") ) +
  labs(x="Length", y="Count")+
  theme(axis.title.y = element_blank())+
  theme(axis.title.x = element_blank())+
  theme(axis.text.x = element_text(size = 25))+
  theme(axis.text.y = element_text(size = 25))+
  theme(panel.background = element_rect(fill = "white", colour = "black", size = 1.5))+
  guides(fill = F)

plot(g)

#Fig.1F H3K27Ac vs Length (Axis labels were added using PowerPoint after this plot)
g <- ggplot(df, aes(x=W/100000, y=H/10000, fill=L, alpha = 0.5))+
  theme_classic(base_size = 20)+
  theme(legend.position = "none", strip.background = element_blank(),
        strip.text = element_blank())+
  scale_fill_manual(values = c("firebrick2","royalblue2") ) +
  labs(x="Length", y="Norm. H3K27Ac signal")+
  theme(axis.title.y = element_blank())+
  theme(axis.title.x = element_blank())+
  theme(axis.text.x = element_text(size = 25))+
  theme(axis.text.y = element_text(size = 25))+
  geom_point(size=2,shape=21)

plot(g)

cor(w,H) #Peason correlation

#Fig.1G RelA vs Length (Axis labels were added using PowerPoint after this plot)
g <- ggplot(df, aes(x=W/100000, y=R/10000, fill=L, alpha = 0.5))+
  theme_classic(base_size = 20)+
  theme(legend.position = "none", strip.background = element_blank(),
        strip.text = element_blank())+
  scale_fill_manual(values = c("firebrick2","royalblue2") ) +
  labs(x="Length", y="Norm. RelA signal")+
  theme(axis.title.y = element_blank())+
  theme(axis.title.x = element_blank())+
  theme(axis.text.x = element_text(size = 25))+
  theme(axis.text.y = element_text(size = 25))+
  geom_point(size=2,shape=21)

plot(g)

cor(w,R) #Peason correlation

#Fig.1H H3K27Ac vs RelA (Axis labels were added using PowerPoint after this plot)
g <- ggplot(df, aes(x=R/10000, y=H/10000, fill=L, alpha = 0.5))+
  theme_classic(base_size = 20)+
  theme(legend.position = "none", strip.background = element_blank(),
        strip.text = element_blank())+
  scale_fill_manual(values = c("firebrick2","royalblue2") ) +
  labs(x="Norm. RelA signal", y="Norm. H3K27Ac signal")+
  theme(axis.title.y = element_blank())+
  theme(axis.title.x = element_blank())+
  theme(axis.text.x = element_text(size = 25))+
  theme(axis.text.y = element_text(size = 25))+
  geom_point(size=2,shape=21)

plot(g)

cor(H,R) #Peason correlation

#SFig.3A ATAC vs Length (Axis labels were added using PowerPoint after this plot)
g <- ggplot(df, aes(x=W/100000, y=A/10000, fill=L, alpha = 0.5))+
  theme_classic(base_size = 20)+
  theme(legend.position = "none", strip.background = element_blank(),
        strip.text = element_blank())+
  scale_fill_manual(values = c("firebrick2","royalblue2") ) +
  labs(x="Length", y="Norm. ATAC signal")+
  theme(axis.title.y = element_blank())+
  theme(axis.title.x = element_blank())+
  theme(axis.text.x = element_text(size = 25))+
  theme(axis.text.y = element_text(size = 25))+
  geom_point(size=2,shape=21)

plot(g)

cor(w,A) #Peason correlation