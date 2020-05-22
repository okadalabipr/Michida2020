library(tidyverse)
library(ggplot2)
library(optimx)

motif_in_peak <- read.table("path to RelA-PU.1 overlaps in anti-IgM 060 min condition obtained using Get_RelA_PU1_overlap.sh")
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
    En_chr$overlap_num = v
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
    En_chr$overlap_num = v
    TE_sg <- rbind(TE_sg,En_chr)
  }
}


#Fig5A-----------------------------------------------------------------------------------------------------
stats_se <- quantile(SE_sg$overlap_num,c(1/10,1/4,1/2,3/4,9/10)) %>% as.numeric()
stats_te <- quantile(TE_sg$overlap_num,c(1/10,1/4,1/2,3/4,9/10)) %>% as.numeric()


df <- data.frame(rbind(stats_se,stats_te))
df$L <- c("Gained SE","Gained TE")

df$mean <- c(mean(SE_sg$overlap_num),mean(TE_sg$overlap_num))


p <- ggplot(df, aes(x=L, lower=X2, upper=X4, middle=X3, ymin=X1, ymax=X5, fill=L))+
  theme_classic(base_size = 20)+
  theme(legend.position = "none")+
  scale_fill_manual(values = c("firebrick3","royalblue") ) +
  #scale_y_continuous(limits=c(0,15), breaks=seq(0,15,5))+
  labs(x="", y="Number of RelA ChIP peak - PU.1 motif overlaps")+
  geom_errorbar(width=0.5)+
  theme(axis.text.x = element_text(size = 20))+
  theme(axis.text.y = element_text(size = 20))+
  geom_boxplot(stat="identity")

p <- p + geom_point(data = df,aes(x = L,y = mean),fill = "white", size = 3, shape = 23)

plot(p)

t.test(SE_sg$overlap_num,TE_sg$overlap_num,var.equal = F)$p.value #t-test

mean(SE_sg$overlap_num) #Mean value of gained SEs
mean(TE_sg$overlap_num) #Mean Value of gained TEs

#non-zero percentage
100*(which(SE_sg$overlap_num > 0) %>% length())/nrow(SE_sg)
100*(which(TE_sg$overlap_num > 0) %>% length())/nrow(TE_sg)

#Mathematical model fitting to individual genes for Fig.5D-E-----------------------------------------------------------
#model
exp_x   = c(0, 23.5, 30.9, 44, 100) # NF-kB activity from Shinohara et al., 2014, Science.
sim_x  = seq(0, max(exp_x), length.out = 1000)

model  = function(parS, sim_x) (sim_x^parS[2] / (sim_x^parS[2] + parS[3]*(parS[1]^parS[2]+sim_x^parS[2])))*parS[4]

parStart <- c(50, 1, 1, 1)

df <- rbind(SE_sg, TE_sg)

skipped <- c()
Nopt <- c()
for(i in 1:nrow(df)){
  obs = as.numeric(df[i,25:28]-1) #25:28 specifies FC_0.01 ~ FC_10.
  obs = append(0,obs)
  resid = function(p) sum((obs-model(p,exp_x))^2)
  opt = try(optimx(par = parStart, fn = resid,control=list(all.methods = T), 
                   lower = rep(1e-5,4), upper = c(Inf,12,Inf,10*obs[5])))
  if (class(opt) == "try-error"){
    skipped <- append(skipped,i)
    Nopt <- append(Nopt, -1)
    next
  }
  
  par <- opt[which.min(opt$value),1:4] %>% as.numeric()
  Nopt <- append(Nopt,par[2])
}

df$Nopt <- Nopt

SE_sg <- df[1:nrow(SE_sg),]
TE_sg <- df[(nrow(SE_sg)+1):nrow(df),]

#Fig.5B-----------------------------------------------------------------------------------------------------------
SE_with_overlap <- filter(SE_sg, overlap_num > 0)
SE_without_overlap <- filter(SE_sg, overlap_num == 0)

TE_with_overlap <- filter(TE_sg,  overlap_num > 0)
TE_without_overlap <- filter(TE_sg,  overlap_num == 0)

SE_with_overlap$length <- SE_with_overlap$end - SE_with_overlap$start
SE_without_overlap$length <- SE_without_overlap$end - SE_without_overlap$start
TE_with_overlap$length <- TE_with_overlap$end - TE_with_overlap$start
TE_without_overlap$length <- TE_without_overlap$end - TE_without_overlap$start

stats_se_with <- quantile(SE_with_overlap$length,c(1/10,1/4,1/2,3/4,9/10)) %>% as.numeric()
stats_se_without <- quantile(SE_without_overlap$length,c(1/10,1/4,1/2,3/4,9/10)) %>% as.numeric()

stats_te_with <- quantile(TE_with_overlap$length,c(1/10,1/4,1/2,3/4,9/10)) %>% as.numeric()
stats_te_without <- quantile(TE_without_overlap$length,c(1/10,1/4,1/2,3/4,9/10)) %>% as.numeric()

df_SE <- data.frame(rbind(stats_se_with,stats_se_without))
df_SE$L <- c("SE with overlap","SE without overlap")

df_TE <- data.frame(rbind(stats_te_with,stats_te_without))
df_TE$L <- c("TE with overlap","TE without overlap")

df_SE$mean <- c(mean(SE_with_overlap$length),mean(SE_without_overlap$length))
df_TE$mean <- c(mean(TE_with_overlap$length),mean(TE_without_overlap$length))

#Fig.5B (Gained SE)
p_SE <- ggplot(df_SE, aes(x=L, lower=X2, upper=X4, middle=X3, ymin=X1, ymax=X5, fill=L))+
  theme_classic(base_size = 20)+
  theme(legend.position = "none")+
  scale_fill_manual(values = c("salmon", rgb(1, 200/255, 190/255))) +
  #scale_y_continuous(limits=c(0,6), breaks=seq(0,6,2))+
  labs(x="", y="")+
  geom_errorbar(width=0.5)+
  theme(axis.text.x = element_text(size = 20))+
  theme(axis.text.y = element_text(size = 20))+
  geom_boxplot(stat="identity")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

p_SE <- p_SE + geom_point(data = df_SE,aes(x = L,y = mean),fill = "white", size = 3, shape = 23)

plot(p_SE)

mean(SE_with_overlap$length);mean(SE_without_overlap$length) #Mean value
t.test(SE_with_overlap$length, SE_without_overlap$length,var.equal = F)$p.value #t-test

#Fig.5B (Gained TE)
p_TE <- ggplot(df_TE, aes(x=L, lower=X2, upper=X4, middle=X3, ymin=X1, ymax=X5, fill=L))+
  theme_classic(base_size = 20)+
  theme(legend.position = "none")+
  scale_fill_manual(values = c("cyan3", rgb(120/255, 240/255, 242/255))) +
  #scale_y_continuous(limits=c(0,6), breaks=TEq(0,6,2))+
  labs(x="", y="")+
  geom_errorbar(width=0.5)+
  theme(axis.text.x = element_text(size = 20))+
  theme(axis.text.y = element_text(size = 20))+
  geom_boxplot(stat="identity")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

p_TE <- p_TE + geom_point(data = df_TE,aes(x = L,y = mean),fill = "white", size = 3, shape = 23)

plot(p_TE)
mean(TE_with_overlap$length);mean(TE_without_overlap$length) #Mean value
t.test(TE_with_overlap$length, TE_without_overlap$length,var.equal = F)$p.value #t-test

#Fig.5C--------------------------------------------------------------------------------------------------------
stats_se_with <- quantile(SE_with_overlap$log2FC_10,c(1/10,1/4,1/2,3/4,9/10)) %>% as.numeric()
stats_se_without <- quantile(SE_without_overlap$log2FC_10,c(1/10,1/4,1/2,3/4,9/10)) %>% as.numeric()

stats_te_with <- quantile(TE_with_overlap$log2FC_10,c(1/10,1/4,1/2,3/4,9/10)) %>% as.numeric()
stats_te_without <- quantile(TE_without_overlap$log2FC_10,c(1/10,1/4,1/2,3/4,9/10)) %>% as.numeric()

df <- data.frame(rbind(stats_se_with,stats_se_without,stats_te_with,stats_te_without))
df$L <- c("SE with overlap","SE without overlap","TE with overlap","TE without overlap")

df$mean <- c(mean(SE_with_overlap$log2FC_10),mean(SE_without_overlap$log2FC_10),
             mean(TE_with_overlap$log2FC_10),mean(TE_without_overlap$log2FC_10))

p <- ggplot(df, aes(x=L, lower=X2, upper=X4, middle=X3, ymin=X1, ymax=X5, fill=L))+
  theme_classic(base_size = 20)+
  theme(legend.position = "none")+
  scale_fill_manual(values = c("salmon", rgb(1, 200/255, 190/255), 
                               "cyan3", rgb(120/255, 240/255, 242/255))) +
  #scale_y_continuous(limits=c(0,10), breaks=seq(0,10,2))+
  labs(x="", y="Fold change in mRNA level (log2)")+
  geom_errorbar(width=0.5)+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 20))+
  theme(axis.text.y = element_text(size = 20))+
  geom_boxplot(stat="identity")

p <- p + geom_point(data = df,aes(x = L,y = mean),fill = "white", size = 3, shape = 23)

plot(p)

df$mean #Mean values

#Get p.adjust values
p1 <- t.test(SE_with_overlap$log2FC_10, SE_without_overlap$log2FC_10,var.equal = F)$p.value
p2 <- t.test(SE_without_overlap$log2FC_10, TE_without_overlap$log2FC_10,var.equal = F)$p.value
p3 <- t.test(SE_with_overlap$log2FC_10, TE_with_overlap$log2FC_10,var.equal = F)$p.value
p4 <- t.test(SE_with_overlap$log2FC_10, TE_without_overlap$log2FC_10,var.equal = F)$p.value
p5 <- t.test(SE_without_overlap$log2FC_10, TE_with_overlap$log2FC_10,var.equal = F)$p.value
p6 <- t.test(TE_with_overlap$log2FC_10, TE_without_overlap$log2FC_10,var.equal = F)$p.value
q <- p.adjust(c(p1,p2,p3,p4,p5,p6), method = "BH");q

#Fig.5D--------------------------------------------------------------------------------------------------
stats_se_with <- quantile(SE_with_overlap$Nopt,c(1/10,1/4,1/2,3/4,9/10)) %>% as.numeric()
stats_se_without <- quantile(SE_without_overlap$Nopt,c(1/10,1/4,1/2,3/4,9/10)) %>% as.numeric()

stats_te_with <- quantile(TE_with_overlap$Nopt,c(1/10,1/4,1/2,3/4,9/10)) %>% as.numeric()
stats_te_without <- quantile(TE_without_overlap$Nopt,c(1/10,1/4,1/2,3/4,9/10)) %>% as.numeric()

df <- data.frame(rbind(stats_se_with,stats_se_without,stats_te_with,stats_te_without))
df$L <- c("SE with overlap","SE without overlap","TE with overlap","TE without overlap")

df$mean <- c(mean(SE_with_overlap$Nopt),mean(SE_without_overlap$Nopt),
             mean(TE_with_overlap$Nopt),mean(TE_without_overlap$Nopt))

p <- ggplot(df, aes(x=L, lower=X2, upper=X4, middle=X3, ymin=X1, ymax=X5, fill=L))+
  theme_classic(base_size = 20)+
  theme(legend.position = "none")+
  scale_fill_manual(values = c("salmon", rgb(1, 200/255, 190/255), 
                               "cyan3", rgb(120/255, 240/255, 242/255))) +
  scale_y_continuous(limits=c(0,12), breaks=seq(0,12,2))+
  labs(x="", y="Nopt")+
  geom_errorbar(width=0.5)+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 20))+
  theme(axis.text.y = element_text(size = 20))+
  geom_boxplot(stat="identity")

p <- p + geom_point(data = df,aes(x = L,y = mean),fill = "white", size = 3, shape = 23)

plot(p)

df$mean #Mean values

#Get p.adjust values
p1 <- t.test(SE_with_overlap$Nopt, SE_without_overlap$Nopt,var.equal = F)$p.value
p2 <- t.test(SE_without_overlap$Nopt, TE_without_overlap$Nopt,var.equal = F)$p.value
p3 <- t.test(SE_with_overlap$Nopt, TE_with_overlap$Nopt,var.equal = F)$p.value
p4 <- t.test(SE_with_overlap$Nopt, TE_without_overlap$Nopt,var.equal = F)$p.value
p5 <- t.test(SE_without_overlap$Nopt, TE_with_overlap$Nopt,var.equal = F)$p.value
p6 <- t.test(TE_with_overlap$Nopt, TE_without_overlap$Nopt,var.equal = F)$p.value
q <- p.adjust(c(p1,p2,p3,p4,p5,p6), method = "BH");q

#Fig.5E---------------------------------------------------------------------------------------
SE_with_overlap <- filter(SE_sg, overlap_num > 0)
SE_with_overlap <- filter(SE_with_overlap, Nopt != 12)
SE_with_overlap <- filter(SE_with_overlap, Nopt > 0.00001)
SE_without_overlap <- filter(SE_sg, overlap_num == 0)
SE_without_overlap <- filter(SE_without_overlap, Nopt != 12)
SE_without_overlap <- filter(SE_without_overlap, Nopt > 0.00001)

TE_with_overlap <- filter(TE_sg,  overlap_num > 0)
TE_with_overlap <- filter(TE_with_overlap, Nopt != 12)
TE_with_overlap <- filter(TE_with_overlap, Nopt > 0.00001)
TE_without_overlap <- filter(TE_sg,  overlap_num == 0)
TE_without_overlap <- filter(TE_without_overlap, Nopt != 12)
TE_without_overlap <- filter(TE_without_overlap, Nopt > 0.00001)

stats_se_with <- quantile(SE_with_overlap$Nopt,c(1/10,1/4,1/2,3/4,9/10)) %>% as.numeric()
stats_se_without <- quantile(SE_without_overlap$Nopt,c(1/10,1/4,1/2,3/4,9/10)) %>% as.numeric()

stats_te_with <- quantile(TE_with_overlap$Nopt,c(1/10,1/4,1/2,3/4,9/10)) %>% as.numeric()
stats_te_without <- quantile(TE_without_overlap$Nopt,c(1/10,1/4,1/2,3/4,9/10)) %>% as.numeric()

df <- data.frame(rbind(stats_se_with,stats_se_without,stats_te_with,stats_te_without))
df$L <- c("SE with overlap","SE without overlap","TE with overlap","TE without overlap")

df$mean <- c(mean(SE_with_overlap$Nopt),mean(SE_without_overlap$Nopt),
             mean(TE_with_overlap$Nopt),mean(TE_without_overlap$Nopt))

p <- ggplot(df, aes(x=L, lower=X2, upper=X4, middle=X3, ymin=X1, ymax=X5, fill=L))+
  theme_classic(base_size = 20)+
  theme(legend.position = "none")+
  scale_fill_manual(values = c("salmon", rgb(1, 200/255, 190/255), 
                               "cyan3", rgb(120/255, 240/255, 242/255))) +
  scale_y_continuous(limits=c(0,10), breaks=seq(0,10,2))+
  labs(x="", y="Nopt")+
  geom_errorbar(width=0.5)+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 20))+
  theme(axis.text.y = element_text(size = 20))+
  geom_boxplot(stat="identity")

p <- p + geom_point(data = df,aes(x = L,y = mean),fill = "white", size = 3, shape = 23)

plot(p)

df$mean #Mean values

#Get p.adjust values
p1 <- t.test(SE_with_overlap$Nopt, SE_without_overlap$Nopt,var.equal = F)$p.value
p2 <- t.test(SE_without_overlap$Nopt, TE_without_overlap$Nopt,var.equal = F)$p.value
p3 <- t.test(SE_with_overlap$Nopt, TE_with_overlap$Nopt,var.equal = F)$p.value
p4 <- t.test(SE_with_overlap$Nopt, TE_without_overlap$Nopt,var.equal = F)$p.value
p5 <- t.test(SE_without_overlap$Nopt, TE_with_overlap$Nopt,var.equal = F)$p.value
p6 <- t.test(TE_with_overlap$Nopt, TE_without_overlap$Nopt,var.equal = F)$p.value
q <- p.adjust(c(p1,p2,p3,p4,p5,p6), method = "BH");q

#SFig.6B-------------------------------------------------------------------------------------------------
SE_with_overlap <- filter(SE_sg, overlap_num > 0)
SE_without_overlap <- filter(SE_sg, overlap_num == 0)

TE_with_overlap <- filter(TE_sg,  overlap_num > 0)
TE_without_overlap <- filter(TE_sg,  overlap_num == 0)

get_param <- function(data){
  x   = c(0, 23.5, 30.9, 44, 100) 
  xx  = seq(0, 100, length.out = 1000)
  parStart <- c(50, 1, 1, 1)
  obs <- apply(data[,25:28]-1,2,mean) #25:28 specifies FC_0.01 ~ FC_10
  obs <- append(0,obs)
  sd <- apply(data[,25:28]-1,2,sd)
  sd <- append(0,sd)
  resid = function(p) sum((obs-model(p,x))^2)
  opt = optimx(par = parStart, fn = resid, control=list(all.methods = T),
               lower = rep(1e-5,4), upper = c(Inf,12,Inf,10*obs[5]))
  par <- opt[which.min(opt$value),1:4]
  return(list(par, obs, sd))
}

SE_with <- get_param(SE_with_overlap);SE_with[[1]]
SE_without <- get_param(SE_without_overlap);SE_without[[1]]

TE_with <- get_param(TE_with_overlap);TE_with[[1]]
TE_without <- get_param(TE_without_overlap);TE_without[[1]]

get_sim <- function(fit_res,iro){
  x   = c(0, 23.5, 30.9, 44, 100) # NF-kB activity from Shinohara et al., 2014, Science
  xx  = seq(0, 100, length.out = 1000)
  df1 = data.frame(x=x, v=fit_res[[2]], sd = fit_res[[3]])
  df2 = data.frame(x=xx, v=model(as.numeric(fit_res[[1]]),xx))
  
  g = ggplot(df2, aes(x = x, y = v))+
    geom_line(color = iro, size = 1.5)
  g = g + 
    geom_point(data = df1, aes(x = x, y = v), color = iro, size = 6)+
    #scale_y_continuous(limits = c(-0.5,2.0),breaks = seq(-0.5,2.0,0.5))+
    ggtitle(sprintf("Km = %s, N = %s",round(fit_res[[1]][1],digits = 1),round(fit_res[[1]][2],digits = 1)))+
    #geom_errorbar(data = df1_SE, aes(ymin = v - sd, ymax = v + sd),width = 3,color = iro ,alpha=0.5,size=1)+
    theme_bw()+
    theme(axis.text.x = element_text(size = 30),axis.text.y = element_text(size=30))+
    theme(axis.title.x = element_blank(),axis.title.y = element_blank())+
    theme(panel.grid = element_blank(), strip.background = element_blank(),strip.text = element_blank())+
    theme(panel.background = element_rect(fill = "white", colour = "black", size = 2.0))
  
  plot(g)
}

get_sim(SE_with, "salmon")
get_sim(SE_without, rgb(1, 200/255, 190/255))
get_sim(TE_with, "cyan3")
get_sim(TE_without, rgb(120/255, 240/255, 242/255))

