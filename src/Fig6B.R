library(ggplot2)
library(gridExtra)

#Fig.6B left pannel------------------------------------------------------------------------------------------------
model <- function(parS, sim_x) (sim_x^parS[2] / (sim_x^parS[2] + parS[3]*(parS[1]^parS[2]+sim_x^parS[2])))*parS[4]

#The parameters in this paper. (Km, N, K1, K2)
parm_SE <- c(53.8399221785168,	4.22343199706818,	3.37923332573206,	7.97302213444357)
parm_TE <- c(262.284020468864,	1.28724329750798,	2.01829647337324,	10.4313326529265)

#NFkB <- rnorm(1000, mean = 50, sd = 10) #Make hypothetical single cell nuclear NFkB distribution

NFkB <- readRDS("hypo_NFkB.obj")#Read hypothetical single cell nuclear NFkB distribution used in this paper
                                #"hypo_NFkB.obj" is provided in our repository.

#Estimated mRNA distribution
df <- data.frame(value = c(model(parm_SE,NFkB),model(parm_TE,NFkB)),
                 type = c(rep("SE",length(NFkB)),rep("TE",length(NFkB))))

mRNA_dist <- ggplot(df,aes(x = value, fill = type))+
  geom_histogram(bins = 100)+
  facet_wrap(~type,ncol = 1,scales = "free_y")+
  theme_bw()+
  theme(panel.grid = element_blank(), strip.background = element_blank(),strip.text = element_blank())+
  labs(x="mRNA level (A.U.)", y="Count")+
  theme(axis.title.y = element_text(size = 25))+
  theme(axis.title.x = element_text(size = 25))+
  theme(axis.text.x = element_text(size = 20))+
  theme(axis.text.y = element_text(size = 20))+
  scale_fill_manual(values =c("firebrick3","royalblue"))+
  theme(plot.title = element_blank())+
  theme(panel.background = element_rect(fill = "white", colour = "black", size = 1.2))+
  guides(alpha = F, fill = F)


plot(mRNA_dist)

#Fig.6B right pannel-------------------------------------------------------------------------------------------
TPM <- read.table("cleaned_TPM.txt", header=T, stringsAsFactors=F,sep="\t") #This file is obtained using Make_data_table.R.

#These two files are obtained using Filtering_genes.R.
SE <- read.csv("filtered_Gained_SE.csv",stringsAsFactors = F)
TE <- read.csv("filtered_Gained_TE.csv",stringsAsFactors = F)

choose_dose <- c(1,grep("^IgM_01000_*",colnames(TPM))) #Choose anti-IgM 1 ug/ml condition to plot


SE_df <- NULL
for(i in 1:nrow(SE)){
  index <- which(TPM$Geneid == SE$gene[i])
  SE_df <- rbind(SE_df,TPM[index,choose_dose])
}

TE_df <- NULL
for(i in 1:nrow(TE)){
  index <- which(TPM$Geneid == TE$gene[i])
  TE_df <- rbind(TE_df,TPM[index,choose_dose])
}

SE_mean <- apply(SE_df[,2:ncol(SE_df)],2,mean)
TE_mean <- apply(TE_df[,2:ncol(TE_df)],2,mean)


df <- data.frame(v = c(SE_mean,TE_mean),type = c(rep("SE",length(SE_mean)),rep("TE",length(TE_mean))))

g <- ggplot(df,aes(x = v, fill = type))+
  geom_histogram(bins = 100)+
  facet_wrap(~type,ncol = 1, scales = "free_y")+
  theme_bw()+
  theme(panel.grid = element_blank(), strip.background = element_blank(),strip.text = element_blank())+
  scale_fill_manual(values = c("firebrick3","royalblue") ) +
  labs(x="CV", y="Count")+
  theme(axis.title.y = element_text(size = 25))+
  theme(axis.title.x = element_text(size = 25))+
  theme(axis.text.x = element_text(size = 20))+
  theme(axis.text.y = element_text(size = 20))+
  theme(panel.background = element_rect(fill = "white", colour = "black", size = 1.2))+
  guides(fill = F)

plot(g)

sd(SE_mean)/mean(SE_mean) #CV of Gained SEs
sd(TE_mean)/mean(TE_mean) #CV of Gained TEs

