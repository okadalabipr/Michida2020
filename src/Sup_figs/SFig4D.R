library(ggplot2)
library(tidyverse)
library(optimx)

#These two files are obtained using Filtering_genes.R. Please specify their path depending on your environment.
SE <- read.csv("filtered_Gained_SE.csv",stringsAsFactors = F)
TE <- read.csv("filtered_Gained_TE.csv",stringsAsFactors = F)

df <- rbind(
  SE[which(SE$gene == "Irf4")[1],],
  SE[which(SE$gene == "Smad3")[1],],
  SE[which(SE$gene == "Fam43a")[1],],
  TE[which(TE$gene == "Rela")[1],],
  TE[which(TE$gene == "Ran")[1],],
  TE[which(TE$gene == "Eif3d")[1],]
)

exp_x <- c(0, 23.5, 30.9, 44, 100) # NF-kB activity from Shinohara et al., 2014, Science
sim_x <- seq(0, max(exp_x), length.out = 1000)

par_df <- NULL

for(i in 1:nrow(df)){
  
  model <- function(parS, sim_x) (sim_x^df$peak_num[i] / (sim_x^df$peak_num[i] + parS[2]*(parS[1]^df$peak_num[i]+sim_x^df$peak_num[i])))*parS[3]
  parStart <- c(50, 1, 1)
  
  obs <- as.numeric(df[i,25:28]-1) #25:28 specifies FC_0.01 ~ FC_10
  obs <- append(0,obs)
  
  resid <- function(p) sum((obs-model(p,exp_x))^2)
  
  opt <- optimx(par = parStart, fn = resid, control=list(all.methods = T), 
                   lower = rep(1e-5,3), upper = c(Inf,Inf,10*obs[5]))
  
  par <- opt[which.min(opt$value),1:3] %>% as.numeric()
  
  par_df <- rbind(par_df,par)
}

par_df <- as.data.frame(par_df)
rownames(par_df) <- seq(1,nrow(par_df),1)
colnames(par_df) <- c("Km", "K1", "K2")

df <- cbind(df,par_df)


plot_SFig4D <- function(gene_name,iro){
  row_n <- which(df$gene == gene_name)
  plot_model <- function(parS, sim_x) (sim_x^df$peak_num[row_n] / (sim_x^df$peak_num[row_n] + parS[2]*(parS[1]^df$peak_num[row_n]+sim_x^df$peak_num[row_n])))*parS[3]
  
  obs = as.numeric(df[row_n,25:28]-1) #25:28 specifies FC_0.01 ~ FC_10
  obs = append(0,obs)
  
  par <- df[row_n,(ncol(df)-2):ncol(df)] %>% as.numeric()
  
  df1 = data.frame( x=exp_x, v=obs)
  df2 = data.frame( x=sim_x, v=plot_model(par,sim_x) )
  
  g = ggplot(df2, aes(x = x, y = v))+
    geom_line(color = iro, size = 1.5)
  
  g = g + geom_point(data = df1, aes(x = x, y = v),color = iro,size = 6)+
    theme_classic()+
    labs(x = "Nuclear NFkB level (A.U.)", y = "Fold change in mRNA level - 1")+
    theme(axis.text.x = element_text(size = 20),axis.text.y = element_text(size=20))+
    theme(axis.title.x = element_text(size = 20),axis.title.y = element_text(size=20))+
    theme(panel.grid = element_blank(), strip.background = element_blank(),strip.text = element_blank())+
    theme(panel.background = element_rect(fill = "white", colour = "black", size = 2.0))
  
  plot(g)
}

#Plot the graphs
plot_SFig4D("Irf4", "firebrick3")
plot_SFig4D("Smad3", "firebrick3")
plot_SFig4D("Fam43a", "firebrick3")
plot_SFig4D("Rela", "royalblue")
plot_SFig4D("Ran", "royalblue")
plot_SFig4D("Eif3d", "royalblue")