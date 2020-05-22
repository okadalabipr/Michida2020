library(tidyverse)
library(optimx)

#These two files are obtained using Filtering_genes.R. Please specify their path depending on your environment.
SE <- read.csv("filtered_Gained_SE.csv",stringsAsFactors = F)
TE <- read.csv("filtered_Gained_TE.csv",stringsAsFactors = F)

#Fig.3B-----------------------------------------------------------------------------------------------------
stats_se <- quantile(SE$peak_num,c(1/10,1/4,1/2,3/4,9/10)) %>% as.numeric()
stats_te <- quantile(TE$peak_num,c(1/10,1/4,1/2,3/4,9/10)) %>% as.numeric()

df <- data.frame(rbind(stats_se,stats_te))
df$L <- c("Gained SE","Gained TE")

df$mean <- c(mean(SE$peak_num),mean(TE$peak_num))

p <- ggplot(df, aes(x=L, lower=X2, upper=X4, middle=X3, ymin=X1, ymax=X5, fill=L))+
  theme_classic(base_size = 20)+
  theme(legend.position = "none")+
  scale_fill_manual(values = c("firebrick3","royalblue") ) +
  scale_y_continuous(limits=c(0,8), breaks=seq(0,8,2))+
  labs(x="", y="Number of RelA ChIP peaks")+
  geom_errorbar(width=0.5)+
  theme(axis.title.y = element_text(size = 20))+
  theme(axis.text.x = element_text(size = 20))+
  theme(axis.text.y = element_text(size = 20))+
  geom_boxplot(stat="identity")

p <- p + geom_point(data = df,aes(x = L,y = mean),fill = "white", size = 3, shape = 23)

plot(p)

t.test(SE$peak_num, TE$peak_num, var.equal = F)$p.value #t-test

mean(SE$peak_num) #Mean value of Gained SE
mean(TE$peak_num) #Mean value of Gained TE

#SFig.4A-------------------------------------------------------------------------------------------
stats_se <- quantile(SE$motif_num,c(1/10,1/4,1/2,3/4,9/10)) %>% as.numeric()
stats_te <- quantile(TE$motif_num,c(1/10,1/4,1/2,3/4,9/10)) %>% as.numeric()

df <- data.frame(rbind(stats_se,stats_te))
df$L <- c("Gained SE","Gained TE")

df$mean <- c(mean(SE$motif_num),mean(TE$motif_num))

p <- ggplot(df, aes(x=L, lower=X2, upper=X4, middle=X3, ymin=X1, ymax=X5, fill=L))+
  theme_classic(base_size = 20)+
  theme(legend.position = "none")+
  scale_fill_manual(values = c("firebrick3","royalblue") ) +
  scale_y_continuous(limits=c(0,12), breaks=seq(0,12,3))+
  labs(x="", y="Number of RelA ChIP peaks")+
  geom_errorbar(width=0.5)+
  theme(axis.title.y = element_text(size = 20))+
  theme(axis.text.x = element_text(size = 20))+
  theme(axis.text.y = element_text(size = 20))+
  geom_boxplot(stat="identity")

p <- p + geom_point(data = df,aes(x = L,y = mean),fill = "white", size = 3, shape = 23)

plot(p)

t.test(SE$motif_num, TE$motif_num, var.equal = F)$p.value #t-test

mean(SE$motif_num) #Mean value of Gained SE
mean(TE$motif_num) #Mean value of Gained TE

#Fig.3D------------------------------------------------------------------------------------------------

#model
exp_x   = c(0, 23.5, 30.9, 44, 100) # NF-kB activity from Shinohara et al., 2014, Science.
sim_x  = seq(0, max(exp_x), length.out = 1000)

model  = function(parS, sim_x) (sim_x^parS[2] / (sim_x^parS[2] + parS[3]*(parS[1]^parS[2]+sim_x^parS[2])))*parS[4]

parStart <- c(50, 1, 1, 1)


#Gained SE
obs_SE <- apply(SE[,25:28]-1,2,mean) #25:28 specifies FC_0.01 ~ FC_10 it can change depending on your process.
obs_SE <- append(0,obs_SE)
SE_sd <- apply(SE[,25:28]-1,2,sd)
SE_sd <- append(0,SE_sd)

resid_SE = function(p) sum((obs_SE-model(p,exp_x))^2) #Objective function of Gained SE

opt_SE = optimx(par = parStart, fn = resid_SE, control=list(all.methods = T),
                lower = rep(1e-5,4), upper = c(Inf,Inf,Inf,10*obs_SE[5]))

par_SE <- opt_SE[which.min(opt_SE$value),1:4]

df1_SE = data.frame( x=exp_x, v=obs_SE, sd = SE_sd)
df2_SE = data.frame( x=sim_x, v=model(as.numeric(par_SE),sim_x))

#Gained TE
obs_TE <- apply(TE[,25:28]-1,2,mean)
obs_TE <- append(0,obs_TE)
TE_sd <- apply(TE[,25:28]-1,2,sd)
TE_sd <- append(0,TE_sd)

resid_TE = function(p) sum( (obs_TE-model(p,exp_x))^2 ) #Objective function of Gained TE

opt_TE = optimx(par = parStart, fn = resid_TE,  control=list(all.methods = T),
                lower = rep(1e-5,4), upper = c(Inf,Inf,Inf,10*obs_TE[5]))

par_TE <- opt_TE[which.min(opt_SE$value),1:4]


df1_TE = data.frame( x=exp_x, v=obs_TE, sd = TE_sd)
df2_TE = data.frame( x=sim_x, v=model(as.numeric(par_TE),sim_x))

#Plot the graph
g = ggplot(df2_SE, aes(x = x, y = v))+
  geom_line(color = "firebrick3",size = 1.5)

g = g + 
  geom_point(data = df1_SE, aes(x = x, y = v), color = "firebrick3", size = 6)+
  #You can plot SFig.4C canceling the comment out below.
  #geom_errorbar(data = df1_SE, aes(ymin = v - sd, ymax = v + sd),width = 3,color = "firebrick3",alpha=0.5,size=1)+
  theme_bw()

g = g + geom_line(data = df2_TE, aes(x = x, y = v), color = "royalblue", size = 1.5)

g = g + 
  geom_point(data = df1_TE, aes(x = x, y = v),color = "royalblue", size = 6)+
  #You can plot SFig.4C canceling the comment out below.
  #geom_errorbar(data = df1_TE, aes(ymin = v - sd, ymax = v + sd),width = 3,color = "royalblue",alpha=0.5,size=1)+
  theme_bw()+
  xlab("Nuclear NFkB level (A.U.)")+
  ylab("Fold change in mRNA level - 1")+
  theme(axis.text.x = element_text(size = 30),axis.text.y = element_text(size=30))+
  theme(axis.title.x = element_text(size = 30),axis.title.y = element_text(size = 30))+
  theme(panel.grid = element_blank(), strip.background = element_blank(),strip.text = element_blank())+
  theme(plot.title = element_text(size = 25))+
  theme(panel.background = element_rect(fill = "white", colour = "black", size = 1.5))

plot(g)

par_SE
par_TE