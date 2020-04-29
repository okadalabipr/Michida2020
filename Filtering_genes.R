library(ggplot2)
library(tidyverse)

#These two files were obtained using Make_data_table.R. Please specify their path depending on your environment.
SE_sg <- read.csv("cleaned_Gained_SE.csv",stringsAsFactors = F)

TE_sg <- read.csv("cleaned_Gained_TE.csv", stringsAsFactors = F)

df <- rbind(SE_sg,TE_sg)

df <- df %>% distinct(gene, .keep_all = TRUE) %>% select(gene,TPM_0,TPM_0.01,TPM_0.1,TPM_1,TPM_10,
                                                         FC_0.01,FC_0.1,FC_1,FC_10)

#See histgram in each step of filtering
get_histgram <- function(df,x,BINS,line,XLAB){
  g <- ggplot(df,aes(x = x))+
    geom_histogram(bins = BINS)+
    theme_bw()+
    xlab(XLAB)+
    geom_vline(xintercept = line, linetype="dashed")+
    theme(axis.text.x = element_text(size=15),axis.text.y = element_text(size=15))+
    theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20))
  plot(g)
}

get_histgram(df,log2(df$TPM_0),100,0,"log2TPM")

df_TPM_OK <- filter(df, log2(TPM_0) >= 0)
df_TPM_OUT <- filter(df, log2(TPM_0) < 0)

get_histgram(df_TPM_OK,df_TPM_OK$FC_10,100,0,"FC_10")

df_TPM_OK <- filter(df_TPM_OK, FC_10 < 350) #In the case of this paper, the extream value was excluded.

get_histgram(df_TPM_OK,log2(df_TPM_OK$FC_10),100,0,"log2FC_10")

df_TPM_FC_tmp <- filter(df_TPM_OK, FC_10 > 1)

get_histgram(df_TPM_FC_tmp,log2(df_TPM_FC_tmp$FC_10),100,0.2,"log2FC_10")

df_TPM_FC_OK <- filter(df_TPM_FC_tmp, log2(FC_10) >= 0.2)

df_TPM_FC_OUT <- filter(df_TPM_OK, log2(FC_10) < 0.2)

filtered_df <- df_TPM_FC_OK

#Enhancers whose target genes passed this filtering will move on to downstream analysis.
v <- c()
for(i in 1:nrow(filtered_df)){
  v <- append(v,which(filtered_df$gene[i] == SE_sg$gene))
}

SE_sg_filtered <- SE_sg[v,]

v <- c()
for(i in 1:nrow(filtered_df)){
  v <- append(v,which(filtered_df$gene[i] == TE_sg$gene))
}

TE_sg_filtered <- TE_sg[v,]

#These two files are required for the downstream analysis.
write.csv(SE_sg_filtered,"filtered_Gained_SE.csv",row.names = F)
write.csv(TE_sg_filtered,"filtered_Gained_TE.csv",row.names = F)
