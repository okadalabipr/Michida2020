library(ggplot2)
library(tidyverse)

count <- read.table("path to your output file of Count_RelA_PU1_overlap_for_scATAC.sh", 
                    header = F, stringsAsFactors = F)$V1

SE <- count[2:(grep("TE",count)-1)] %>% as.numeric()
TE <- count[(grep("TE",count)+1):length(count)] %>% as.numeric()

SE_ratio <- SE/280 #280 is number of SEs in the catalog.
TE_ratio <- TE/8876 #8876 is number of TEs in the catalog.

N <- length(SE_ration) + length(TE_ratio)

df <- data.frame(ratio = c(SE_ratio, TE_ratio), label = c(rep("SE",N), rep("TE",N)))

#Plot the graph
g <- ggplot(df, aes(x = label, y = ratio, color = label, alpha = 0.5))+
  geom_jitter()+
  theme_classic(base_size = 20)+
  theme(legend.position = "none")+
  scale_color_manual(values =c("firebrick2","royalblue2"))+
  theme(axis.text.y = element_text(size = 20))+
  theme(axis.text.x = element_text(size = 20))+
  labs(x = "", y="overlaps / enhancer")

plot(g)

t.test(SE_ratio, TE_ratio,var.equal = F)$p.value #t-test
