library(ggplot2)
library(tidyverse)

TXT <- scan("path to your ${TXT} output by SFig1D.sh",skip = 1, nlines = 4, what = "" , sep = "\t")

df <- data.frame(bin = c(seq(1.0,20000.0,1.0),seq(1.0,20000.0,1.0)),
                 value = c(TXT[(match("SE",TXT)+1):(match("TE",TXT)-2)] %>% as.numeric(),
                                TXT[(match("TE",TXT)+1):length(TXT)] %>% as.numeric()),
                      type = c(rep("SE",20000),rep("TE",20000)))

g <- ggplot(df,aes(x = bin, y = value, color = type))+
  geom_line()+
  theme_classic(base_size = 20)+
  theme(legend.position = "none", strip.background = element_blank(),
        strip.text = element_blank())+
  scale_color_manual(values = c("firebrick2","royalblue2") ) +
  labs(x="Genomic region", y="Norm. coverage")+
  theme(axis.title.y = element_text(size = 30))+
  theme(axis.title.x = element_text(size = 30))+
  theme(axis.text.x = element_text(size = 25))+
  theme(axis.text.y = element_text(size = 25))+
  scale_x_continuous(breaks = c(0,5000,15000,20000),labels=c("-5kb", "Start", "End", "+5kb"))


plot(g)