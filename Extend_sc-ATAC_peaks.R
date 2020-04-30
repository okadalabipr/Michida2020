setwd("dir_scATAC/peakcall")
library(tidyverse)

dirs <- list.dirs()

for(i in 1:(length(dirs)-1)){
  setwd(dirs[i+1])
  dir_name <- substr(dirs[i+1], 3, 8)
  tsv <- read.table(sprintf("%s_peaks.tsv",dir_name),header = T)
  delete <- c()
  for(j in 1:(nrow(tsv)-1)){
    if(tsv$chr[j+1] == tsv$chr[j] & tsv$start[j+1] == tsv$start[j] & tsv$end[j+1] == tsv$end[j]){
      delete <- append(delete,j+1)
    }
  }
  tsv <- tsv[-delete,]
  
  center <- floor((tsv$start + tsv$end)/2)
  fixed_start <- center - 250
  fixed_end <- center + 250
  
  bed <- data.frame(chr = tsv$chr, start = fixed_start, end = fixed_end)
  write.table(bed, sprintf("%s_fixed_peaks.bed",dir_name), sep = "\t", col.names = F, quote = F, row.names = F)
  setwd("../")
}
