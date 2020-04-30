setwd("dir_scATAC/peakcall")
library(tidyverse)

get_exact_motif <- function(E, dir_name){
  
  motif_tbl <- read.table(sprintf("%s_PU1_in_ATAC_in_%s.txt", dir_name, E),header=T, stringsAsFactors=F, sep="\t")
  
  peak_tbl <- read.table(sprintf("%s_ATAC_in_%s.hb", dir_name, E),header=F, stringsAsFactors=F, sep="\t")[,1:4]
  colnames(peak_tbl) <- c("PositionID","chr","start","end")
  
  motif_in_peak <- dplyr::inner_join(peak_tbl,motif_tbl,by="PositionID")
  
  motif_in_peak$center = (motif_in_peak$start + motif_in_peak$end)/2
  
  motif_in_peak$motif_position = motif_in_peak$center + motif_in_peak$Offset
  
  motif_in_peak$Sequence = as.character(motif_in_peak$Sequence)
  
  motif_in_peak$nchar = nchar(motif_in_peak$Sequence)
  
  motif_in_peak_positive <- dplyr::filter(motif_in_peak,Strand == "+")
  
  motif_in_peak_negative <- dplyr::filter(motif_in_peak,Strand == "-")
  
  strand_positive <- function(x){
    x$motif_s <- x$motif_position
    x$motif_e <- x$motif_position + (x$nchar - 1)
    return(x)
  }
  
  strand_negative <- function(x){
    x$motif_s <- x$motif_position - (x$nchar - 1)
    x$motif_e <- x$motif_position
    return(x)
  }
  
  
  motif_in_peak_positive <- strand_positive(motif_in_peak_positive)
  
  motif_in_peak_negative <- strand_negative(motif_in_peak_negative)
  
  motif_in_peak <- rbind(motif_in_peak_positive,motif_in_peak_negative)
  
  
  motif_in_peak$motif_s = floor(motif_in_peak$motif_s) + 1
  motif_in_peak$motif_e = floor(motif_in_peak$motif_e) + 1
  
  output <- data.frame(chr = motif_in_peak$chr, start = motif_in_peak$motif_s, end = motif_in_peak$motif_e)
  
  write.table(output, sprintf("%s_PU1_in_ATAC_in_%s.bed", dir_name, E), 
              row.names = F, col.names = F, sep = "\t", quote = F)
  
}

dirs <- list.dirs()
dirs <- dirs[2:length(dirs)] #The first item is current directory ("./") and skip it.

for(i in 1:length(dirs)){
  setwd(dirs[i])
  dir_name <- substr(dirs[i], 3, nchar(dirs[i]))
  get_exact_motif("SE", dir_name)
  get_exact_motif("TE", dir_name)
  setwd("../")
}

