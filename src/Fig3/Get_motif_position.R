library(tidyverse)
library(ggplot2)

motif_tbl <- read.table("path to your findMotigsGenome.pl output file of Get_NFkB_motifs_in_enhancers.sh",
                        header=T, stringsAsFactors=F, sep="\t")

peak_tbl <- read.table("path to HOMER BED file used as the motif search region in findMotifsGenome.pl 
                       in Get_NFkB_motifs_in_enhancers.sh")[,1:4]

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

write.table(motif_in_peak[,c(2,13,14)],"path to output file (BED format)",quote = F, row.names = F, sep = "\t", col.names = F)
