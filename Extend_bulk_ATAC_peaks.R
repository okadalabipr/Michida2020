tbl = read.table("path to your peak summits file of bulk ATAC-seq obtained using peak_call_for_bulk.sh")[,1:3]

colnames(tbl) <- c("chr","start","end")

extended_tbl = data.frame(chr = tbl$chr, start = tbl$start-249, end = tbl$end+250)

#This file is required for downstream analysis.
write.table(extended_tbl,"extended.bed",quote = F, row.names = F, sep = "\t", col.names = F)