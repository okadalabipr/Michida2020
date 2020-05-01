library(clusterProfiler)
library(org.Mm.eg.db)

SE <- read.csv("path to ${SE_signal} made by Get_signal.sh",stringsAsFactors = F, sep = "\t")
TE <- read.csv("path to ${TE_signal} made by Get_signal.sh",stringsAsFactors = F, sep = "\t")

SE <- unique(SE$Entrez.ID)
TE <- unique(TE$Entrez.ID)

clusters = list(SE = SE, TE = TE)

ck <- compareCluster(geneCluster = clusters , fun = "enrichGO", OrgDb = "org.Mm.eg.db",
                     ont           = "BP",
                     pAdjustMethod = "BH",
                     pvalueCutoff  = 0.01,
                     qvalueCutoff  = 0.05,
                     readable      = TRUE)

dotplot(ck)