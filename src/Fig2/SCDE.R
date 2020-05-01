#I'm sorry for this very messy source code...
library(scde)
inputff <- read.table("path to your count data obtained with featureCounts", header=T, stringsAsFactors=F,sep="\t")
rownames(inputff) <- inputff$Geneid

ref = c("RamDA_00000",
        "RamDA_00010",
        "RamDA_00100",
        "RamDA_01000",
        "RamDA_10000")

v_min = c()
v_max = c()

for (i in ref){
  v_min = c(v_min,min(grep(i,colnames(inputff))))
  v_max = c(v_max,max(grep(i,colnames(inputff))))
} 

colnames(inputff) <- gsub("RamDA_00000", "", colnames(inputff))
colnames(inputff) <- gsub("RamDA_00010", "", colnames(inputff))
colnames(inputff) <- gsub("RamDA_00100", "", colnames(inputff))
colnames(inputff) <- gsub("RamDA_01000", "", colnames(inputff))
colnames(inputff) <- gsub("RamDA_10000", "", colnames(inputff))
colnames(inputff) <- gsub("_Aligned.sortedByCoord.out.bam", "", colnames(inputff))
colnames(inputff)[v_min[1]:v_max[1]] <- gsub("_", "IgM_00000_", colnames(inputff)[v_min[1]:v_max[1]])
colnames(inputff)[v_min[2]:v_max[2]] <- gsub("_", "IgM_00010_", colnames(inputff)[v_min[2]:v_max[2]])
colnames(inputff)[v_min[3]:v_max[3]] <- gsub("_", "IgM_00100_", colnames(inputff)[v_min[3]:v_max[3]])
colnames(inputff)[v_min[4]:v_max[4]] <- gsub("_", "IgM_01000_", colnames(inputff)[v_min[4]:v_max[4]])
colnames(inputff)[v_min[5]:v_max[5]] <- gsub("_", "IgM_10000_", colnames(inputff)[v_min[5]:v_max[5]])

sg <- factor( c(rep("IgM_00000", length( grep("^IgM_00000_*", colnames(inputff)) ) ),
                rep("IgM_00010", length( grep("^IgM_00010_*", colnames(inputff)) ) ),
                rep("IgM_00100", length( grep("^IgM_00100_*", colnames(inputff)) ) ),
                rep("IgM_01000", length( grep("^IgM_01000_*", colnames(inputff)) ) ),
                rep("IgM_10000", length( grep("^IgM_10000_*", colnames(inputff)) ) )
) )

conc <- c(0.01, 0.1, 1, 10)

for (c in conc) {
  
  tag <- sprintf("^IgM_%05d_*",c*1000)
  ctl <- inputff[,grep("^IgM_00000_*", colnames(inputff))]
  trt <- inputff[,grep(tag, colnames(inputff))]
  
  df <- cbind( trt, ctl )
  
  cd <- clean.counts(df, min.lib.size=1000, min.reads=1, min.detected=1)
  
  treat = length(grep(sprintf("^IgM_%05d_*",c*1000),colnames(cd)))
  
  control = length(grep("^IgM_00000_*",colnames(cd)))
  
  sg <- as.factor(c(rep("00Trt", treat), rep("01Ctl", control)))
  
  v = c()
  for (i in colnames(df)){
    if (length(grep(i,colnames(cd))) == 1){
      v = c(v,i)
    }else{
    }
  }
  

  names(sg) <- v
  table(sg)
  
  
  
  # error model
  o.ifm <- scde.error.models(counts=cd, groups=sg, n.cores=1,
                             threshold.segmentation=T,
                             save.crossfit.plots=F,
                             save.model.plots=T, verbose=1)
  saveRDS(o.ifm, paste0("error_model_0_", c, ".obj"))
  
  # filter out cells that don't show positive correlation with
  # the expected expression magnitudes (very poor fits)
  valid.cells <- o.ifm$corr.a > 0
  table(valid.cells)
  
  o.ifm <- o.ifm[valid.cells, ]
  
  # estimate gene expression prior
  o.prior <- scde.expression.prior(models=o.ifm, counts=cd,
                                   length.out=400, show.plot=FALSE)
  
  
  # define two groups of cells
  nctl <- length( grep("^IgM_00000_*", rownames(o.ifm)) )
  ntrt <- nrow(o.ifm) - nctl
  groups <- factor( c( rep("00Trt", ntrt), rep("01Ctl", nctl) ) )
  names(groups) <- row.names(o.ifm)
  # run differential expression tests on all genes.
  ediff <- scde.expression.difference(o.ifm, cd, o.prior, groups=groups, return.posteriors = T,
                                      n.randomizations=100, n.cores=1, verbose=1)
  
  
  # top upregulated genes (tail would show top downregulated ones)
  test_results = ediff$results
  head(test_results[order(test_results$Z, decreasing=TRUE), ])
  
  
  # output data
  # scde.test.gene.expression.difference("Cd83", models = o.ifm, counts = cd,
  #                                      return.details = F, prior = o.prior)
  
  # calculate p-value
  p.values <- 2*pnorm(abs(test_results$Z),lower.tail=F)
  p.values.adj <- 2*pnorm(abs(test_results$cZ),lower.tail=F) 
  
  # significant genes 
  significant.genes <- which(p.values.adj<10)
  length(significant.genes)
  
  ord <- order(p.values.adj[significant.genes]) # order by p-value
  de <- cbind(ediff$results[significant.genes,1:3],p.values.adj[significant.genes])[ord,]
  colnames(de) <- c("Lower bound","log2 fold change","Upper bound","p-value")
  write.csv(de, paste0("quantified_", c, ".csv"), row.names=T)
}

# merging data
A <- read.csv( sprintf("quantified_%s.csv",conc[1]),  header=T, stringsAsFactors=F )
B <- read.csv( sprintf("quantified_%s.csv",conc[2]),  header=T, stringsAsFactors=F )
C <- read.csv( sprintf("quantified_%s.csv",conc[3]),  header=T, stringsAsFactors=F )
D <- read.csv( sprintf("quantified_%s.csv",conc[4]),  header=T, stringsAsFactors=F )

target <- intersect(intersect(intersect(A$X, B$X), C$X), D$X)

A2 <- A[match(target, A$X),]
B2 <- B[match(target, B$X),]
C2 <- C[match(target, C$X),]
D2 <- D[match(target, D$X),]

E <- cbind(A2, B2[,-1], C2[,-1], D2[,-1])
colnames(E)[2:17] <- c("Lower_0.01", "log2FC_0.01", "Upper_0.01", "pvalue_0.01",
                       "Lower_0.1", "log2FC_0.1", "Upper_0.1", "pvalue_0.1",
                       "Lower_1", "log2FC_1", "Upper_1", "pvalue_1",
                       "Lower_10", "log2FC_10", "Upper_10", "pvalue_10"
)

write.csv(E, "RamDA.csv", row.names = F) #This file is required for the downstream analysis.