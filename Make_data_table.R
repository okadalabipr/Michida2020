#This code makes data table which is required for the downstream analyses such as mathmatical model analysis.
library(scde)
library(tidyverse)

#Import gene counts data and conduct the same cleaning as SCDE.R did--------------------------------------------------
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


dose_0 <- inputff[,grep("^IgM_00000_*", colnames(inputff))]
dose_0.01 <- inputff[,grep("^IgM_00010_*", colnames(inputff))]
dose_0.1 <- inputff[,grep("^IgM_00100_*", colnames(inputff))]
dose_1 <- inputff[,grep("^IgM_01000_*", colnames(inputff))]
dose_10 <- inputff[,grep("^IgM_10000_*", colnames(inputff))]

cd_0.01 <- clean.counts(cbind(dose_0.01,dose_0), min.lib.size=1000, min.reads=1, min.detected=1)
cd_0.1 <- clean.counts(cbind(dose_0.1,dose_0), min.lib.size=1000, min.reads=1, min.detected=1)
cd_1 <- clean.counts(cbind(dose_1,dose_0), min.lib.size=1000, min.reads=1, min.detected=1)
cd_10 <- clean.counts(cbind(dose_10,dose_0), min.lib.size=1000, min.reads=1, min.detected=1)

cd_0.01$Geneid <- rownames(cd_0.01)
cd_0.1$Geneid <- rownames(cd_0.1)
cd_1$Geneid <- rownames(cd_1)
cd_10$Geneid <- rownames(cd_10)

target <- intersect(intersect(intersect(cd_0.01$Geneid, cd_0.1$Geneid), cd_1$Geneid), cd_10$Geneid)

tmp_0.01 <- cd_0.01[match(target, cd_0.01$Geneid),]
tmp_0.1 <- cd_0.1[match(target, cd_0.1$Geneid),]
tmp_1 <- cd_1[match(target, cd_1$Geneid),]
tmp_10 <- cd_10[match(target, cd_10$Geneid),]

Last_IgM_00010 <- max(grep("^IgM_00010_*", colnames(tmp_0.01)))
Last_IgM_00100 <- max(grep("^IgM_00100_*", colnames(tmp_0.1)))
Last_IgM_01000 <- max(grep("^IgM_01000_*", colnames(tmp_1)))
Last_IgM_10000 <- max(grep("^IgM_10000_*", colnames(tmp_10)))

dose_0 <- tmp_0.01[,-(1:Last_IgM_00010)]
dose_0.01 <- tmp_0.01[,-((Last_IgM_00010+1):ncol(tmp_0.01))]
dose_0.1 <- tmp_0.1[,-((Last_IgM_00100+1):ncol(tmp_0.1))]
dose_1 <- tmp_1[,-((Last_IgM_01000+1):ncol(tmp_1))]
dose_10 <- tmp_10[,-((Last_IgM_10000+1):ncol(tmp_10))]

dose_0 <- dose_0[,c(ncol(dose_0),1:(ncol(dose_0)-1))]

cd <- cbind(dose_0,dose_0.01,dose_0.1,dose_1,dose_10)

cd$Length <- inputff$Length[match(target, inputff$Geneid)]

#Calculate mean TPM values for single cells for each gene-------------------------------------------------------
TPM = cd

#TPM
#divide by each transcript length
len = TPM$Length


for(i in 2:(ncol(TPM)-1)){
  TPM[,i] = (TPM[,i]/len)*1000
}

#calculate library size
SUM = apply(TPM[,c(2:(ncol(TPM)-1))],2,sum)

#divide by library size 
for(i in 1:length(SUM)){
  TPM[,i+1] = (TPM[,i+1]/SUM[i])*1000000
}


rownames(TPM) <- seq(1,nrow(TPM),1)

#write.table(TPM,"cleaned_TPM.txt",row.names = F,quote = F, sep = "\t") #This file is required for downstream analysis.


TPM_mean_0 <- TPM[,2:max(grep("^IgM_00000_*", colnames(TPM)))] %>% rowMeans()

TPM_mean_0.01 <- TPM[,min(grep("^IgM_00010_*", colnames(TPM))):max(grep("^IgM_00010_*", colnames(TPM)))] %>% rowMeans()

TPM_mean_0.1 <- TPM[,min(grep("^IgM_00100_*", colnames(TPM))):max(grep("^IgM_00100_*", colnames(TPM)))] %>% rowMeans()

TPM_mean_1 <- TPM[,min(grep("^IgM_01000_*", colnames(TPM))):max(grep("^IgM_01000_*", colnames(TPM)))] %>% rowMeans()

TPM_mean_10 <- TPM[,min(grep("^IgM_10000_*", colnames(TPM))):max(grep("^IgM_10000_*", colnames(TPM)))] %>% rowMeans()

TPM_mean <- data.frame(geneid = TPM$Geneid, Length = TPM$Length, TPM_0 = TPM_mean_0, 
                       TPM_0.01 = TPM_mean_0.01, TPM_0.1 = TPM_mean_0.1,
                       TPM_1 = TPM_mean_1, TPM_10 = TPM_mean_10)

#Join peak count data and motif count data--------------------------------------------------------------------------
#These four files are obtained using RelA_peak_motif_count.R. Please specify their path depending on your environment.
SE_sg_peak <- read.csv("SE_sg_peak_count.csv",stringsAsFactors = F)
SE_sg_motif <- read.csv("SE_sg_motif_count.csv", stringsAsFactors = F)
TE_sg_peak <- read.csv("TE_sg_peak_count.csv",stringsAsFactors = F)
TE_sg_motif <- read.csv("TE_sg_motif_count.csv", stringsAsFactors = F)

#Add column of motif num to SE_sg_peak and TE_sg_peak
v <- c()
for(i in 1:nrow(SE_sg_peak)){
  index <- which(SE_sg_peak$ID[i] == SE_sg_motif$ID)
  if(length(index) != 0){
    v <- append(v,index)
  }
}

SE_sg_peak$motif_num <- SE_sg_motif$peak_num[v]

v <- c()
for(i in 1:nrow(TE_sg_peak)){
  index <- which(TE_sg_peak$ID[i] == TE_sg_motif$ID)
  if(length(index) != 0){
    v <- append(v,index)
  }
}

TE_sg_peak$motif_num <- TE_sg_motif$peak_num[v]

#Convert log2FC into FC
SE_sg_peak$FC_0.01 <- 2^SE_sg_peak$log2FC_0.01
SE_sg_peak$FC_0.1 <- 2^SE_sg_peak$log2FC_0.1
SE_sg_peak$FC_1 <- 2^SE_sg_peak$log2FC_1
SE_sg_peak$FC_10 <- 2^SE_sg_peak$log2FC_10

TE_sg_peak$FC_0.01 <- 2^TE_sg_peak$log2FC_0.01
TE_sg_peak$FC_0.1 <- 2^TE_sg_peak$log2FC_0.1
TE_sg_peak$FC_1 <- 2^TE_sg_peak$log2FC_1
TE_sg_peak$FC_10 <- 2^TE_sg_peak$log2FC_10

#Get mean TPM of each gene
v <- c()
for(i in 1:nrow(SE_sg_peak)){
  tmp <- match(SE_sg_peak$gene[i],TPM_mean$geneid)
  v <- append(v,tmp)
}

TPM_mean_SE_sg <- TPM_mean[v,]

SE_sg_peak <- cbind(SE_sg_peak,TPM_mean_SE_sg[,3:ncol(TPM_mean_SE_sg)])

rownames(SE_sg_peak) <- seq(1,nrow(SE_sg_peak),1)

v <- c()
for(i in 1:nrow(TE_sg_peak)){
  tmp <- match(TE_sg_peak$gene[i],TPM_mean$geneid)
  v <- append(v,tmp)
}

TPM_mean_TE_sg <- TPM_mean[v,]

TE_sg_peak <- cbind(TE_sg_peak,TPM_mean_TE_sg[,3:ncol(TPM_mean_TE_sg)])

rownames(TE_sg_peak) <- seq(1,nrow(TE_sg_peak),1)

#These files are required for downstream analyses.
write.csv(SE_sg_peak,"cleaned_Gained_SE.csv",row.names = F)
write.csv(TE_sg_peak,"cleaned_Gained_TE.csv",row.names = F)
