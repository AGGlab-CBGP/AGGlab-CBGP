#!bin/R

#Alberto Gonzalez Delgado
#11/2023
#Centro de Biotecnologia y Genomica de Plantas (UPM-INIA/CSIC)

library(limma)
library(dplyr)

#```
#vst normalization: In all cases, the transformation is scaled such that for large counts, 
#it becomes asymptotically (for large values) equal to the logarithm to base 2 of normalized counts.
#so it is like it was applied log2 transformation
#```

# Variables
##ZT
ZT<-"03"
##gene
gene_id<-"Solyc03g081240"

#Read data
data<-read.table("documents/02.1.vst_normalized_counts.tsv",sep='\t',header=TRUE)
row_names<-data$gene_id
#Select ZT in which DEG analysis will be executed
filtered_data<-data %>% select(matches(paste0("CT",ZT))) 
rownames(filtered_data)<-row_names 
number_genes<-nrow(filtered_data)
## Experimental design
experimental.design<-model.matrix(~ -1+factor(c(1,1,2,2,3,3,4,4,5,5,6,6,7,7,8,8,9,9,10,10,11,11,12,12)))
colnames(experimental.design)<-c("LD_EID1","LD_LE","LD_LNK2","LD_MM","ND_EID1","ND_LE","ND_LNK2","ND_MM","SD_EID1","SD_LE","SD_LNK2","SD_MM")

#linear fit (calculate mean between replicates)
linear.fit<-lmFit(filtered_data,experimental.design)

#Specify contrasts: 1st treatment, 2nd control
## If data in log2 -> -
## Else -> /

## The names must be EXACTLY the same indicated in experimental.design
contrast.matrix<-makeContrasts(LD_MM - LD_LE,
                               ND_MM - ND_LE,
                               SD_MM - SD_LE,
                               levels=c("LD_EID1","LD_LE","LD_LNK2","LD_MM","ND_EID1","ND_LE","ND_LNK2","ND_MM","SD_EID1","SD_LE","SD_LNK2","SD_MM"))
#Calculo FC
contrast.linear.fit<-contrasts.fit(linear.fit,contrast.matrix)

#Calculo qvalue
contrast.results<-eBayes(contrast.linear.fit)
 

#LD
LD<-topTable(contrast.results,number=number_genes,coef=1,sort.by="logFC")

#ND
ND<-topTable(contrast.results,number=number_genes,coef=2,sort.by="logFC")

#SD
SD<-topTable(contrast.results,number=number_genes,coef=3,sort.by="logFC")


#Print results of a gene
subset(LD, grepl(paste0(gene_id), rownames(LD)))
subset(ND, grepl(paste0(gene_id), rownames(ND)))
subset(SD, grepl(paste0(gene_id), rownames(SD)))       
