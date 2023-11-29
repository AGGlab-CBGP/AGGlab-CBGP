#!bin/R

#Alberto Gonzalez Delgado
#11/2023
#Centro de Biotecnologia y Genomica de Plantas (UPM-INIA/CSIC)

library(DESeq2)

#Importar datos
library(tidyverse)
data<-read.csv("02.1.raw_counts.tsv",sep='\t')

# Variables
##ZT
ZT<-"03"
##gene id
gene_id<-"Solyc05g053850"


##Filtrar por ZT(disminuye tiempo de computacion)
filtered_data<-data %>% select(matches(paste0("CT",ZT)))
rownames(filtered_data)<-data$gene_id

#Crear una colData (nombre de fila= columnaID en counts_data, columna=tratamiento o condition)
colData<-data.frame("condition"=gsub("_CT[0-1][0-9]_r[0-9]","",names(filtered_data)))
rownames(colData)<-names(filtered_data)

#Comprobando que las condiciones están en colData.Tienen que estar en el mismo orden o DESeq2 da error
all(colnames(filtered_data) %in% rownames(colData))
all(colnames(filtered_data)==rownames(colData))

#Generar DESeqDataset object
dds<-DESeqDataSetFromMatrix(countData=filtered_data, colData=colData,
                       design= ~condition)

#Prefiltro: Eliminar genes con menos de 10 counts en total
keep<-rowSums(counts(dds))>=10
dds<-dds[keep,]#Factor level
dds$condition<-relevel(dds$condition,ref="ND_LE")

#DESeq
dds<-DESeq(dds)
res<-results(dds,alpha=0.01)summary(res)

#Extraer resultados de una comparación
res_MM_ND<-results(dds,contrast = c("condition","ND_MM","ND_LE"),alpha=0.05)
summary(res_MM_ND)

#MAplot
plotMA(res_MM_ND)

#Extraer resultados de un gen 
subset(res_MM_ND,grepl(paste0(gene_id),rownames(res_MM_ND)))

#Extraer listas de genes

##Upregulated
indexes<-indices <- which(res_MM_ND$log2FoldChange > 0 & res_MM_ND$padj < 0.05)
upregulated<-rownames(res_MM_ND)[indexes]

##Downregulated
indexes<-indices <- which(res_MM_ND$log2FoldChange < 0 & res_MM_ND$padj < 0.05)
downregulated<-rownames(res_MM_ND)[indexes]
