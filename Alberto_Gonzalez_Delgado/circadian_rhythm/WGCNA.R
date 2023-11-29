#!/bin/R

#Alberto Gonzalez Delgado
#11/2023
#Centro de Biotecnologia y Genomica de Plantas (UPM-INIA/CSIC)

library(WGCNA)
library(DESeq2)
library(GEOquery)
library(tidyverse)
library(CorLevelPlot)
library(gridExtra)
library(ggplot2)
library(rgl)
allowWGCNAThreads()
setwd("C:/Users/admin/Desktop/projects")

#1. Read data ---------------------------------------
data<-read.delim("documents/02.1.raw_counts.tsv")
colnames(data)
data<-data[,c(1,8:ncol(data))]

#Generate metadata 
titles=colnames(data[2:ncol(data)])
genotypes<-c(rep("EID1",24),rep("LE",24),rep("LNK2",24),rep("MM",24))
genotypes<-c(rep(genotypes,3))
photoperiods<-c(rep("LD",96),rep("ND",96),rep("SD",96))
time<-rep(seq(1, 23, by = 2), each = 2, times = 12)

phenoData<-data.frame(samples=titles,genotype=genotypes,photoperiod=photoperiods,replicate=replicates)

head(data)

#2.  Prepare data ---------------------------------------
#data %>% gather(key='samples',value='counts',-gene_id) %>%
#  inner_join(., phenoData,by=c('samples'='samples')) %>%
#  select(1,2,3) %>%
# spread(key='samples',value='counts')
data<- data %>% column_to_rownames(var='gene_id')

#3. QC - outlier detection ----------------------------

gsg<-goodSamplesGenes(t(data)) #Data has to be transposed (genes as columns, samples as rows)
summary(gsg)

##### Both samples and genes
gsg$allOK # TRUE: All the samples and genes have passed the filter
          # FALSE: some outliers have been detected

##### Genes
table(gsg$goodGenes) # TRUE: All the samples and genes have passed the filter
                     # FALSE: some outliers have been detected
##### Samples
table(gsg$goodSamples) # TRUE: All the samples and genes have passed the filter
                       # FALSE: some outliers have been detected

#Remove genes detected as outlier
data<-data[gsg$goodGenes==TRUE,]


#Method 1: Detect outlier samples by hierarchical model 
htree<-hclust(dist(t(data)),method='average') #Data has to be transposed (genes as columns, samples as rows)
pdf("hierarchical_clustering.pdf",width=15,height=8)
plot(htree,cex=0.4)  
dev.off()

#Method 2: PCA

pca<- prcomp(t(data)) #Data has to be transposed (genes as columns, samples as rows)
pca.data<-pca$x #Data about principal components are in x
pca.var<-pca$sdev^2 #Standard variation
pca.var.percent<-round(pca.var/sum(pca.var)*100,digits=2) #Percentage of variance explained
pca.data<-as.data.frame(pca.data)

###Plot 2D PCA
pdf("PCA_samples.pdf",width=10,height=10)
ggplot(pca.data,aes(PC1,PC2))+
  geom_point()+
  geom_text(label=rownames(pca.data),size=1.25)+
  labs(x=paste0('PC1:',pca.var.percent[1],'%'),
       y=paste0('PC2:',pca.var.percent[2],'%'))
dev.off()

### Plot 3D PCA
x <- pca.data$PC1
y <- pca.data$PC2
z <- pca.data$PC3
plot3d(x, y, z, col="red", size=3, 
       xlab = paste0('PC1: ', pca.var.percent[1], '%'), 
       ylab = paste0('PC2: ', pca.var.percent[2], '%'), 
       zlab = paste0('PC3: ', pca.var.percent[3], '%'))
text3d(x, y, z, texts=rownames(pca.data), adj=c(0.5,0.5), color="black", cex=0.6)
text3d(x, y, z, texts=rownames(pca.data), adj=c(0.5,0.5), color="black", cex=0.6)
rgl.postscript("3D_PCA.pdf", fmt = "pdf")



## Remove outliers from further analysis: Separate ND from others?
#photoperiods<-"ND"#c("SD","ND","LD")
#genotypes<-c("LE","LNK2","EID1","MM")
#time<-seq(1,24,2)
#replicates("r1",r2)
#combinations <- expand.grid(photoperiods = photoperiods,
#                            genotypes = genotypes,
#                            time = time,
#                            replicates = replicates)
#samples.to.be.excluded <- paste(combinations$photoperiods, combinations$genotypes, 
#                      paste("CT", combinations$time, sep=""), combinations$replicates, sep="_")
#data.subset<-data[,!(colnames(data) %in% samples.to.be.excluded)]


#4. Normalization --------------------------------------------------------------------
#Generate DESeq2 dataset

#Filter samples excluded 
#colData<-phenoData %>% filter(!row,names(.) %in% in samples.to.be.excluded)
#names(colData)<-gsub(':ch1','',names(colData))
#names(colData)<-gsub('|\\s','_',names(colData))

#rownames and column names identical
colData<-phenoData
rownames(colData)<-colData$samples
all(rownames(colData) %in% colnames(data)) #(data.subset)
colnames(data)


