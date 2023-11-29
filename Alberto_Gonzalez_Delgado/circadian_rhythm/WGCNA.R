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
all(rownames(colData) %in% colnames(data)) #Are all the conditions? #(data.subset)
all(rownames(colData) == colnames(data)) #They must to be in the same order

# Create dds
dds <- DESeqDataSetFromMatrix(countData = data,
                              colData=colData,
                              design = ~1) # Not specifying model because we need 
                                           #the deseq data set to performe variance stabilizing transformation

# Remove all genes with counts >=15 in more than 75% of samples
#Suggested by WGCNA on RNAseq FAQ

n<-288*0.75
dds75<-dds[rowSums(counts(dds)>=15)>=n,]
nrow(dds75) #17711 genes

#Perform variance stabilization
dds_norm<-vst(dds75)

#Get normalized counts (samples as rows and genes as columns)
norm.counts<-assay(dds_norm) %>% t()

#4. Network Construction ------------------------
#Choose a set of soft-thresholding powers
power<-c(c(1:10),seq(12,50,by=2))

# Call the network topology analysis function
sft<-pickSoftThreshold(norm.counts,
                  powerVector = power,
                  networkType = "signed",
                  verbose=5)

#Select scaled-free topology
sft.data<-sft$fitIndices #We want the network with maximun Rsquared and minimum mean conectivity

p1<-ggplot(sft.data,aes(Power,SFT.R.sq,label=Power))+
  geom_point()+
  geom_text(nudge_y=0.1)+
  geom_hline(yintercept = 0.8, color='red')+
  labs(X='Power',y='Cale free topology model fit, signed R^2')+
  theme_gray()

p2<-ggplot(sft.data,aes(Power,mean.k., label=Power))+
  geom_point()+
  geom_text(nudge_y=0.1)+
  labs(x="Power",y='Mean Connectivity')+
  theme_gray()

grid.arrange(p1,p2,ncol=1) #Normally select higher than R2 0.8
                           #and low mean connectivity
#Create adjacency matrix
norm.counts[]<-sapply(norm.counts,as.numeric)
soft.power<-14 #Theshold selected from previous plot
temp_cor<-cor #Assign correlation function to a temporal variable so we use WGCNA correlation function
cor<-WGCNA::cor #Assign correlation function to WGCNA correlation function

#Memory estimate w.r.t blocksize
bwnet<-blockwiseModules(norm.counts, 
                 maxBlockSize = 15000, # How many genes include in one block keeping in mind the memory
                                # the system has access to (most desktop:5000, 4Gb RAM: 8000-10000,
                                # 16Gb RAM: 20000, 32Gb RAM:30000).
                  TOMType = "signed",
                 power=soft.power,
                 mergeCutHeight=0.25,
                numericLabels = FALSE,
                randomSeed = 24,
                verbose=3)
cor<-temp_cor

#5. Module Eigengenes ------------------

module_eigengenes<-bwnet$MEs
head(module_eigengenes)

#Name of each module
table(bwnet$colors)


#Plot dendrogram
mergedColors = labels2colors(bwnet$colors)
plotDendroAndColors(bwnet$dendrograms[[1]], mergedColors[bwnet$blockGenes[[1]]], "Module colors", dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, guideHang = 0.05)

unmergedColors=labels2colors(bwnet$unmergedColors)
plotDendroAndColors(bwnet$dendrograms[[1]], cbind(mergedColors[bwnet$blockGenes[[1]]],unmergedColors[bwnet$blockGenes[[1]]]),
                    c("unmerged", "merged"),
                    dendroLabels = FALSE,
                    addGuide = TRUE,
                    hang= 0.03,
                    guideHang = 0.05)

table(unname(bwnet$colors))
table(unname(bwnet$unmergedColors))#Tutorial: more colors in unmerged, meaning that some clusters were grouped into a single one, so merged will be used




