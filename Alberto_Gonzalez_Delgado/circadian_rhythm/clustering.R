#!/bin/R

#Alberto Gonzalez Delgado
#11/2023
#Centro de Biotecnologia y Genomica de Plantas (UPM-INIA/CSIC)


BiocManager::install("multiWGCNA")
#library(muliWGCNA)
setwd("C:/Users/admin/Desktop/projects")


#!/bin/R
#https://bioconductor.org/packages/devel/bioc/vignettes/multiWGCNA/inst/doc/autism_full_workflow.html
#https://www.biostars.org/p/433112/

library(WGCNA)
library(multiWGCNA)
library(dplyr)
library(DESeq2)
library(gridExtra)
###################################################################################################

#  The multiWGCNA R package is a WGCNA-based procedure designed for RNA-seq datasets with two     #
#  biologically meaningful variables.                                                             # 

###################################################################################################

#1. Read data
data<-read.delim("documents/02.1.raw_counts.tsv")
oscill_data<-read.delim('oscillation_descriptors.tsv')
gene_list<-oscill_data$CycID
head(data)

rownames(data)<-data$gene_id
data<-data[c(8:ncol(data))]



#Generate metadata -----------------------------------
titles=colnames(data[1:ncol(data)])
genotypes<-c(rep("EID1",24),rep("LE",24),rep("LNK2",24),rep("MM",24))
genotypes<-c(rep(genotypes,3))
photoperiods<-c(rep("LD",96),rep("ND",96),rep("SD",96))
time<-rep(seq(1, 23, by = 2), each = 2, times = 12)
replicates<-rep(c("r1","r2"),144)

phenoData<-data.frame(Sample=titles,genotype=genotypes,photoperiod=photoperiods,replicate=replicates)
rownames(phenoData)<-phenoData$Sample
head(phenoData)
head(data)

#Define conditions for traits: genotype and photoperiod
genotypes<-unique(phenoData[,2])
photoperiods<-unique(phenoData[,3])

#4. Exploratory data analysis -------------------------------------------------------

# Hierarchical clustering 
htree<-hclust(dist(t(data)),method='average') #Data has to be transposed (genes as columns, samples as rows)
pdf("hierarchical_clustering.pdf",width=15,height=8)
plot(htree,cex=0.4)  
dev.off()

#PCA

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


#5. Normalization --------------------------------------------------------------------
#Generate DESeq2 dataset

#Filter samples excluded 
#colData<-phenoData %>% filter(!row,names(.) %in% in samples.to.be.excluded)
#names(colData)<-gsub(':ch1','',names(colData))
#names(colData)<-gsub('|\\s','_',names(colData))

#rownames and column names identical
colData<-phenoData
 
dev.off()rownames(colData)<-colData$Sample
all(rownames(colData) %in% colnames(data)) #Are all the conditions? #(data.subset)
all(rownames(colData) == colnames(data)) #They must to be in the same order

# Create dds
dds <- DESeqDataSetFromMatrix(countData = data,
                              colData=colData,
                              design = ~1) # Not specifying model because we need  
                                           #the deseq data set to performe variance stabilizing transformation

# Remove all genes with counts >=15 in more than 75% of samples (# samples: 288)
#Suggested by WGCNA on RNAseq FAQ
# I skipped this part because we filtered ny all the genes we want to analyse

#n<-288*0.75
#dds75<-dds[rowSums(counts(dds)>=15)>=n,]
#nrow(dds75) #5765 genes

#Perform variance stabilization
#dds_norm<-vst(dds75)
dds_norm<-vst(dds)
#Get normalized counts (samples as rows and genes as columns)
norm.counts<-as.data.frame(assay(dds_norm))

# Remove all genes with counts >=15 in more than 75% of samples (# samples: 288)
#Suggested by WGCNA on RNAseq FAQ
# I skipped this part because we filtered by all the genes we want to analyse

#n<-288*0.75
#dds75<-dds[rowSums(counts(dds)>=15)>=n,]
#nrow(dds75) #5765 genes

norm.counts <- norm.counts %>%
  mutate(gene_id=rownames(norm.counts)) %>%
  filter(gene_id %in% gene_list) %>%
  select(-c(gene_id))


#6. Perform network construction, module eigengene calculation, module-trait correlation --------------------


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
  labs(X='Power',y='Cale free topology model fit, signed R^2')+
  theme_gray()

p2<-ggplot(sft.data,aes(Power,mean.k., label=Power))+
  geom_point()+
  geom_text(nudge_y=0.1)+
  labs(x="Power",y='Mean Connectivity') +
  theme_gray()

grid.arrange(p1,p2,ncol=1) #Normally select higher than R2 0.8
#and low mean connectivity

class(data)
networks = constructNetworks(norm.counts, phenoData, genotypes, photoperiods,
                             networkType = "signed", power = 12,
                             minModuleSize = 40, maxBlockSize = 25000,
                             reassignThreshold = 0, minKMEtoStay = 0.7,
                             mergeCutHeight = 0.10, numericLabels = TRUE,
                             pamRespectsDendro = FALSE, verbose=3)

#7. Compare modules by overlap across conditions
pdf("multiWGCNA_results.pdf")
results=list()
results$overlaps=iterate(networks, overlapComparisons, plot=TRUE)
dev.off()

#Visualize matches
head(results$overlaps$LD_vs_SD$bestMatches)

#8. Perform differential module expression analysis

