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


#4. Normalization --------------------------------------------------------------------
#Generate DESeq2 dataset

#Filter samples excluded 
#colData<-phenoData %>% filter(!row,names(.) %in% in samples.to.be.excluded)
#names(colData)<-gsub(':ch1','',names(colData))
#names(colData)<-gsub('|\\s','_',names(colData))

#rownames and column names identical
colData<-phenoData
rownames(colData)<-colData$Sample
all(rownames(colData) %in% colnames(data)) #Are all the conditions? #(data.subset)
all(rownames(colData) == colnames(data)) #They must to be in the same order

# Create dds
dds <- DESeqDataSetFromMatrix(countData = data,
                              colData=colData,
                              design = ~1) # Not specifying model because we need  #the deseq data set to performe variance stabilizing transformation

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
norm.counts<-as.data.frame(norm.counts)
norm.counts <- norm.counts %>%
  mutate(gene_id=rownames(norm.counts)) %>%
  filter(gene_id %in% gene_list) %>%
  select(-c(gene_id))


#5. Perform network construction, module eigengene calculation, module-trait correlation --------------------

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

#6. Compare modules by overlap across conditions
pdf("multiWGCNA_results.pdf")
results=list()
results$overlaps=iterate(networks, overlapComparisons, plot=TRUE)
dev.off()

#Visualize matches
head(results$overlaps$LD_vs_SD$bestMatches)

#7. Perform differential module expression analysis

