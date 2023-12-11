#!/bin/R

#Alberto Gonzalez Delgado
#11/2023
#Centro de Biotecnologia y Genomica de Plantas (UPM-INIA/CSIC)


BiocManager::install("multiWGCNA")
#library(muliWGCNA)
#setwd("C:/Users/admin/Desktop/projects") #Home-desktop
setwd("D:/Alberto/projects") #Lab-desktop

#!/bin/R
#https://bioconductor.org/packages/devel/bioc/vignettes/multiWGCNA/inst/doc/autism_full_workflow.html
#https://www.biostars.org/p/433112/

library(WGCNA)
library(multiWGCNA)
library(dplyr)
library(DESeq2)
library(gridExtra)
library(plotmics)
set.seed(123)
###################################################################################################

#  The multiWGCNA R package is a WGCNA-based procedure designed for RNA-seq datasets with two     #
#  biologically meaningful variables.                                                             # 

###################################################################################################

#1. Read data
data<-read.delim("documents/02.1.raw_counts.tsv")
oscill_data<-read.delim('../oscillation_descriptors.tsv')
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
htmlwidgets::saveWidget(rgl::rglwidget(), "3D_PCA.html")


#5. Normalization --------------------------------------------------------------------
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
                              design = ~1) # Not specifying model because we need  
                                           #the deseq data set to performe variance stabilizing transformation

# Remove all genes with counts >=15 in more than 75% of samples (# samples: 288)
#Suggested by WGCNA on RNAseq FAQ
# I skipped this part because we filtered ny all the genes we want to analyse

#n<-288*0.75
#dds75<-dds[rowSums(counts(dds)>=15)>=n,]
#nrow(dds75) #5765 genes

#Perform variance stabilization
dds_data<-dds[rownames(dds) %in% gene_list,]
dds_norm<-vst(dds_data)
#Get normalized counts (samples as rows and genes as columns)
norm.counts<-as.data.frame(assay(dds_norm))

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
results$diffModExp = runDME(networks[["combined"]], 
                            colData, 
                            p.adjust="fdr", 
                            refCondition="genotype", 
                            testCondition="photoperiod",
                            write=FALSE,
                            plot=FALSE)

test="ANOVA"
mode="PC1"
testCondition="photoperiod"
refCondition="genotype"
datExpr<-networks[["combined"]]@datExpr
modules <- sort(unique(datExpr$dynamicLabels))
modulePrefix <- name(networks[["combined"]])
pval.dfs = list()
pdf("Differential_expression_module_analysis.pdf",width=15)
for (module in modules) {
  moduleGenes = datExpr$X[datExpr$dynamicLabels == module]
  cleanDatExpr = t(cleanDatExpr(datExpr))
  subset = cleanDatExpr[rownames(cleanDatExpr) %in% moduleGenes,]
  if(mode=="Zscore"){
  #Zscore
  mean = rowMeans(subset)
  stdev = apply(subset, 1, sd)
  zscoreMatrix = (subset - mean)/stdev
  zscoreMatrix = na.omit(zscoreMatrix)
  averageExpression = apply(zscoreMatrix, 2, mean)
  moduleExpression = data.frame(Sample = names(averageExpression), 
                                moduleExpression = averageExpression)
  }
  if (mode == "PC1") {
    PC1 = moduleEigengenes(t(subset), colors = rep("Module", 
                                                   length(moduleGenes)), nPC = 1)$eigengenes
    parts <- strsplit(rownames(PC1), "_")
    samples <- sapply(parts, function(x) paste(toupper(c(x[1], x[4],x[3])), collapse = "_"))
    moduleExpression = data.frame(Sample = rownames(PC1), 
                                  moduleExpression = PC1)
    colnames(moduleExpression) = c("Sample", "moduleExpression")
  }
  mergedData = cbind(moduleExpression, colData[match(moduleExpression$Sample, 
                                                     colData$Sample), -1])
  mergedData = mergedData %>% arrange(eval(parse(text = testCondition)), 
                                      eval(parse(text = refCondition)))
  bargraph <- ggplot(data = mergedData, aes(x = factor(Sample, levels = Sample), y = moduleExpression)) + 
    ylab(mode) + 
    geom_bar(aes(fill = moduleExpression), stat = "identity", color = "black", position = position_dodge(9)) + 
    scale_fill_gradient2(name = mode, low = "blue", mid = "white", high = "red", midpoint = 0) + 
    scale_x_discrete(labels = tolower(paste0(substr(mergedData[,3], 0, 3), "_", substr(mergedData[, 4], 0, 3)))) + 
    theme(axis.ticks.x = element_blank(), axis.text.x = element_text(size = 6,angle=90)) +
    labs(title=paste0(module))+
    facet_wrap(~genotype)  

if (test == "ANOVA") {
    pval.df <- performANOVA(moduleExpression, colData, testCondition, 
                            refCondition)
  }
  if (test == "PERMANOVA") {
    requireNamespace("vegan", quietly = TRUE)
    permanova = vegan::adonis(t(subset) ~ colData[[testCondition]] + 
                                colData[[refCondition]] + colData[[testCondition]] * 
                                colData[[refCondition]], method = "euclidean", permutations = 9999)
    Factors = c(testCondition, refCondition, paste0(testCondition, 
                                                    "*", refCondition))
    p.value = c(permanova$aov.tab$`Pr(>F)`[1:3])
    pval.df = data.frame(Factors, p.value)
  }
  boxplot <- ggplot(data = mergedData, aes(x = eval(parse(text = refCondition)), 
                                           y = moduleExpression, color = eval(parse(text = testCondition)))) + 
    labs(title = paste0(module," p: ", paste(pval.df$Factors, signif(pval.df$p.value, 
                                                             2), sep = "=", collapse = ", ")), y = "Expression", 
         x = refCondition) + geom_boxplot(width = 1/length(unique(mergedData[, 
                                                                             refCondition]))) + scale_color_manual(values = c("red","black" ,
                                                                                                                              "blue")) + guides(color = guide_legend(title = testCondition)) + 
    theme(panel.background = element_blank(), axis.line = element_line(colour = "black"), 
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
          plot.title = element_text(hjust = 0.5))
  
 
  print(boxplot)
  print(bargraph)

  
}
dev.off()

#9 Perform the module preservation analysis
enableWGCNAThreads()

results$preservation=iterate(networks[photoperiods], 
                             preservationComparisons, 
                             write=FALSE, 
                             plot=TRUE, 
                             nPermutations=10)



