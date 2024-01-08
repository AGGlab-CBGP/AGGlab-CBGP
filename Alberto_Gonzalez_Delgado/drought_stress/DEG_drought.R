#!/bin/R

#Alberto Gonzalez Delgado
#12/2023
#Centro de Biotecnologia y Genomica de Plantas (UPM-INIA/CSIC)

library(DESeq2)
library(tidyverse) 
library(ggplot2)
library(pheatmap)
library(clusterProfiler)
library(AnnotationForge)
library(igraph)
library(WGCNA)
library(cluster)
allowWGCNAThreads()
setwd("C:/Users/admin/Desktop/projects/0.3 DEG drought stress")

# Import data ----------------------------------------------------------------
exp_data<-read.delim("05.1.feature_counts.counts.tsv",sep='\t')
rownames(exp_data)<-exp_data$gene_id
exp_data<-exp_data[c(8:ncol(exp_data))]
col_names <- names(exp_data)

annotation<-read.delim("GO_annotation.txt",header=TRUE)
upregulated<- noquote(upregulated)
annotation<-split(annotation$GO,annotation$GID)
enricher(downregulated, TERM2GENE=annotation)

#################### Split data by timepoint ##################################
data_5 <- data.frame(matrix(ncol=0,nrow=nrow(exp_data)))
data_8 <- data.frame(matrix(ncol=0,nrow=nrow(exp_data)))
data_11<- data.frame(matrix(ncol=0,nrow=nrow(exp_data)))
for (name in col_names) {
  time_point <- as.numeric(substr(strsplit(name,"\\.")[[1]][3], 2, nchar(strsplit(name,"\\.")[[1]][3])))
  print(time_point)
  if (time_point == 5) {
    data_5 <- cbind(data_5, setNames(exp_data[, name, drop = FALSE], name))
  } else if (time_point == 8) {
    data_8 <- cbind(data_8, setNames(exp_data[, name, drop = FALSE], name))
  } else if (time_point == 11) {
    data_11 <- cbind(data_11, setNames(exp_data[, name, drop = FALSE], name))
  }
}
##############################################################################
######################## Gene co-expression network ##########################
##############################################################################

#1. Generate metadata ----------------------------------------------------------
names_<-names(exp_data)
probe <- c()
genotype <- c()
treatment <- c()
for (name in names_) {
  split <- strsplit(name, "\\.")[[1]]
  probe <- c(probe, substr(split[1], 2, nchar(split[1])))
  genotype <- c(genotype, ifelse(substr(split[2], 1, 1) == "M", "MM", "Penelli"))
  treatment <- c(treatment, ifelse(substr(split[2], 2, 2) == "D", "Drought", "Control"))
}
colData<-data.frame(Sample=probe,Treatment=treatment,Genotype=genotype)
rownames(colData)<-names(exp_data)

table(genotype) #MM 16 Penelli 16
table(treatment) #Control 16 Drought 16

#Check whether the metadata is generated correctly
all(colnames(exp_data) %in% rownames(colData))
all(colnames(exp_data)==rownames(colData))

#2. Generate DESeq dataset -----------------------------------------------------
dds<-DESeqDataSetFromMatrix(countData=exp_data, colData=colData,
                              design= ~1) # Not specifying model because we need 
                                          #the deseq data set to performe variance stabilizing transformation

#Pre-filtering: remove genes with less than 10 counts in total
keep<-rowSums(counts(dds))>=10
dds<-dds[keep,]

#Variance stabilizing transformation
vst_ <- vst(dds, blind=FALSE)
#7. Gene coexpression network -------------------------------------------------
## Generate network
coexpression.data.t <- t(assay(vst_))
gene.correlation <- cor(coexpression.data.t)
adjacency <- (gene.correlation > 0.91) & (gene.correlation < 1)
rm(gene.correlation)
gene.coexpression.network <- graph.adjacency(adjacency, mode="undirected")
#write.graph(gene.coexpression.network,file="gene_coexpression_network_083.gml",format="gml")
rm(adjacency)
network.degree.distribution <- degree.distribution(gene.coexpression.network)
rm(gene.coexpression.network)
#Asses free-scale fit
degree.histogram <- hist(network.degree.distribution,freq=FALSE,col="blue",xlab="Node degree",
                         ylab="Probability",main="Degree distribution")
fit.scale.free <- power.law.fit(network.degree.distribution)
fit.scale.free[["KS.p"]] #High value: it adjust to a negative potencial

degree.frequencies <- table(network.degree.distribution)
degree.frequencies.no.0 <- degree.frequencies[-1]
log10.degrees.frequencies <- log10(degree.frequencies.no.0)
log10.node.degrees <- log10(as.numeric(names(degree.frequencies.no.0)))
lm.r <- lm(log10.degrees.frequencies ~ log10.node.degrees)
summary(lm.r) 
# 98: 0.870402  | Multiple R-squared:  0.3942, p-value: 0.0004551
# 97: 0.9997389 | Multiple R-squared:  0.4234, p-value: 9.753e-06
# 96: 0.5640501 | Multiple R-squared:  0.5034, p-value: 1.435e-07
# 95: 0.9999849 | Multiple R-squared:  0.5273, p-value: 5.197e-09
# 94: 0.9561293 | Multiple R-squared:  0.5005, p-value: 1.09e-09
# 93: 0.7877362 | Multiple R-squared:  0.4513, p-value: 1.416e-11
# 92: 0.9898965 | Multiple R-squared:  0.4659, p-value: 2.275e-12
# 91: 0.9975117 | Multiple R-squared:  0.5566, p-value: 4.337e-13
# 90: 0.8245205 | Multiple R-squared:  0.5537, p-value: 8.09e-13

# 84: 0.9731735 | Multiple R-squared:  0.5904, p-value: 1.254e-13
# 83: 0.9381306 | Multiple R-squared:  0.6408, p-value: 2.517e-14
# 82: 0.9168867 | Multiple R-squared:  0.5993, p-value: 9.971e-14
# 81: 0.9299475 | Multiple R-squared:  0.6093, p-value: 1.192e-13
# 80: 0.9272208 | Multiple R-squared:  0.6055, p-value: 6.186e-14
gene.correlation <- cor(coexpression.data.t)
similarity.matrix <- 1 - gene.correlation
hierarchical.clustering <- hclust(as.dist(similarity.matrix),method="average") 

hclust.2 <- cutree(hierarchical.clustering,k=2)
hclust.3 <- cutree(hierarchical.clustering,k=3)
hclust.4 <- cutree(hierarchical.clustering,k=4)
hclust.5 <- cutree(hierarchical.clustering,k=5)
hclust.6 <- cutree(hierarchical.clustering,k=6)

sil.hclust.2 <- silhouette(hclust.2,dist=similarity.matrix)
sil.hclust.3 <- silhouette(hclust.3,dist=similarity.matrix)
sil.hclust.4 <- silhouette(hclust.4,dist=similarity.matrix)
sil.hclust.5 <- silhouette(hclust.5,dist=similarity.matrix)
sil.hclust.6 <- silhouette(hclust.6,dist=similarity.matrix)

hclust.sil.values <- c(summary(sil.hclust.2)[["avg.width"]],summary(sil.hclust.3)[["avg.width"]],
                       summary(sil.hclust.4)[["avg.width"]],summary(sil.hclust.5)[["avg.width"]],
                       summary(sil.hclust.6)[["avg.width"]])


pam.2 <- pam(as.dist(similarity.matrix),k=2,diss=TRUE)
pam.3 <- pam(as.dist(similarity.matrix),k=3,diss=TRUE)
pam.4 <- pam(as.dist(similarity.matrix),k=4,diss=TRUE)
pam.5 <- pam(as.dist(similarity.matrix),k=5,diss=TRUE)
pam.6 <- pam(as.dist(similarity.matrix),k=6,diss=TRUE)

sil.pam.2 <- silhouette(pam.2)
sil.pam.3 <- silhouette(pam.3)
sil.pam.4 <- silhouette(pam.4)
sil.pam.5 <- silhouette(pam.5)
sil.pam.6 <- silhouette(pam.6)

pam.sil.values <- c(summary(sil.pam.2)[["avg.width"]],summary(sil.pam.3)[["avg.width"]],
                    summary(sil.pam.4)[["avg.width"]],summary(sil.pam.5)[["avg.width"]],
                    summary(sil.pam.6)[["avg.width"]])
pdf("Silhouette_score_clustering.pdf")
plot(2:6,pam.sil.values,type="o",col="blue",pch=0,ylim=c(0.0,0.8),xlab="Number of clusters",ylab="Silhouette",lwd=3)
lines(2:6,hclust.sil.values,type="o",col="red",pch=1,xlab="",ylab="",lwd=3)
legend("topright",legend=c("PAM","HCLUST"),col=c("blue","red"),pch=c(0,1),lwd=3)
dev.off()

clustering.pam.2 <- pam.2[["clustering"]]
write.table(clustering.pam.2,file="pam_2.txt", quote=FALSE)
###############################################################################
######################### Effect of drought at each day #######################
###############################################################################

###############################  Day 5 #######################################

#1. Generate metadata ----------------------------------------------------------
names_<-names(data_5)
probe <- c()
genotype <- c()
treatment <- c()
for (name in names_) {
  split <- strsplit(name, "\\.")[[1]]
  probe <- c(probe, substr(split[1], 2, nchar(split[1])))
  genotype <- c(genotype, ifelse(substr(split[2], 1, 1) == "M", "MM", "Penelli"))
  treatment <- c(treatment, ifelse(substr(split[2], 2, 2) == "D", "Drought", "Control"))
}
colData_5<-data.frame(Sample=probe,Treatment=treatment,Genotype=genotype)
rownames(colData_5)<-names(data_5)

table(genotype) #MM 16 Penelli 16
table(treatment) #Control 16 Drought 16

#Check whether the metadata is generated correctly
all(colnames(data_5) %in% rownames(colData_5))
all(colnames(data_5)==rownames(colData_5))

#2. Generate DESeq dataset -----------------------------------------------------
dds_5<-DESeqDataSetFromMatrix(countData=data_5, colData=colData_5,
                            design= ~Genotype*Treatment) # Effect of Genotype and 
                                                        #Treatment, as well as their
                                                        #interaction
#Pre-filtering: remove genes with less than 10 counts in total
keep<-rowSums(counts(dds_5))>=10
dds_5<-dds_5[keep,]

#Factor Level
dds_5$Genotype <- relevel(dds_5$Genotype, ref = "Penelli")
dds_5$Treatment <- relevel(dds_5$Treatment, ref = "Control")

#3. Perform Differential Analysis ----------------------------------------------
dds_5 <- DESeq(dds_5)

#4. Effect of drought -------------------------------------------------------- 
#Obtain results
resTreatment <- results(dds_5, contrast=c("Treatment","Drought","Control"))
#Select upregulated genes
indexes<- which(resTreatment$log2FoldChange > 0.5 & resTreatment$padj < 0.05)
upregulated<-rownames(resTreatment)[indexes]
write.table(upregulated,"upregulated_genes_drough_day5.txt"
            ,row.names = FALSE,quote = FALSE,col.names = FALSE)
length(upregulated)
#Select downregulated genes
indexes<-which(resTreatment$log2FoldChange < -0.5 & resTreatment$padj < 0.05)
downregulated<-rownames(resTreatment)[indexes]
length(downregulated)
write.table(downregulated,"downregulated_genes_drough_day5.txt"
            ,row.names = FALSE,quote = FALSE,col.names = FALSE)

#Volcano plot
resTreatment_df <- as.data.frame(resTreatment)
resTreatment_df$regulation <- "not significant"
resTreatment_df$regulation[rownames(resTreatment_df) %in% upregulated] <- "upregulated"
resTreatment_df$regulation[rownames(resTreatment_df) %in% downregulated] <- "downregulated"
pdf("Volcano_Plot_Day5_Drought_Effect.pdf")
ggplot(resTreatment_df, aes(x=log2FoldChange, y=-log10(padj), color=regulation)) +
  geom_point() +
  theme_minimal() +
  labs(title="Effect of Drought", x="log2 fold change", 
       y="-log10 adjusted p-value", color="Regulation") +
  scale_color_manual(values=c("not significant"="black",
                              "upregulated"="blue", 
                              "downregulated"="red"))
dev.off()

#Heatmap plot
rld_5 <- vst(dds_5, blind=FALSE)
sigGenes <- rownames(resTreatment)[which(abs(resTreatment$log2FoldChange) > 0.5 & resTreatment$padj < 0.05)]
df <- assay(rld_5)[sigGenes, ]
pdf("Heatmap_diff_genes.pdf",height =10 ,width=10)
pheatmap(df, scale = "row", clustering_distance_rows = "euclidean", 
         clustering_distance_cols = "euclidean",
         clustering_method = "complete", fontsize=5)
dev.off()

#GO term analysis



#5. Effect of genotype --------------------------------------------------------
resGenotype <- results(dds_5, contrast=c("Genotype","Penelli","MM"))
#Select upregulated genes
indexes<- which(resGenotype$log2FoldChange > 3 & resGenotype$padj < 0.05)
upregulated<-rownames(resGenotype)[indexes]
length(upregulated)
write.table(upregulated,"upregulated_genes_genotype_day5.txt"
            ,row.names = FALSE,quote = FALSE,col.names = FALSE)

#Select downregulated genes
indexes<-which(resGenotype$log2FoldChange < -3 & resGenotype$padj < 0.05)
downregulated<-rownames(resGenotype)[indexes]
length(downregulated)
write.table(downregulated,"downregulated_genes_genotype_day5.txt"
            ,row.names = FALSE,quote = FALSE,col.names = FALSE)

#Volcano plot
resGenotype_df <- as.data.frame(resGenotype)
resGenotype_df$regulation <- "not significant"
resGenotype_df$regulation[rownames(resGenotype_df) %in% upregulated] <- "upregulated"
resGenotype_df$regulation[rownames(resGenotype_df) %in% downregulated] <- "downregulated"
pdf("Volcano_Plot_Day5_Genotype_Effect.pdf")
ggplot(resGenotype_df, aes(x=log2FoldChange, y=-log10(padj), color=regulation)) +
  geom_point() +
  theme_minimal() +
  labs(title="Effect of Genotype", x="log2 fold change", 
       y="-log10 adjusted p-value", color="Regulation") +
  scale_color_manual(values=c("not significant"="black",
                              "upregulated"="blue", 
                              "downregulated"="red"))
dev.off()

#Heatmap plot
rld_5 <- vst(dds_5, blind=FALSE)
sigGenes <- rownames(resGenotype)[which(abs(resGenotype$log2FoldChange) > 3 & resTreatment$padj < 0.05)]
df <- assay(rld_5)[sigGenes, ]
pdf("Heatmap_diff_genes_genotypes.pdf",height =10 ,width=10)
pheatmap(df, scale = "row", clustering_distance_rows = "euclidean", 
         clustering_distance_cols = "euclidean",
         clustering_method = "complete", fontsize=5)
dev.off()

#GO term analysis

#6. Interaction drought:genotype ----------------------------------------------
resInteraction <- results(dds_5, name="GenotypeMM.TreatmentDrought")
#Select upregulated genes
indexes<- which(resInteraction$log2FoldChange > 1 & resInteraction$padj < 0.05)
upregulated<-rownames(resInteraction)[indexes]
length(upregulated)
write.table(upregulated,"upregulated_genes_interaction_day5.txt"
            ,row.names = FALSE,quote = FALSE,col.names = FALSE)

#Select downregulated genes
indexes<-which(resInteraction$log2FoldChange < -1 & resInteraction$padj < 0.05)
downregulated<-rownames(resInteraction)[indexes]
length(downregulated)
write.table(downregulated,"downregulated_genes_interaction_day5.txt"
            ,row.names = FALSE,quote = FALSE,col.names = FALSE)

#Volcano plot
resInteraction_df <- as.data.frame(resInteraction)
resInteraction_df$regulation <- "not significant"
resInteraction_df$regulation[rownames(resInteraction_df) %in% upregulated] <- "upregulated"
resInteraction_df$regulation[rownames(resInteraction_df) %in% downregulated] <- "downregulated"
pdf("Volcano_Plot_Day5_intetaction_Effect.pdf")
ggplot(resInteraction_df, aes(x=log2FoldChange, y=-log10(padj), color=regulation)) +
  geom_point() +
  theme_minimal() +
  labs(title="Effect of Drought and Genotype", x="log2 fold change", 
       y="-log10 adjusted p-value", color="Regulation") +
  scale_color_manual(values=c("not significant"="black",
                              "upregulated"="blue", 
                              "downregulated"="red"))
dev.off()

#Heatmap plot
rld_5 <- vst(dds_5, blind=FALSE)
sigGenes <- rownames(resInteraction)[which(abs(resInteraction$log2FoldChange) > 1 & resTreatment$padj < 0.05)]
df <- assay(rld_5)[sigGenes, ]
pdf("Heatmap_diff_genes_interaction.pdf",height =10 ,width=10)
pheatmap(df, scale = "row", clustering_distance_rows = "euclidean", 
         clustering_distance_cols = "euclidean",
         clustering_method = "complete", fontsize=5)
dev.off()

###############################  Day 8 #######################################

#1. Generate metadata ----------------------------------------------------------
names_<-names(data_8)
probe <- c()
genotype <- c()
treatment <- c()
for (name in names_) {
  split <- strsplit(name, "\\.")[[1]]
  probe <- c(probe, substr(split[1], 2, nchar(split[1])))
  genotype <- c(genotype, ifelse(substr(split[2], 1, 1) == "M", "MM", "Penelli"))
  treatment <- c(treatment, ifelse(substr(split[2], 2, 2) == "D", "Drought", "Control"))
}
colData_8<-data.frame(Sample=probe,Treatment=treatment,Genotype=genotype)
rownames(colData_8)<-names(data_8)

table(genotype) #MM 16 Penelli 16
table(treatment) #Control 16 Drought 16

#Check whether the metadata is generated correctly
all(colnames(data_8) %in% rownames(colData_8))
all(colnames(data_8)==rownames(colData_8))


#2. Generate DESeq dataset -----------------------------------------------------
dds_8<-DESeqDataSetFromMatrix(countData=data_8, colData=colData_8,
                              design= ~Genotype*Treatment) # Effect of Genotype and 
#Treatment, as well as their
#interaction
#Pre-filtering: remove genes with less than 10 counts in total
keep<-rowSums(counts(dds_8))>=10
dds_8<-dds_8[keep,]

#Factor Level
dds_8$Genotype <- relevel(dds_8$Genotype, ref = "Penelli")
dds_8$Treatment <- relevel(dds_8$Treatment, ref = "Control")

#3. Perform Differential Analysis ----------------------------------------------
dds_8 <- DESeq(dds_8)

#4. Effect of drought -------------------------------------------------------- 
#Obtain results
resTreatment <- results(dds_8, contrast=c("Treatment","Drought","Control"))
#Select upregulated genes
indexes<- which(resTreatment$log2FoldChange > 1 & resTreatment$padj < 0.05)
upregulated<-rownames(resTreatment)[indexes]
length(upregulated)
write.table(upregulated,"upregulated_genes_drought_day8.txt"
            ,row.names = FALSE,quote = FALSE,col.names = FALSE)

#Select downregulated genes
indexes<-which(resTreatment$log2FoldChange < -1 & resTreatment$padj < 0.05)
downregulated<-rownames(resTreatment)[indexes]
length(downregulated)
write.table(downregulated,"downregulated_genes_drought_day8.txt"
            ,row.names = FALSE,quote = FALSE,col.names = FALSE)


#Volcano plot
resTreatment_df <- as.data.frame(resTreatment)
resTreatment_df$regulation <- "not significant"
resTreatment_df$regulation[rownames(resTreatment_df) %in% upregulated] <- "upregulated"
resTreatment_df$regulation[rownames(resTreatment_df) %in% downregulated] <- "downregulated"
pdf("Volcano_Plot_Day8_Drought_Effect.pdf")
ggplot(resTreatment_df, aes(x=log2FoldChange, y=-log10(padj), color=regulation)) +
  geom_point() +
  theme_minimal() +
  labs(title="Effect of Drought", x="log2 fold change", 
       y="-log10 adjusted p-value", color="Regulation") +
  scale_color_manual(values=c("not significant"="black",
                              "upregulated"="blue", 
                              "downregulated"="red"))
dev.off()

#Heatmap plot
rld_8 <- vst(dds_8, blind=FALSE)
sigGenes <- rownames(resTreatment)[which(abs(resTreatment$log2FoldChange) > 1 & resTreatment$padj < 0.05)]
df <- assay(rld_8)[sigGenes, ]
pdf("Heatmap_diff_genes_day_8.pdf",height =10 ,width=10)
pheatmap(df, scale = "row", clustering_distance_rows = "euclidean", 
         clustering_distance_cols = "euclidean",
         clustering_method = "complete", fontsize=5)
dev.off()

#GO term analysis



#5. Effect of genotype --------------------------------------------------------
resGenotype <- results(dds_8, contrast=c("Genotype","Penelli","MM"))
#Select upregulated genes
indexes<- which(resGenotype$log2FoldChange > 3 & resGenotype$padj < 0.05)
upregulated<-rownames(resGenotype)[indexes]
length(upregulated)
write.table(upregulated,"upregulated_genes_genotype_day8.txt"
            ,row.names = FALSE,quote = FALSE,col.names = FALSE)

#Select downregulated genes
indexes<-which(resGenotype$log2FoldChange < -3 & resGenotype$padj < 0.05)
downregulated<-rownames(resGenotype)[indexes]
length(downregulated)
write.table(downregulated,"downregulated_genes_genotype_day8.txt"
            ,row.names = FALSE,quote = FALSE,col.names = FALSE)

#Volcano plot
resGenotype_df <- as.data.frame(resGenotype)
resGenotype_df$regulation <- "not significant"
resGenotype_df$regulation[rownames(resGenotype_df) %in% upregulated] <- "upregulated"
resGenotype_df$regulation[rownames(resGenotype_df) %in% downregulated] <- "downregulated"
pdf("Volcano_Plot_Day8_Genotype_Effect.pdf")
ggplot(resGenotype_df, aes(x=log2FoldChange, y=-log10(padj), color=regulation)) +
  geom_point() +
  theme_minimal() +
  labs(title="Effect of Genotype", x="log2 fold change", 
       y="-log10 adjusted p-value", color="Regulation") +
  scale_color_manual(values=c("not significant"="black",
                              "upregulated"="blue", 
                              "downregulated"="red"))
dev.off()

#Heatmap plot
rld_8 <- vst(dds_8, blind=FALSE)
sigGenes <- rownames(resGenotype)[which(abs(resGenotype$log2FoldChange) > 3 & resTreatment$padj < 0.05)]
df <- assay(rld_8)[sigGenes, ]
pdf("Heatmap_diff_genes_genotypes_day8.pdf",height =10 ,width=10)
pheatmap(df, scale = "row", clustering_distance_rows = "euclidean", 
         clustering_distance_cols = "euclidean",
         clustering_method = "complete", fontsize=5)
dev.off()

#GO term analysis

#6. Interaction drought:genotype ----------------------------------------------
resInteraction <- results(dds_8, name="GenotypeMM.TreatmentDrought")
#Select upregulated genes
indexes<- which(resInteraction$log2FoldChange > 1 & resInteraction$padj < 0.05)
upregulated<-rownames(resInteraction)[indexes]
length(upregulated)
write.table(upregulated,"upregulated_genes_interaction_day8.txt"
            ,row.names = FALSE,quote = FALSE,col.names = FALSE)

#Select downregulated genes
indexes<-which(resInteraction$log2FoldChange < -1 & resInteraction$padj < 0.05)
downregulated<-rownames(resInteraction)[indexes]
length(downregulated)
write.table(downregulated,"downregulated_genes_interaction_day8.txt"
            ,row.names = FALSE,quote = FALSE,col.names = FALSE)

#Volcano plot
resInteraction_df <- as.data.frame(resInteraction)
resInteraction_df$regulation <- "not significant"
resInteraction_df$regulation[rownames(resInteraction_df) %in% upregulated] <- "upregulated"
resInteraction_df$regulation[rownames(resInteraction_df) %in% downregulated] <- "downregulated"
pdf("Volcano_Plot_Day8_intetaction_Effect.pdf")
ggplot(resInteraction_df, aes(x=log2FoldChange, y=-log10(padj), color=regulation)) +
  geom_point() +
  theme_minimal() +
  labs(title="Effect of Drought and Genotype", x="log2 fold change", 
       y="-log10 adjusted p-value", color="Regulation") +
  scale_color_manual(values=c("not significant"="black",
                              "upregulated"="blue", 
                              "downregulated"="red"))
dev.off()

#Heatmap plot
rld_8 <- vst(dds_8, blind=FALSE)
sigGenes <- rownames(resInteraction)[which(abs(resInteraction$log2FoldChange) > 1 & resTreatment$padj < 0.05)]
df <- assay(rld_8)[sigGenes, ]
pdf("Heatmap_diff_genes_interaction_day8.pdf",height =10 ,width=10)
pheatmap(df, scale = "row", clustering_distance_rows = "euclidean", 
         clustering_distance_cols = "euclidean",
         clustering_method = "complete", fontsize=5)
dev.off()

#7. Gene coexpression network -------------------------------------------------

#Copy from cluster


###############################  Day 11 #######################################

#1. Generate metadata ----------------------------------------------------------
names_<-names(data_11)
probe <- c()
genotype <- c()
treatment <- c()
for (name in names_) {
  split <- strsplit(name, "\\.")[[1]]
  probe <- c(probe, substr(split[1], 2, nchar(split[1])))
  genotype <- c(genotype, ifelse(substr(split[2], 1, 1) == "M", "MM", "Penelli"))
  treatment <- c(treatment, ifelse(substr(split[2], 2, 2) == "D", "Drought", "Control"))
}
colData_11<-data.frame(Sample=probe,Treatment=treatment,Genotype=genotype)
rownames(colData_11)<-names(data_11)

table(genotype) #MM 15 Penelli 16
table(treatment) #Control 15 Drought 16

#Check whether the metadata is generated correctly
all(colnames(data_11) %in% rownames(colData_11))
all(colnames(data_11)==rownames(colData_11))


#2. Generate DESeq dataset -----------------------------------------------------
dds_11<-DESeqDataSetFromMatrix(countData=data_11, colData=colData_11,
                              design= ~Genotype*Treatment) # Effect of Genotype and 
#Treatment, as well as their
#interaction
#Pre-filtering: remove genes with less than 10 counts in total
keep<-rowSums(counts(dds_11))>=10
dds_11<-dds_11[keep,]

#Factor Level
dds_11$Genotype <- relevel(dds_11$Genotype, ref = "Penelli")
dds_11$Treatment <- relevel(dds_11$Treatment, ref = "Control")

#3. Perform Differential Analysis ----------------------------------------------
dds_11 <- DESeq(dds_11)

#4. Effect of drought -------------------------------------------------------- 
#Obtain results
resTreatment <- results(dds_11, contrast=c("Treatment","Drought","Control"))
#Select upregulated genes
indexes<- which(resTreatment$log2FoldChange > 2.5 & resTreatment$padj < 0.05)
upregulated<-rownames(resTreatment)[indexes]
length(upregulated)
write.table(upregulated,"upregulated_genes_drought_day11.txt"
            ,row.names = FALSE,quote = FALSE,col.names = FALSE)

#Select downregulated genes
indexes<-which(resTreatment$log2FoldChange < -2.5 & resTreatment$padj < 0.05)
downregulated<-rownames(resTreatment)[indexes]
length(downregulated)
write.table(downregulated,"downregulated_genes_drought_day11.txt"
            ,row.names = FALSE,quote = FALSE,col.names = FALSE)

#Volcano plot
resTreatment_df <- as.data.frame(resTreatment)
resTreatment_df$regulation <- "not significant"
resTreatment_df$regulation[rownames(resTreatment_df) %in% upregulated] <- "upregulated"
resTreatment_df$regulation[rownames(resTreatment_df) %in% downregulated] <- "downregulated"
pdf("Volcano_Plot_Day11_Drought_Effect.pdf")
ggplot(resTreatment_df, aes(x=log2FoldChange, y=-log10(padj), color=regulation)) +
  geom_point() +
  theme_minimal() +
  labs(title="Effect of Drought", x="log2 fold change", 
       y="-log10 adjusted p-value", color="Regulation") +
  scale_color_manual(values=c("not significant"="black",
                              "upregulated"="blue", 
                              "downregulated"="red"))
dev.off()

#Heatmap plot
rld_11 <- vst(dds_11, blind=FALSE)
sigGenes <- rownames(resTreatment)[which(abs(resTreatment$log2FoldChange) > 2.5 & resTreatment$padj < 0.05)]
df <- assay(rld_11)[sigGenes, ]
pdf("Heatmap_diff_genes_day_11.pdf",height =10 ,width=10)
pheatmap(df, scale = "row", clustering_distance_rows = "euclidean", 
         clustering_distance_cols = "euclidean",
         clustering_method = "complete", fontsize=5)
dev.off()

#GO term analysis



#5. Effect of genotype --------------------------------------------------------
resGenotype <- results(dds_11, contrast=c("Genotype","Penelli","MM"))
#Select upregulated genes
indexes<- which(resGenotype$log2FoldChange > 3.5 & resGenotype$padj < 0.05)
upregulated<-rownames(resGenotype)[indexes]
length(upregulated)
write.table(upregulated,"upregulated_genes_genotype_day11.txt"
            ,row.names = FALSE,quote = FALSE,col.names = FALSE)

#Select downregulated genes
indexes<-which(resGenotype$log2FoldChange < -3.5 & resGenotype$padj < 0.05)
downregulated<-rownames(resGenotype)[indexes]
length(downregulated)
write.table(downregulated,"downregulated_genes_genotype_day11.txt"
            ,row.names = FALSE,quote = FALSE,col.names = FALSE)

#Volcano plot
resGenotype_df <- as.data.frame(resGenotype)
resGenotype_df$regulation <- "not significant"
resGenotype_df$regulation[rownames(resGenotype_df) %in% upregulated] <- "upregulated"
resGenotype_df$regulation[rownames(resGenotype_df) %in% downregulated] <- "downregulated"
pdf("Volcano_Plot_Day11_Genotype_Effect.pdf")
ggplot(resGenotype_df, aes(x=log2FoldChange, y=-log10(padj), color=regulation)) +
  geom_point() +
  theme_minimal() +
  labs(title="Effect of Genotype", x="log2 fold change", 
       y="-log10 adjusted p-value", color="Regulation") +
  scale_color_manual(values=c("not significant"="black",
                              "upregulated"="blue", 
                              "downregulated"="red"))
dev.off()

#Heatmap plot
rld_11 <- vst(dds_11, blind=FALSE)
sigGenes <- rownames(resGenotype)[which(abs(resGenotype$log2FoldChange) > 3.5 & resTreatment$padj < 0.05)]
df <- assay(rld_8)[sigGenes, ]
pdf("Heatmap_diff_genes_genotypes_day11.pdf",height =10 ,width=10)
pheatmap(df, scale = "row", clustering_distance_rows = "euclidean", 
         clustering_distance_cols = "euclidean",
         clustering_method = "complete", fontsize=5)
dev.off()

#GO term analysis

#6. Interaction drought:genotype ----------------------------------------------
resInteraction <- results(dds_11, name="GenotypeMM.TreatmentDrought")
#Select upregulated genes
indexes<- which(resInteraction$log2FoldChange > 2.5 & resInteraction$padj < 0.05)
upregulated<-rownames(resInteraction)[indexes]
length(upregulated)
write.table(upregulated,"upregulated_genes_interaction_day11.txt"
            ,row.names = FALSE,quote = FALSE,col.names = FALSE)

#Select downregulated genes
indexes<-which(resInteraction$log2FoldChange < -2.5 & resInteraction$padj < 0.05)
downregulated<-rownames(resInteraction)[indexes]
length(downregulated)
write.table(downregulated,"downregulated_genes_interaction_day11.txt"
            ,row.names = FALSE,quote = FALSE,col.names = FALSE)

#Volcano plot
resInteraction_df <- as.data.frame(resInteraction)
resInteraction_df$regulation <- "not significant"
resInteraction_df$regulation[rownames(resInteraction_df) %in% upregulated] <- "upregulated"
resInteraction_df$regulation[rownames(resInteraction_df) %in% downregulated] <- "downregulated"
pdf("Volcano_Plot_Day11_intetaction_Effect.pdf")
ggplot(resInteraction_df, aes(x=log2FoldChange, y=-log10(padj), color=regulation)) +
  geom_point() +
  theme_minimal() +
  labs(title="Effect of Drought and Genotype", x="log2 fold change", 
       y="-log10 adjusted p-value", color="Regulation") +
  scale_color_manual(values=c("not significant"="black",
                              "upregulated"="blue", 
                              "downregulated"="red"))
dev.off()

#Heatmap plot
rld_11 <- vst(dds_11, blind=FALSE)
sigGenes <- rownames(resInteraction)[which(abs(resInteraction$log2FoldChange) > 2.5 & resTreatment$padj < 0.05)]
df <- assay(rld_11)[sigGenes, ]
pdf("Heatmap_diff_genes_interaction_day11.pdf",height =10 ,width=10)
pheatmap(df, scale = "row", clustering_distance_rows = "euclidean", 
         clustering_distance_cols = "euclidean",
         clustering_method = "complete", fontsize=5)
dev.off()

###############################################################################
################################# Effect of time ##############################
###############################################################################

############################# Split data by genotype ##########################

MM <- data.frame(matrix(ncol=0,nrow=nrow(exp_data)))
Penelli <- data.frame(matrix(ncol=0,nrow=nrow(exp_data)))
for (name in col_names) {
  gen <- substr(strsplit(name,"\\.")[[1]][2], 1, 1)
  print(gen)
  if (gen == "M") {
    MM <- cbind(MM, setNames(exp_data[, name, drop = FALSE], name))
  } else if (gen == "P") {
    Penelli <- cbind(Penelli, setNames(exp_data[, name, drop = FALSE], name))
  }
}


###############################  MM #######################################

#1. Generate metadata ----------------------------------------------------------
names_<-names(MM)
probe <- c()
treatment <- c()
time<-c()
for (name in names_) {
  split <- strsplit(name, "\\.")[[1]]
  probe <- c(probe, substr(split[1], 2, nchar(split[1])))
  treatment <- c(treatment, ifelse(substr(split[2], 2, 2) == "D", "Drought", "Control"))
  time<-c(time,as.numeric(substr(strsplit(name,"\\.")[[1]][3], 2, nchar(strsplit(name,"\\.")[[1]][3]))))
}
colData_MM<-data.frame(Sample=probe,Treatment=treatment,Time=time)
rownames(colData_MM)<-names(MM)

table(time) #D5:16; D8:16; D11;15
table(treatment) #Control 23 Drought 24

#Check whether the metadata is generated correctly
all(colnames(MM) %in% rownames(colData_MM))
all(colnames(MM)==rownames(colData_MM))


#2. Generate DESeq dataset -----------------------------------------------------
colData_MM$Time <- factor(colData_MM$Time) #Converting time to categorical variable
dds_MM<-DESeqDataSetFromMatrix(countData=MM, colData=colData_MM,
                               design= ~Treatment*Time) # Effect of time and effect of time
                                                        # in each treatment
#Pre-filtering: remove genes with less than 10 counts in total
keep<-rowSums(counts(dds_MM))>=10
dds_MM<-dds_MM[keep,]

#Factor Level
dds_MM$Time <- relevel(dds_MM$Time, ref = "5")
dds_MM$Treatment <- relevel(dds_MM$Treatment, ref = "Control")

#3. Perform Differential Analysis ----------------------------------------------
dds_MM <- DESeq(dds_MM)

#4. Effect of time and drougth 5-8 ----------------------------------------------------- 
#Obtain results
resultsNames(dds_MM)
resTime <- results(dds_MM, name="TreatmentDrought.Time8")
#Select upregulated genes
indexes<- which(resTime$log2FoldChange > 1 & resTime$padj < 0.05)
upregulated<-rownames(resTime)[indexes]
write.table(upregulated,"upregulated_genes_time_day5_8_MM.txt"
            ,row.names = FALSE,quote = FALSE,col.names = FALSE)
length(upregulated)
#Select downregulated genes
indexes<-which(resTime$log2FoldChange < -1 & resTime$padj < 0.05)
downregulated<-rownames(resTime)[indexes]
length(downregulated)
write.table(downregulated,"downregulated_genes__time_day5_8_MM.txt"
            ,row.names = FALSE,quote = FALSE,col.names = FALSE)

#Volcano plot
resTime_df <- as.data.frame(resTime)
resTime_df$regulation <- "not significant"
resTime_df$regulation[rownames(resTime_df) %in% upregulated] <- "upregulated"
resTime_df$regulation[rownames(resTime_df) %in% downregulated] <- "downregulated"
pdf("Volcano_Plot_time_day5_8_MM.pdf")
ggplot(resTime_df, aes(x=log2FoldChange, y=-log10(padj), color=regulation)) +
  geom_point() +
  theme_minimal() +
  labs(title="Effect of Drought Days 5-8", x="log2 fold change", 
       y="-log10 adjusted p-value", color="Regulation") +
  scale_color_manual(values=c("not significant"="black",
                              "upregulated"="blue", 
                              "downregulated"="red"))
dev.off()

#Heatmap plot
rld_MM <- vst(dds_MM, blind=FALSE)
sigGenes <- rownames(resTime)[which(abs(resTime$log2FoldChange) > 1 & resTime$padj < 0.05)]
df <- assay(rld_MM)[sigGenes, ]
pdf("Heatmap_diff_time_day5_8_MM.pdf",height =10 ,width=10)
pheatmap(df, scale = "row", clustering_distance_rows = "euclidean", 
         clustering_distance_cols = "euclidean",
         clustering_method = "complete", fontsize=5)
dev.off()

#GO term analysis

ds_MM <- DESeq(dds_MM)

#4. Effect of time and drougth 5-11 ----------------------------------------------------- 
#Obtain results
resultsNames(dds_MM)
resTime <- results(dds_MM, name="TreatmentDrought.Time11")
#Select upregulated genes
indexes<- which(resTime$log2FoldChange > 2.5 & resTime$padj < 0.05)
upregulated<-rownames(resTime)[indexes]
write.table(upregulated,"upregulated_genes_time_day5_11_MM.txt"
            ,row.names = FALSE,quote = FALSE,col.names = FALSE)
length(upregulated)
#Select downregulated genes
indexes<-which(resTime$log2FoldChange < -2.5 & resTime$padj < 0.05)
downregulated<-rownames(resTime)[indexes]
length(downregulated)
write.table(downregulated,"downregulated_genes__time_day5_11_MM.txt"
            ,row.names = FALSE,quote = FALSE,col.names = FALSE)

#Volcano plot
resTime_df <- as.data.frame(resTime)
resTime_df$regulation <- "not significant"
resTime_df$regulation[rownames(resTime_df) %in% upregulated] <- "upregulated"
resTime_df$regulation[rownames(resTime_df) %in% downregulated] <- "downregulated"
pdf("Volcano_Plot_time_day5_11.pdf")
ggplot(resTime_df, aes(x=log2FoldChange, y=-log10(padj), color=regulation)) +
  geom_point() +
  theme_minimal() +
  labs(title="Effect of Drought Days 5-11_MM", x="log2 fold change", 
       y="-log10 adjusted p-value", color="Regulation") +
  scale_color_manual(values=c("not significant"="black",
                              "upregulated"="blue", 
                              "downregulated"="red"))
dev.off()

#Heatmap plot
rld_MM <- vst(dds_MM, blind=FALSE)
sigGenes <- rownames(resTime)[which(abs(resTime$log2FoldChange) > 2.5 & resTime$padj < 0.05)]
df <- assay(rld_MM)[sigGenes, ]
pdf("Heatmap_diff_time_day5_11_MM.pdf",height =10 ,width=10)
pheatmap(df, scale = "row", clustering_distance_rows = "euclidean", 
         clustering_distance_cols = "euclidean",
         clustering_method = "complete", fontsize=5)
dev.off()

#GO term analysis



###############################  Penelli #######################################

#1. Generate metadata ----------------------------------------------------------
names_<-names(Penelli)
probe <- c()
treatment <- c()
time<-c()
for (name in names_) {
  split <- strsplit(name, "\\.")[[1]]
  probe <- c(probe, substr(split[1], 2, nchar(split[1])))
  treatment <- c(treatment, ifelse(substr(split[2], 2, 2) == "D", "Drought", "Control"))
  time<-c(time,as.numeric(substr(strsplit(name,"\\.")[[1]][3], 2, nchar(strsplit(name,"\\.")[[1]][3]))))
}
colData_Penelli<-data.frame(Sample=probe,Treatment=treatment,Time=time)
rownames(colData_Penelli)<-names(Penelli)

table(time) #D5:16; D8:16; D11;16
table(treatment) #Control 24 Drought 24

#Check whether the metadata is generated correctly
all(colnames(Penelli) %in% rownames(colData_Penelli))
all(colnames(Penelli)==rownames(colData_Penelli))


#2. Generate DESeq dataset -----------------------------------------------------
colData_Penelli$Time <- factor(colData_Penelli$Time) #Converting time to categorical variable
dds_Penelli<-DESeqDataSetFromMatrix(countData=Penelli, colData=colData_Penelli,
                               design= ~Treatment*Time) # Effect of time and effect of time
# in each treatment
#Pre-filtering: remove genes with less than 10 counts in total
keep<-rowSums(counts(dds_Penelli))>=10
dds_Penelli<-dds_Penelli[keep,]

#Factor Level
dds_Penelli$Time <- relevel(dds_Penelli$Time, ref = "5")
dds_Penelli$Treatment <- relevel(dds_Penelli$Treatment, ref = "Control")

#3. Perform Differential Analysis ----------------------------------------------
dds_Penelli <- DESeq(dds_Penelli)

#4. Effect of time and drougth 5-8 ----------------------------------------------------- 
#Obtain results
resultsNames(dds_Penelli)
resTime <- results(dds_Penelli, name="TreatmentDrought.Time8")
#Select upregulated genes
indexes<- which(resTime$log2FoldChange > 0.5 & resTime$padj < 0.05)
upregulated<-rownames(resTime)[indexes]
write.table(upregulated,"upregulated_genes_time_day5_8_Penelli.txt"
            ,row.names = FALSE,quote = FALSE,col.names = FALSE)
length(upregulated)
#Select downregulated genes
indexes<-which(resTime$log2FoldChange < -0.5 & resTime$padj < 0.05)
downregulated<-rownames(resTime)[indexes]
length(downregulated)
write.table(downregulated,"downregulated_genes__time_day5_8_Penelli.txt"
            ,row.names = FALSE,quote = FALSE,col.names = FALSE)

#Volcano plot
resTime_df <- as.data.frame(resTime)
resTime_df$regulation <- "not significant"
resTime_df$regulation[rownames(resTime_df) %in% upregulated] <- "upregulated"
resTime_df$regulation[rownames(resTime_df) %in% downregulated] <- "downregulated"
pdf("Volcano_Plot_time_day5_8_Penelli.pdf")
ggplot(resTime_df, aes(x=log2FoldChange, y=-log10(padj), color=regulation)) +
  geom_point() +
  theme_minimal() +
  labs(title="Effect of Drought Days 5-8", x="log2 fold change", 
       y="-log10 adjusted p-value", color="Regulation") +
  scale_color_manual(values=c("not significant"="black",
                              "upregulated"="blue", 
                              "downregulated"="red"))
dev.off()

#Heatmap plot
rld_Penelli <- vst(dds_Penelli, blind=FALSE)
sigGenes <- rownames(resTime)[which(abs(resTime$log2FoldChange) > 0.5 & resTime$padj < 0.05)]
df <- assay(rld_Penelli)[sigGenes, ]
pdf("Heatmap_diff_time_day5_8_Penelli.pdf",height =10 ,width=10)
pheatmap(df, scale = "row", clustering_distance_rows = "euclidean", 
         clustering_distance_cols = "euclidean",
         clustering_method = "complete", fontsize=5)
dev.off()

#GO term analysis

#4. Effect of time and drougth 5-11 ----------------------------------------------------- 
#Obtain results
resultsNames(dds_Penelli)
resTime <- results(dds_Penelli, name="TreatmentDrought.Time11")
#Select upregulated genes
indexes<- which(resTime$log2FoldChange > 2.5 & resTime$padj < 0.05)
upregulated<-rownames(resTime)[indexes]
write.table(upregulated,"upregulated_genes_time_day5_11_Penelli.txt"
            ,row.names = FALSE,quote = FALSE,col.names = FALSE)
length(upregulated)
#Select downregulated genes
indexes<-which(resTime$log2FoldChange < -2.5 & resTime$padj < 0.05)
downregulated<-rownames(resTime)[indexes]
length(downregulated)
write.table(downregulated,"downregulated_genes__time_day5_11_Penelli.txt"
            ,row.names = FALSE,quote = FALSE,col.names = FALSE)

#Volcano plot
resTime_df <- as.data.frame(resTime)
resTime_df$regulation <- "not significant"
resTime_df$regulation[rownames(resTime_df) %in% upregulated] <- "upregulated"
resTime_df$regulation[rownames(resTime_df) %in% downregulated] <- "downregulated"
pdf("Volcano_Plot_time_day5_11_Penelli.pdf")
ggplot(resTime_df, aes(x=log2FoldChange, y=-log10(padj), color=regulation)) +
  geom_point() +
  theme_minimal() +
  labs(title="Effect of Drought Days 5-11", x="log2 fold change", 
       y="-log10 adjusted p-value", color="Regulation") +
  scale_color_manual(values=c("not significant"="black",
                              "upregulated"="blue", 
                              "downregulated"="red"))
dev.off()

#Heatmap plot
rld_Penelli <- vst(dds_Penelli, blind=FALSE)
sigGenes <- rownames(resTime)[which(abs(resTime$log2FoldChange) > 2.5 & resTime$padj < 0.05)]
df <- assay(rld_Penelli)[sigGenes, ]
pdf("Heatmap_diff_time_day5_11_Penelli.pdf",height =10 ,width=10)
pheatmap(df, scale = "row", clustering_distance_rows = "euclidean", 
         clustering_distance_cols = "euclidean",
         clustering_method = "complete", fontsize=5)
dev.off()

#GO term analysis
