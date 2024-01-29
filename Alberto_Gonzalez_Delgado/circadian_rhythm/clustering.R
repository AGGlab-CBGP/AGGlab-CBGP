#!/bin/R

#Alberto Gonzalez Delgado
#01/2024
#Centro de Biotecnologia y Genomica de Plantas (UPM-INIA/CSIC)

#setwd("C:/Users/admin/Desktop/projects") #Home-desktop
setwd("D:/Alberto/projects") #Lab-desktop

# Import packages
library(rgl)
library(DESeq2)
library(WGCNA)
library(sva)
library(dplyr)
library(ggplot2)
library(gridExtra)
library(stringr)

set.seed(123) # For reproducibility
#1. Read data --------------------------------------------------------------------------------------------------
data<-read.delim("documents/02.1.raw_counts.tsv")
oscill_data<-read.delim('./oscillation_descriptors.tsv')
mean_counts<-read.delim('./0.1.Circadian_Rhythm/genot_mean/vst_counts_mean_genotype.tsv',sep=',')
gene_list<-oscill_data$CycID
head(data)

rownames(data)<-data$gene_id
data<-data[c(8:ncol(data))]

#Generate metadata
titles=colnames(data[1:ncol(data)])
genotypes<-c(rep("EID1",24),rep("LE",24),rep("LNK2",24),rep("MM",24))
genotypes<-c(rep(genotypes,3))
photoperiods<-c(rep("LD",96),rep("ND",96),rep("SD",96))
time<-rep(seq(1, 23, by = 2), each = 2, times = 12)
replicates<-rep(c("r1","r2"),144)
colData<-data.frame(Sample=titles,photoperiod=photoperiods,genotype=genotypes,replicate=replicates)
rownames(colData)<-colData$Sample


#Define conditions for trait: photoperiod
photoperiods<-unique(colData[,2])
rownames(colData)<-colData$Sample
all(rownames(colData) %in% colnames(data)) #Are all the conditions? #(data.subset)
all(rownames(colData) == colnames(data)) #They must to be in the same order

##2. Remove batch effects ----------------------------------------------------------------------------------
batch<-c(rep(0,96),rep(1,96),rep(0,96))
colData$batch<-batch
batch = colData$batch
combat_edata = ComBat_seq(as.matrix(data), batch=batch)

#3. Normalization --------------------------------------------------------------------
#Generate DESeq2 dataset
dds <- DESeqDataSetFromMatrix(countData = combat_edata,
                              colData=colData,
                              design = ~1) # Not specifying model because we need  
#the deseq data set to performe variance stabilizing transformation

#5. Preliminary exploratory data analysis ------------------------------------------------------------------
#Detect outliers samples and genes
filtered<-data[rownames(data) %in% gene_list,]
gsg <- goodSamplesGenes(t(filtered))
summary(gsg)
gsg$allOK
table(gsg$goodGenes)
table(gsg$goodSamples)

dds_data<-dds[rownames(dds) %in% gene_list,]
#Perform Variance Stabilizing Transformation
dds_norm<-vst(dds_data)
norm.counts_<-as.data.frame(assay(dds_norm))

#7. Exploratory data analysis -------------------------------------------------------
  
  # Hierarchical clustering 
htree<-hclust(dist(t(norm.counts_)),method='average') #Data has to be transposed (genes as columns, samples as rows)
pdf("hierarchical_clustering.pdf",width=15,height=8)
plot(htree,cex=0.4)  
dev.off()

#PCA
pca<- prcomp(t(norm.counts_)) #Data has to be transposed (genes as columns, samples as rows)
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

#6. Perform network construction, module eigengene calculation, module-trait correlation --------------------

#Choose a set of soft-thresholding powers
traits<-unique(colData[,2])
datExpr<-as.data.frame(norm.counts_)

power<-c(c(1:11),seq(12,50,by=2))

# Call the network topology analysis function
sft<-pickSoftThreshold(as.data.frame(datExpr,
                       powerVector = power,
                       networkType = "signed",
                       verbose=5)
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
print(grid.arrange(p1,p2,ncol=1)) #Power=1

powers<-c(1,1,26)
#9. Construct networks per photoperiod
datExpr<-as.data.frame(norm.counts_)
sampleTable<-colData
traits<-unique(sampleTable[,2])
clusters<-matrix(nrow = nrow(datExpr))
rownames(clusters)<-rownames(datExpr)
clusters<-as.data.frame(clusters) %>% select(-V1)
networks<-list()
count<-1
temp_cor <- cor
cor <- WGCNA::cor
for(trait in traits){
  my_net = blockwiseModules(t(datExpr %>% select(matches(paste0(trait,".*"))))
                            , alphaLevel = 0.05, plot = T,
                            networkType = "signed", power = 1,
                            minModuleSize = 40, maxBlockSize = 25000,
                            reassignThreshold = 0, minKMEtoStay = 0.65,
                            mergeCutHeight = 0.10, numericLabels = TRUE,
                            pamRespectsDendro = FALSE,verbose=3)
  dynamicLabels=paste(paste0(trait), "_", str_pad(my_net$colors, 3, pad="0"), sep="")
  
  names(dynamicLabels)<-names(my_net$colors)
  temp<-as.data.frame(dynamicLabels)
  names(temp)<-paste0(trait)
  clusters <- merge(clusters, temp, by = "row.names", all = TRUE)
  rownames(clusters) <- clusters$Row.names
  clusters <- clusters%>% select(-Row.names)
  networks[[trait]]<-my_net
  count<-count+1
}
cor <- temp_cor
library(tibble)
rownames(mean_counts)<-mean_counts$gene_id
pdf("vst_counts_modules.pdf")
#arreglalo
for (photoperiod in photoperiods) {
  mod_data_list <- list()
  subset <- mean_counts %>% dplyr::select(matches(photoperiod,'.*'))
  modules <- unique(unname(unlist(clusters[photoperiod])))
  for(i in 1:(length(modules)-1)){
    module <- modules[i]
    genes<-rownames(clusters[clusters[photoperiod]==module,][photoperiod])
    mod_data <- subset[rownames(subset) %in% genes,]
    mod_data$gene_id<-rownames(mod_data)
    mod_data <- mod_data %>%
      pivot_longer(cols =-gene_id, names_to = "sample", values_to = "counts")%>%
      separate(sample, into = c("photoperiod", "CT", "replicate"), sep = "_")%>%
      dplyr::filter(photoperiod==photoperiod)
    
 
    mod_data_long <- mod_data %>% 
      mutate(time = as.numeric(str_remove(CT, "CT")))
    
    
    mod_data_long$module <- i-1
    
    mod_data_list[[paste(module)]] <- mod_data_long
    all_mod_data <- bind_rows(mod_data_list)
    all_mod_data$module_group <- as.integer(as.numeric(as.factor(all_mod_data$module)) - 1) %/% 1 + 1
  }
  p<-ggplot(all_mod_data, aes(x = time, y = counts, group =gene_id)) +
    geom_line(color = "lightgrey", alpha = 0.25) +
    stat_summary(aes(group = module), fun = mean, geom = "line", color = "black", size = 1) +
    facet_wrap(~module) +
    labs(x = "Day", y = "vst", title = paste0(photoperiod)) +
    theme_bw()
  print(p)
}
dev.off()

palette_ <- c("#000000", "#FFFF00", "#1CE6FF", "#FF34FF", "#FF4A46", "#008941", "darkgreen", "#A30059",
              "darkorange", "red", "#0020A6", "#63FFAC", "#B79762", "#004D43", "#8FB0FF", "#997D87",
              "#5A0087", "#809693", "orange", "#1B4400", "#4FC601", "#3B5DFF", "#4A3B53", "#FF2F90")


distr<-list()
for(photoperiod in photoperiods){
  modules <- unique(unname(unlist(clusters[photoperiod])))
  dist_mod<-list()
  for(i in 1:(length(modules)-1)){
    print(modules[i])
    module <- modules[i]
    gene_list<-rownames(clusters[clusters[photoperiod]==module,][photoperiod])
    filtered_data<- oscill_data %>% filter(CycID %in% gene_list) %>% select(matches(".*phase")) %>% select(matches(paste0(photoperiod)))
    dist_mod[[i]]<-filtered_data %>%unlist() %>% unname()
  }
  distr[[photoperiod]]<-dist_mod
}


df_phase_LD <- bind_rows(lapply(distr$LD, as.data.frame), .id = "module")

df_phase_ND <- bind_rows(lapply(distr$ND, as.data.frame), .id = "module") 

df_phase_SD <- bind_rows(lapply(distr$SD, as.data.frame), .id = "module") 

# Order modules by photoperiod
df_phase_ND$photoperiod <- "ND"
df_phase_SD$photoperiod <- "SD"
df_phase_LD$photoperiod <- "LD"
df_phase <- bind_rows(df_phase_ND, df_phase_SD, df_phase_LD)
df_phase$module <- factor(df_phase$module, levels = unique(df_phase$module))
names(df_phase)[2]<-"value"
mean_phase_by_module <- df_phase %>%
  group_by(module,photoperiod) %>%
  dplyr::summarize(mean_value = median(value)) %>%
  dplyr::arrange(mean_value)
print(paste(mean_phase_by_module[mean_phase_by_module$photoperiod=='LD',]$module,collapse=","))
print(paste(mean_phase_by_module[mean_phase_by_module$photoperiod=='ND',]$module,collapse=","))
print(paste(mean_phase_by_module[mean_phase_by_module$photoperiod=='SD',]$module,collapse=","))
# Plot
order_LD <- c(3,9,7,4,6,2,5,1,8)
order_ND <- c(7,9,4,6,2,5,8,1,3)
order_SD <- c(4,7,6,2,5,1,8,3)
LD<-df_phase[df_phase$photoperiod=='LD',]
LD <- merge(LD, mean_phase_by_module, by = c("module", "photoperiod"))
LD <- LD %>%
  mutate(
    adjusted_value = ifelse(value > mean_value + 12, value - 24,
                            ifelse(value < mean_value -12, value + 24, value))
  )

ND<-df_phase[df_phase$photoperiod=='ND',]
ND <- merge(ND, mean_phase_by_module, by = c("module", "photoperiod"))
ND <- ND %>%
  mutate(
    adjusted_value = ifelse(value > mean_value + 12, value - 24,
                            ifelse(value < mean_value -12, value + 24, value))
  )
SD<-df_phase[df_phase$photoperiod=='SD',]
SD <- merge(SD, mean_phase_by_module, by = c("module", "photoperiod"))
SD <- SD %>%
  mutate(
    adjusted_value = ifelse(value > mean_value + 12, value - 24,
                            ifelse(value < mean_value -12, value + 24, value))
  )
LD$module<-factor(LD$module, levels=order_LD)
ND$module<-factor(ND$module, levels=order_ND)
SD$module<-factor(SD$module, levels=order_SD)

library(ggsignif)
pdf("phase_distribution.pdf",width = 5,height = 4)
p1<-ggplot(LD, aes(x = module, y = adjusted_value)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(alpha = 0.1,size=0.1) +
  theme_minimal() +
  labs(x = "Module", y = "Phase (h)",title="LD") +
  theme(axis.text.x = element_text( hjust = 1))+
  geom_signif(comparisons = list(c(1,2),c(2,3),c(3,4),c(4,5),
                                 c(5,6),c(6,7),c(7,8),c(8,9)),map_signif_level = TRUE)+
  scale_y_continuous(breaks = seq(0, 25, by = 5)) +
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = 18, ymax = Inf, fill = "darkgrey", alpha = 0.5)
p2<-ggplot(ND, aes(x = module, y = adjusted_value)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(alpha = 0.1,size=0.1) +
  theme_minimal() +
  labs(x = "Module", y = "Phase (h)",title="ND") +
  theme(axis.text.x = element_text( hjust = 1))+
  geom_signif(comparisons = list(c(1,2),c(2,3),c(3,4),c(4,5),
                                 c(5,6),c(6,7),c(7,8),c(8,9)),map_signif_level = TRUE)+
  scale_y_continuous(breaks = seq(0, 25, by = 5))+
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = 12, ymax = Inf, fill = "darkgrey", alpha = 0.5)
p3<-ggplot(SD, aes(x = module, y = adjusted_value)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(alpha = 0.1,size=0.1) +
  theme_minimal() +
  labs(x = "Module", y = "Phase (h)",title="SD") +
  theme(axis.text.x = element_text( hjust = 1))+
  geom_signif(comparisons = list(c(1,2),c(2,3),c(3,4),c(4,5),
                                 c(5,6),c(6,7),c(7,8)),map_signif_level = TRUE)+
  scale_y_continuous(breaks = seq(0, 25, by = 5))+
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = 6, ymax = Inf, fill = "darkgrey", alpha = 0.5)
print(p1)
print(p2)
print(p3)
dev.off()


#Sankey plot -------------------------------------------------------------------------------------
network_info<-data.frame(gene=character(),module=character(),photoperiod=character())
for(photoperiod in photoperiods){
  modules <- unique(unname(unlist(clusters[photoperiod])))
  dist_mod<-list()
  for(i in 1:(length(modules)-1)){
    module <- modules[i]
    gene_list<-rownames(clusters[clusters[photoperiod]==module,][photoperiod])
    temp_df<-data.frame(gene=gene_list,module=rep(module,length(gene_list)),photoperiod=rep(photoperiod,length(gene_list)))
    network_info<-rows_append(network_info,temp_df)
  }
}
grouped<-network_info%>%group_by(gene)
sankey_info<- data.frame(gene=character(),
                         node=character(),x=numeric(),
                         next_node=character(), next_x=numeric())
for (name in unique(network_info$gene)) {
  group <- filter(grouped, gene == name)
  if (nrow(group) > 1) {
    for (i in 1:(nrow(group)-1)) {
      new_row <- data.frame(gene=name, node=group$module[i+1], x=group$photoperiod[i+1],
                            next_node=group$module[i], next_x=group$photoperiod[i])
      sankey_info <- rbind(sankey_info, new_row)
    }
  }
}
labels<-unique(network_info$module)
colors<-c("#000000", "#FFFF00", "#1CE6FF", "#FF34FF", "#3B5DFF", "#008941", "darkgreen", "#A30059",
          "darkorange", "red", "#0020A6", "#63FFAC", "#B79762", "#004D43", "#8FB0FF", "#997D87",
          "#5A0087",#LD
          "#000000", "#FFFF00", "#1CE6FF", "#FF34FF", "#3B5DFF", "#008941", "darkgreen", "#A30059",
          "darkorange", "red", "#0020A6", "#63FFAC", "#B79762", "#004D43", #ND
          "#000000", "#FFFF00", "#1CE6FF", "#FF34FF", "#3B5DFF", "#008941", "darkgreen", "#A30059",
          "darkorange", "red","#0020A6")
label_dict <- setNames(0:length(labels), labels)
sankey_info$node <- label_dict[sankey_info$node]
sankey_info$next_node <- label_dict[sankey_info$next_node]
link <- aggregate(x ~ node + next_node, sankey_info, length)
names(link) <- c("source","target","value")
link<-as.list(link)
#install.packages("plotly")
library(plotly)
p<-plot_ly(type="sankey",orientation="h",
           node=list(
             label=labels,
             color=colors,
             pad=15,
             thickness=20,
             line=list(color="black",width=0.7)
           ),
           link=link
)
htmlwidgets::saveWidget(p, "sankey.html")
###################################################################################################
#11.4 Distribution of descriptors
descriptor="amplitude"
library(multcompView)
plot_descriptors<-function(descriptor) { 
  distr<-list()
  for(photoperiod in photoperiods){
    modules <- unique(unname(unlist(clusters[photoperiod])))
    dist_mod<-list()
    for(i in 1:(length(modules)-1)){
      module <- modules[i]
      gene_list_<-rownames(clusters[clusters[photoperiod]==module,][photoperiod])
      filtered_data<- oscill_data %>% filter(CycID %in% gene_list_) %>% select(matches(paste0(".*",descriptor))) %>% select(matches(paste0(photoperiod))) %>% unlist() %>% unname()
      dist_mod[[as.numeric(i)]]<-filtered_data
    }
    distr[[photoperiod]]<-dist_mod
  }
  
  df_LD <- bind_rows(lapply(distr$LD, as.data.frame), .id = "module")
  
  df_ND <- bind_rows(lapply(distr$ND, as.data.frame), .id = "module") 
  
  df_SD <- bind_rows(lapply(distr$SD, as.data.frame), .id = "module") 
  
  # Order modules by photoperiod
  df_ND$photoperiod <- "ND"
  df_SD$photoperiod <- "SD"
  df_LD$photoperiod <- "LD"
  df <- bind_rows(df_ND, df_SD, df_LD)
  df$module <- factor(df$module, levels = unique(df$module))
  # Plot
  order_LD <- c(3,9,7,4,6,10,2,5,1,8)
  order_ND <- c(7,9,4,6,2,10,5,8,1,3)
  order_SD <- c(4,7,6,2,9,5,1,8,3)
  LD<-df[df$photoperiod=='LD',]
  ND<-df[df$photoperiod=='ND',]
  SD<-df[df$photoperiod=='SD',]
  LD$module<-factor(LD$module, levels=order_LD)
  ND$module<-factor(ND$module, levels=order_ND)
  SD$module<-factor(SD$module, levels=order_SD)
  names(LD)[2]<-"value"
  names(ND)[2]<-"value"
  names(SD)[2]<-"value"
  # Multianova
  anova <- aov(value ~ module, data = LD)
  tukey<-TukeyHSD(anova)
  cld <- multcompLetters4(anova, tukey)
  LD_tuk<- LD %>% group_by(module) %>%
    summarise(mean=mean(value), quant = quantile(value, probs = 0.75)) %>%
    arrange(desc(mean))
  cld <- as.data.frame.list(cld$module)
  LD_tuk$cld <- cld$Letters
  
  anova <- aov(value ~ module, data = ND)
  tukey<-TukeyHSD(anova)
  cld <- multcompLetters4(anova, tukey)
  ND_tuk<- ND %>% group_by(module) %>%
    summarise(mean=mean(value), quant = quantile(value, probs = 0.75)) %>%
    arrange(desc(mean))
  cld <- as.data.frame.list(cld$module)
  ND_tuk$cld <- cld$Letters
  
  anova <- aov(value ~ module, data = SD)
  tukey<-TukeyHSD(anova)
  cld <- multcompLetters4(anova, tukey)
  
  SD_tuk<- SD %>% group_by(module) %>%
    summarise(mean=mean(value), quant = quantile(value, probs = 0.75)) %>%
    arrange(desc(mean))
  cld <- as.data.frame.list(cld$module)
  SD_tuk$cld <- cld$Letters
  
  
  p1<-ggplot(LD, aes(x = module, y = value)) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(alpha = 0.1,size=0.1) +
    theme_minimal() +
    labs(x = "Module", y = paste0(descriptor),title="LD") +
    theme(axis.text.x = element_text( hjust = 1))+
    annotate("rect", xmin = which(order_LD==5), xmax = length(order_LD)+0.5, ymin = -Inf, ymax = Inf,
             fill = "darkgrey", alpha = 0.25)  +
    geom_text(data = LD_tuk, aes(x = module, y = quant, label = cld))
  
  p2<-ggplot(ND, aes(x = module, y = value)) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(alpha = 0.1,size=0.1) +
    theme_minimal() +
    labs(x = "Module", y = paste0(descriptor),title="ND") +
    theme(axis.text.x = element_text( hjust = 1))+
    annotate("rect", xmin = which(order_ND==2), xmax = length(order_ND)+0.5, ymin = -Inf, ymax = Inf,
             fill = "darkgrey", alpha = 0.25)   +
    geom_text(data = ND_tuk, aes(x = module, y = quant, label = cld))
  
  p3<-ggplot(SD, aes(x = module, y = value)) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(alpha = 0.1,size=0.1) +
    theme_minimal() +
    labs(x = "Module", y = paste0(descriptor),title="SD") +
    theme(axis.text.x = element_text( hjust = 1))+
    annotate("rect", xmin = which(order_SD==7), xmax = length(order_SD)+0.5, ymin = -Inf, ymax = Inf,
             fill = "darkgrey", alpha = 0.25)   +
    geom_text(data = SD_tuk, aes(x = module, y = quant, label = cld))
  
  pdf(paste0("0.1.Circadian_Rhythm/distribution/",descriptor,"_distribution.pdf"),width = 5,height = 4)
  print(p1)
  print(p2)
  print(p3)
  dev.off()
}

descriptors<-c("BH.Q","period","amplitude","AUC_shape","AUC_exp","width","time_increasing",
               "time_decreasing","skewness","kurtosis","prominence","height","residuals",
               "start","end")
plot_descriptors("amplitude")
for(descriptor in descriptors){
  plot_descriptors(descriptor)
}

#################################################################################################################################################
#1. GO terms enrichment -----------------------------------------------------------------------------------------------------------------
library(ComplexHeatmap)
library(topGO)
#2.Generate GO annotation -----------------------------------------------------
geneID2GO<-readMappings(file = "GO_annotation.txt")
#GO2geneID <- inverseList(geneID2GO)
geneNames <- names(geneID2GO)


#Function: GO enrichment from gene list -----------------------------------------------
GO_analysis<-function(gene_list,ntop,type,output){
  pdf(paste0(output,".pdf"),width=12)
  regulated_genes <- factor(as.integer(geneNames %in% gene_list))
  names(regulated_genes) <- geneNames
  GOdata <- new("topGOdata", ontology = paste0(type),
                allGenes = regulated_genes,
                annot = annFUN.gene2GO,
                #nodeSize = 5, #Prunning (5-10)
                gene2GO = geneID2GO)
  resultFis <- runTest(GOdata, algorithm = "classic", statistic = "fisher")
  pvalFis <- score(resultFis)
  resultKS <- runTest(GOdata, algorithm = "classic", statistic = "ks")
  test.stat <- new("weightCount", testStatistic = GOFisherTest, name = "Fisher test")
  resultWeight <- getSigGroups(GOdata, test.stat)
  goEnrichment <- GenTable(GOdata, classic = resultFis, KS = resultKS, weight = resultWeight,
                           orderBy = "classic", ranksOf = "classic", topNodes = 100,numChar=99)
  write.table(goEnrichment,paste0(output,".csv"))
  showSigOfNodes(GOdata, score(resultFis), useInfo = 'all',firstSigNodes = 5)
  goEnrichment$pval <- as.numeric(goEnrichment$classic)
  goEnrichment$size<-as.numeric(goEnrichment$Significant)
  goEnrichment <- goEnrichment[goEnrichment$pval < 0.05,] 
  
  goEnrichment <- goEnrichment[,c("GO.ID","Term","pval","size")]
  ggdata <- goEnrichment[1:ntop,]
  ggdata$Term <- factor(ggdata$Term, levels = rev(ggdata$Term)) 
  gg1 <- ggplot(ggdata,
                aes(x = Term, y = -log10(pval), size = size, fill = -log10(pval))) +
    
    expand_limits(y = 1) +
    geom_point(shape = 21) +
    scale_size(range = c(2.5,12.5)) +
    scale_fill_continuous(low = 'royalblue', high = 'red4') +
    
    xlab('') + ylab('Enrichment score') +
    labs(
      title = paste0('GO',type),
      subtitle = paste0("Top ", ntop," terms ordered by t Fisher's exact test p-value"),
      caption = 'Cut-off lines drawn at equivalents of p=0.05, p=0.01, p=0.001') +
    
    geom_hline(yintercept = c(-log10(0.05), -log10(0.01), -log10(0.001)),
               linetype = c("dotted", "longdash", "solid"),
               colour = c("black", "black", "black"),
               size = c(0.5, 1.5, 3)) +
    
    theme_bw(base_size = 24) +
    theme(
      legend.position = 'right',
      legend.background = element_rect(),
      plot.title = element_text(angle = 0, size = 16, face = 'bold', vjust = 1),
      plot.subtitle = element_text(angle = 0, size = 14, face = 'bold', vjust = 1),
      plot.caption = element_text(angle = 0, size = 12, face = 'bold', vjust = 1),
      
      axis.text.x = element_text(angle = 0, size = 12, face = 'bold', hjust = 1.10),
      axis.text.y = element_text(angle = 0, size = 12, face = 'bold', vjust = 0.5),
      axis.title = element_text(size = 12, face = 'bold'),
      axis.title.x = element_text(size = 12, face = 'bold'),
      axis.title.y = element_text(size = 12, face = 'bold'),
      axis.line = element_line(colour = 'black'),
      
      #Legend
      legend.key = element_blank(), 
      legend.key.size = unit(1, "cm"), 
      legend.text = element_text(size = 14, face = "bold"), # Text size
      title = element_text(size = 14, face = "bold")) +
    
    coord_flip()
  print(gg1)
  dev.off()
}
class(df)
heatmap_plot<-function(norm.counts_,gene_list,output){
  df<-norm.counts_[gene_list,]
  z<-t(apply(df,1,scale))
  colnames(z)<-colnames(df)
  pdf(paste0(output,".pdf"))
  p<-Heatmap(z,cluster_rows=F,cluster_columns=F,column_labels=colnames(z),name="z-score", 
             column_names_gp = grid::gpar(fontsize = 3),row_names_gp = grid::gpar(fontsize = 5))
  print(p)
  dev.off()
}

for(photoperiod in photoperiods){
  modules <- unique(unname(unlist(clusters[photoperiod])))
  dist_mod<-list()
  for(i in 1:(length(modules)-1)){
    module <- modules[i]
    gene_list_<-rownames(clusters[clusters[photoperiod]==module,][photoperiod])
    GO_analysis(gene_list_,20,"BP",paste0("BP_",photoperiod,"_module_",module))
    GO_analysis(gene_list_,20,"MF",paste0("MF_",photoperiod,"_module_",module))
    heatmap_plot(mean_counts%>%dplyr::select(-gene_id),gene_list_,paste0("heatmap_",photoperiod,"_module_",module))
  }

}


