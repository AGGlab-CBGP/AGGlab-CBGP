arser_data<-read.table("TFM/02.2.ARSER_parameters.avg_expr.tsv",sep='\t',header=TRUE)

folder <- 'centered_data/peaks'
files <- list.files(folder, full.names = TRUE)
AUC_shape<-data.frame(gene=gene_list);AUC_exp<-data.frame(gene=gene_list);skewness<-data.frame(gene=gene_list);kurtosis<-data.frame(gene=gene_list)
width<-data.frame(gene=gene_list);time_increasing<-data.frame(gene=gene_list);time_decreasing<-data.frame(gene=gene_list)
height<-data.frame(gene=gene_list);prominence<-data.frame(gene=gene_list);residuals<-data.frame(gene=gene_list);
centers<-data.frame(gene=gene_list);coefs<-data.frame(gene=gene_list);start<-data.frame(gene=gene_list);
end<-data.frame(gene=gene_list)
for (file in files) {
  genotype<-unlist(strsplit(unlist(strsplit(basename(file),".csv"))[1],"_"))[2]
  photoperiod<-unlist(strsplit(unlist(strsplit(basename(file),".csv"))[1],"_"))[3]
  read<-read.table(file,header=TRUE,sep=',')
  
  #AUC shape (normalized exp)
  temp_AUC<-data.frame(nrow=rep(0,length(gene_list)))
  temp_AUC$X1<-read$gene
  temp_AUC$X2<-read$AUC_shape
  temp_AUC <- temp_AUC  %>% select(-c("nrow"))
  colnames(temp_AUC)<-c("gene",paste0(genotype,"_",photoperiod))
  AUC_shape<-merge(AUC_shape,temp_AUC,by = intersect("gene","gene"))
  
  #AUC exp
  temp_AUC_exp<-data.frame(nrow=rep(0,length(gene_list)))
  temp_AUC_exp$X1<-read$gene
  temp_AUC_exp$X2<-read$AUC_exp
  temp_AUC_exp <- temp_AUC_exp  %>% select(-c("nrow"))
  colnames(temp_AUC_exp)<-c("gene",paste0(genotype,"_",photoperiod))
  AUC_exp<-merge(AUC_exp,temp_AUC_exp,by = intersect("gene","gene"))
  
  #width
  temp_width<-data.frame(nrow=rep(0,length(gene_list)))
  temp_width$X1<-read$gene
  temp_width$X2<-read$widths
  temp_width <- temp_width  %>% select(-c("nrow"))
  colnames(temp_width)<-c("gene",paste0(genotype,"_",photoperiod))
  width<-merge(width,temp_width,by = intersect("gene","gene"))
  
  #time increasing
  temp_time_increasing<-data.frame(nrow=rep(0,length(gene_list)))
  temp_time_increasing$X1<-read$gene
  temp_time_increasing$X2<-read$widths
  temp_time_increasing <- temp_time_increasing  %>% select(-c("nrow"))
  colnames(temp_time_increasing)<-c("gene",paste0(genotype,"_",photoperiod))
  time_increasing<-merge(time_increasing,temp_time_increasing,by = intersect("gene","gene"))
  
  #time decreasing
  temp_time_decreasing<-data.frame(nrow=rep(0,length(gene_list)))
  temp_time_decreasing$X1<-read$gene
  temp_time_decreasing$X2<-read$widths
  temp_time_decreasing <- temp_time_decreasing  %>% select(-c("nrow"))
  colnames(temp_time_decreasing)<-c("gene",paste0(genotype,"_",photoperiod))
  time_decreasing<-merge(time_decreasing,temp_time_decreasing,by = intersect("gene","gene"))
  
  #heights
  temp_height<-data.frame(nrow=rep(0,length(gene_list)))
  temp_height$X1<-read$gene
  temp_height$X2<-read$heights
  temp_height <- temp_height  %>% select(-c("nrow"))
  colnames(temp_height)<-c("gene",paste0(genotype,"_",photoperiod))
  height<-merge(height,temp_height,by = intersect("gene","gene"))
  
  #prominence
  temp_prominence<-data.frame(nrow=rep(0,length(gene_list)))
  temp_prominence$X1<-read$gene
  temp_prominence$X2<-read$prominences
  temp_prominence <- temp_prominence  %>% select(-c("nrow"))
  colnames(temp_prominence)<-c("gene",paste0(genotype,"_",photoperiod))
  prominence<-merge(prominence,temp_prominence,by = intersect("gene","gene"))
  
  #skewness
  temp_sknwss<-data.frame(nrow=rep(0,length(gene_list)))
  temp_sknwss$X1<-read$gene
  temp_sknwss$X2<-read$skewness
  temp_sknwss <- temp_sknwss %>% select(-c("nrow"))
  colnames(temp_sknwss)<-c("gene",paste0(genotype,"_",photoperiod))
  skewness<-merge(skewness,temp_sknwss,by = intersect("gene","gene"))
  
  #kurtosis
  temp_kurtosis<-data.frame(nrow=rep(0,length(gene_list)))
  temp_kurtosis$X1<-read$gene
  temp_kurtosis$X2<-read$kurtosis
  temp_kurtosis <- temp_kurtosis %>% select(-c("nrow"))
  colnames(temp_kurtosis)<-c("gene",paste0(genotype,"_",photoperiod))
  kurtosis<-merge(kurtosis,temp_kurtosis,by = intersect("gene","gene"))
  
  #residual
  temp_residuals<-data.frame(nrow=rep(0,length(gene_list)))
  temp_residuals$X1<-read$gene
  temp_residuals$X2<-read$residuals
  temp_residuals <- temp_residuals %>% select(-c("nrow"))
  colnames(temp_residuals)<-c("gene",paste0(genotype,"_",photoperiod))
  residuals<-merge(residuals,temp_residuals,by = intersect("gene","gene"))
  
  #center
  temp_centers<-data.frame(nrow=rep(0,length(gene_list)))
  temp_centers$X1<-read$gene
  temp_centers$X2<-read$centers
  temp_centers <- temp_centers %>% select(-c("nrow"))
  colnames(temp_centers)<-c("gene",paste0(genotype,"_",photoperiod))
  centers<-merge(centers,temp_centers,by = intersect("gene","gene"))
  
  #coef
  temp_coefs<-data.frame(nrow=rep(0,length(gene_list)))
  temp_coefs$X1<-read$gene
  temp_coefs$X2<-read$coefs
  temp_coefs <- temp_coefs %>% select(-c("nrow"))
  colnames(temp_coefs)<-c("gene",paste0(genotype,"_",photoperiod))
  coefs<-merge(coefs,temp_coefs,by = intersect("gene","gene"))
  
  #start
  temp_starts<-data.frame(nrow=rep(0,length(gene_list)))
  temp_starts$X1<-read$gene
  temp_starts$X2<-read$start
  temp_starts <- temp_starts %>% select(-c("nrow"))
  colnames(temp_starts)<-c("gene",paste0(genotype,"_",photoperiod))
  start<-merge(start,temp_starts,by = intersect("gene","gene"))
  
  #end
  temp_ends<-data.frame(nrow=rep(0,length(gene_list)))
  temp_ends$X1<-read$gene
  temp_ends$X2<-read$end
  temp_ends <- temp_ends %>% select(-c("nrow"))
  colnames(temp_ends)<-c("gene",paste0(genotype,"_",photoperiod))
  end<-merge(end,temp_ends,by = intersect("gene","gene"))
}

arser_data<- arser_data %>% filter(CycID %in% gene_list)
new_data <- data.frame(CycID=arser_data$CycID)
descriptors<-list(AUC_shape,AUC_exp,width,time_increasing,time_decreasing,skewness,kurtosis,
                  prominence,height,centers,coefs,residuals,start,end)
i<-1
while(i<length(colnames(arser_data))){
  group_name<-unlist(strsplit(colnames(arser_data)[i+1],"\\.[a-z]*"))[1]
  start_col<-i+6
  end_col<-start_col+2

  new_col_names<-c(paste0(group_name,".AUC_shape"),paste0(group_name,".AUC_exp"),
                   paste0(group_name,".width"),paste0(group_name,".time_increasing"),
                   paste0(group_name,".time_descreasing"),paste0(group_name,".skewness"),
                   paste0(group_name,".kurtosis"),paste0(group_name,".prominence"),paste0(group_name,".height"),
                   paste0(group_name,".centers"),paste0(group_name,".coefs"),paste0(group_name,".residuals"),
                   paste0(group_name,".start"),paste0(group_name,".end"))
  
  for(a in 1:5){
    name_col<-colnames(arser_data)[i+a]
    new_data<-cbind(new_data,arser_data[name_col])
  }
  for(r in 1:14){
    new_data[new_col_names[r]]<-descriptors[[r]][group_name]
  }
  
  
  i<-i+5
}

#Mutate NaN to NA
new_data <- new_data %>%
  mutate_all(~ ifelse(is.nan(.), "NA", .))

write.table(new_data,"oscillation_descriptors.tsv",sep='\t',row.names = FALSE)
