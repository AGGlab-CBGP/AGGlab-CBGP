#!/bin/R

library(dplyr)

data<-read.table("03.1.centered_CTs.tsv",sep=',',header=TRUE)
gene_list<-unique(data$gene)
length(gene_list)
conds<-unique(data$cond)
genots<-unique(data$genot)
CTs<-unique(data$CT)

for(con in conds){
  for(geno in genots){
    print(geno)
    print(con)
    genot_cond<-paste0(geno,"_",con)
    subdata<-data %>% filter(genot==geno, cond==con)
    exp<-list()
    time<-list()
    for(gen in gene_list){
      r1<-subdata %>% filter(gene==gen,rep=='r1') %>% select(vst) %>% unname() %>% unlist()
      r2<-subdata %>% filter(gene==gen,rep=='r2') %>% select(vst) %>% unname() %>% unlist()
      time_<-subdata %>% filter(gene==gen,rep=='r2') %>% select(newct) %>% unname() %>% unlist()
      exp[[gen]]<-c(r1,r2)
      time[[gen]]<-c(time_)
    }
    df <- data.frame(exp)
    df_time<-data.frame(time)
    write.csv(df, paste0(geno,"_",con,".csv"))
    write.csv(df_time, paste0(geno,"_",con,"_time.csv"))
  }
}
