#!/bin/R


#Alberto Gonzalez Delgado
#09/2023
#Centro de Biotecnologia y Genomica de Plantas (UPM-INIA/CSIC)

library(dplyr)
library(stringr)
library(purrr)

setwd("C:/Users/admin/Desktop/projects")
######################################## Separate genes in phase groups ##############################################

oscillation_data<-read.csv("TFM/02.2.ARSER_parameters.avg_expr.tsv",sep='\t',header=TRUE)

#Round phases
oscillation_data<- oscillation_data %>%
  mutate(across(ends_with(".adjphase"), round, digits = 0))
# Filter only oscillatory genes
filtered_data <- oscillation_data %>%
  filter(CycID %in% gene_list)


lists_genes<- vector("list", 24)
for (i in 1:24) {
  # Generate temporal vector
  temp_cycid <- c()
  
  # Search only in adjphase columns
  for (col_name in names(filtered_data[grep(".adjphase$", names(filtered_data))])) {
    # Save genes with phase == i 
    temp_cycid <- c(temp_cycid, filtered_data$CycID[filtered_data[[col_name]] == i])
  }
  
  #Eliminate duplicates
  lists_genes[[i]] <- unique(temp_cycid)
}

length(unlist(lists_genes[3]))

######################################################################### Generate lists .bed files ##############################################################################

##This file contains the position of the CDS
data_table <- read.delim("./input/02.1.vst_normalized_counts.tsv", header = TRUE)

for(i in 1:length(lists_genes)){
  #Access each list (each phase)
  gene_list<-lists_genes[[i]]
  #Filter only gene_list data 
  filtered_data <- data_table %>% filter(gene_id %in% gene_list)
  #Mutate start 1.5 kb and end to the original start 
  filtered_data <- filtered_data %>%
    mutate(original_start = start,
           original_end = end) %>%
    mutate(
      start = ifelse(strand == "+", start - 1500, original_end),
      end = ifelse(strand == "+", original_start, end + 1500)
    ) %>%
    select(-original_start, -original_end)
  #Generate bed data
  bed_data <- data.frame(
    ID = 1:nrow(filtered_data),
    chr = filtered_data$chr,
    start = filtered_data$start,
    end = filtered_data$end,
    strand = filtered_data$strand
  )
  
  #Generate .bed file
  file_name <- paste0("./motifs/phase_", i, ".bed")
  write.table(bed_data, file = file_name, sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
  
}


