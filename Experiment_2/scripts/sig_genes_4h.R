rm(list=ls())

library(janitor)
library(xlsx)
library(readxl)
library(reshape2)
library(tidyverse)


setwd("/Users/lkoval/IEHS Dropbox/Rager Lab/Lauren_Koval/LK_Lab_Notebook/Projects/Transcriptomic_Similarity/Sensitivity_Analysis/Experiment_2/input/4h")


#Create dataframe to store all relevent data from all input files (gene, baseMean, adjusted pvalue, log2 fold change, and the exposure condition)
df4 <- data.frame(gene=character(),
                  baseMean=numeric(),
                  padj=numeric(),
                  log2FoldChange=numeric(),
                  expo_cond=character())


#Get list of input files from working directory
files_4 <- list.files()

#Read in the statistical results of all exposure conditions and store in previously established dataframe
for(f in files_4){
  cond <- str_split(f,"_")[[1]][5]
  temp <- read_csv(f)
  temp <- temp %>% dplyr::rename(gene=...1) %>% select(gene,baseMean,padj,log2FoldChange) %>% mutate(cond=cond)
  df4 <- rbind(df4,temp)
}  

#Change to main directory
setwd("/Users/lkoval/IEHS Dropbox/Rager Lab/Lauren_Koval/LK_Lab_Notebook/Projects/Transcriptomic_Similarity/Sensitivity_Analysis/Experiment_2")



get_sig_genes <- function(df,filt){
  #Filter for genes that differentially expressed  for any exposure condition
  full_df <- as.data.frame(df)
  full_df <- full_df %>% filter(padj<0.1 & abs(log2FoldChange)>=filt)
  
  
  #make mapping of terminology we want to use
  conds <- unique(full_df$cond)
  
  mapping <- data.frame(cond=conds,
                        expo_cond=c("Flaming Eucalyptus",
                                    "Smoldering Eucalyptus",
                                    "LPS",
                                    "Flaming Peat",
                                    "Smoldering Peat",
                                    "Flaming Pine",
                                    "Flaming Pine Needles",
                                    "Smoldering Pine Needles",
                                    "Smoldering Pine",
                                    "Flaming Red Oak",
                                    "Smoldering Red Oak"))
  
  mapping <- mapping %>% mutate(group_comp=paste0(expo_cond," vs. Saline"))
  
  #map in correct terminology to df of DEGs
  final <- merge(full_df, mapping, by="cond", all.x = TRUE)
  final <- final %>% mutate(timepoint="4h") %>% select(expo_cond,group_comp,timepoint,gene,baseMean,padj,log2FoldChange)
  return(final)
}




fc2.5 <-get_sig_genes(df4, log2(2.5))

#write out results
write_csv(fc2.5, "output/DEGs_4h.csv")




