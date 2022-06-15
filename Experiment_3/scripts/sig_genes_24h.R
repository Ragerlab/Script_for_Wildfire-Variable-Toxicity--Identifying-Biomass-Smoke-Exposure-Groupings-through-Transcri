rm(list=ls())

library(janitor)
library(xlsx)
library(readxl)
library(reshape2)
library(tidyverse)


setwd("Experiment_3/input/24h")


#Create dataframe to store all relevent data from all input files (gene, baseMean, adjusted pvalue, log2 fold change, and the exposure condition)
df24 <- data.frame(gene=character(),
                  baseMean=numeric(),
                  padj=numeric(),
                  log2FoldChange=numeric(),
                  expo_cond=character())


#Get list of input files from working directory
files_24 <- list.files()

#Read in the statistical results of all exposure conditions and store in previously established dataframe
for(f in files_24){
  cond <- str_split(f,"_")[[1]][5]
  temp <- read_csv(f)
  temp <- temp %>% dplyr::rename(gene=...1) %>% select(gene,baseMean,padj,log2FoldChange) %>% mutate(cond=cond)
  df24 <- rbind(df24,temp)
}  

#Change to main directory
setwd("Experiment_3")



get_sig_genes <- function(df,filt){
  #Filter for genes that differentially expressed for any exposure condition
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
  final <- final %>% mutate(timepoint="24h") %>% select(expo_cond,group_comp,timepoint,gene,baseMean,padj,log2FoldChange)
  return(final)
}

fc2 <-get_sig_genes(df24, log2(2))
fc1.5 <-get_sig_genes(df24, log2(1.5))


#write out results
write_csv(fc2, "output/DEGs_24h_FC_2.csv")
write_csv(fc1.5, "output/DEGs_24h_FC_1_5.csv")




