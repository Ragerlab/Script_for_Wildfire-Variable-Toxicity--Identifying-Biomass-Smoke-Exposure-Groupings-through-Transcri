rm(list=ls())

library(vegan)
library(janitor)
library(xlsx)
library(readxl)
library(reshape2)
library(tidyverse)
library(cluster)
library(pheatmap)
library(factoextra)


setwd("/Users/lkoval/IEHS Dropbox/Rager Lab/Lauren_Koval/LK_Lab_Notebook/Projects/Transcriptomic_Similarity/Sensitivity_Analysis/Experiment_3/input/24h")


#Create dataframe to store all relevent data from all input files (gene, adjusted pvalue, log2 fold change, and the exposure condition)
df24 <- data.frame(gene=character(),
                  padj=numeric(),
                  log2FoldChange=numeric(),
                  expo_cond=character())


#Get list of input files from working directory
files_24 <- list.files()

#Read in the statistical results of all exposure conditions and store in previously established dataframe
for(f in files_24){
  cond <- str_split(f,"_")[[1]][5]
  temp <- read_csv(f)
  temp <- temp %>% dplyr::rename(gene=...1) %>% select(gene,padj,log2FoldChange) %>% mutate(expo_cond=cond)
  df24 <- rbind(df24,temp)
}  


setwd("/Users/lkoval/IEHS Dropbox/Rager Lab/Lauren_Koval/LK_Lab_Notebook/Projects/Transcriptomic_Similarity/Sensitivity_Analysis/Experiment_3")




#Filter for all rows that show a gene is significantly differentially expressed, here defined as a pval <0.1 and FC >=2.5, then pull a list of all unique genes that pass this filter
full_df <- as.data.frame(df24)


sig <- full_df %>% filter(expo_cond!="LPS") %>%  filter(padj<0.1 & abs(log2FoldChange)>=log2(1.5))
sig <- unique(sig$gene)


#Filter the original dataset for the list of unique genes, then create a column for a binary expression indicated by a cutpoint of padj<0.1 and a fold change of 2.5
#Format the data for clustering by exposure condition such that the rows are the exposure conditions,columns are the genes, and the values are
#the binary expression indicators.
df_sig <- full_df %>%
  filter(gene %in% all_of(sig)) %>% 
  mutate(padj=if_else(is.na(padj), 100,padj)) %>% 
  mutate(expr=if_else(padj<0.1 & abs(log2FoldChange)>=log2(1.5) ,1,0)) %>%
  select(!c("padj","log2FoldChange")) %>% 
  pivot_wider(names_from = gene, values_from = expr) %>%
  column_to_rownames("expo_cond")

# Create a jaccard distance matrix for exposure conditions
samp_comp <- mutate_all(df_sig, function(x) as.numeric(as.character(x)))
D_samp <- vegdist(as.matrix(samp_comp),method = "jaccard")

# make a table of jaccard similarities
sim_mat <- as.data.frame(as.matrix(D_samp))
sim_mat <- sim_mat %>% mutate_all(function(x) 1-x)

#read in file with FC >= 2.5 exposure condition groupings and specify where splits should go in heatmap
ord <- read_csv("input/24h_similarity_lower_tri.csv")
seprows <- c(4,5,7,10)


#Pull the filtered binary data and prepare it for clustering on genes
gene_comp <- samp_comp
gene_comp <- as.data.frame(t(gene_comp))

#Create a jaccard distance matrix for genes
D_gene <- vegdist(as.matrix(gene_comp),method = "jaccard")

#Determine the optimal number of clusters of genes
# Sys.time()
# wss_gene <- fviz_nbclust(gene_comp, diss = D_gene, method = "wss", FUN=hcut, hc_func="diana", k.max=100)
# wssplot_gene <- wss_gene$data
# wssplot_gene <- wssplot_gene %>% ggplot(aes(x=clusters, y=y))+geom_point()+geom_line(group=1)+labs(x="number of clusters k", y="total within sum of squares")
# ggsave("figures/24h_optimal_clusters_gene_wss_FC_1_5.png", height=10, width=16, unit="in",plot=wssplot_gene)
# Sys.time()
# silho_gene <- fviz_nbclust(gene_comp, diss = D_gene, method = "silhouette", FUN=hcut, hc_func="diana", k.max=100)
# silhoplot_gene <- silho_gene$data
# silhoplot_gene <- silhoplot_gene %>% ggplot(aes(x=clusters, y=y))+geom_point()+geom_line(group=1)+labs(x="number of clusters k", y="average silhouette width")
# ggsave("figures/24h_optimal_clusters_gene_silhouette_FC_1_5.png", height=10, width=16, unit="in",plot=silhoplot_gene)
# Sys.time()


#Run hierarchical clustering for genes
cluster_gene <-diana(D_gene, diss=TRUE)
plot(cluster_gene)

#Select 55 clusters and get cluster assignment for each gene
ncut_gene <- 55
cluster_assignments_gene <- cutree(cluster_gene, k = ncut_gene)

#Number of genes in each cluster
k_gene <-table(cluster_assignments_gene)

#Add the cluster assignments to genes and arrange the dataframe by the cluster
gene_comp$cluster <- cluster_assignments_gene
gene_comp <- gene_comp %>% arrange(cluster)


#specify new column order in exposure condition dataframe to reflect clustered genes
new_col_order <- rownames(gene_comp)
samp_comp <- samp_comp[new_col_order]

#Organize data such that the clusters are the same as the FC>=2.5 groupings
samp_comp <- samp_comp %>% rownames_to_column("cond") %>% arrange(factor(cond, levels=ord$ind))

hm_data <- as.data.frame(samp_comp)
hm_data <- hm_data %>% adorn_totals(where = "col", ... =!c("cond")) %>%  mutate(new_rownames=paste0(cond,"_",Total,"DEGs")) %>% column_to_rownames("cond") 
new_rownames <- hm_data$new_rownames
hm_data <- hm_data %>% select(!c("new_rownames","Total"))



#Make heatmap
pheatmap(hm_data, main="24h (padj<0.1 & |FC|>=1.5)",
         cluster_rows=FALSE,
         cluster_cols = FALSE,
         show_rownames = TRUE,
         show_colnames = FALSE,
         color = c("gray75","purple2"),
         breaks = c(0, 0.5, 1),
         legend = FALSE,
         cellheight = 50,
         cellwidth = 1.7,
         fontsize = 10,
         angle_col = 45,
         gaps_row = seprows,
         labels_row = new_rownames,
         filename = "figures/24h_omics_heatmap_FC_1_5.png",
         height = 24,
         width = 36,
         border_color = NA)



# format table of jaccard similarities between exposure conditions and include column of the number of overlapping genes that are significant for each
sim_mat <- sim_mat %>% rownames_to_column("cond") %>% arrange(factor(cond,levels = rownames(hm_data))) %>% column_to_rownames("cond") %>% select(rownames(hm_data))

sim_mat[upper.tri(sim_mat,diag=FALSE)] <- 17

sim_mat_long <- sim_mat %>% rownames_to_column("ind") %>% pivot_longer(!ind, names_to = "comp",values_to = "jaccard_similarity") %>% mutate(merge_col=paste0(ind,"_",comp))

degs <- as.data.frame(t(hm_data))
conds <- colnames(degs)

res <- data.frame(ind=character(),
                  comp=character(),
                  inxn=character())

for(i in 1:length(conds)){
  ind <- conds[i]
  for(j in 1:length(conds)){
    comp <- conds[j]
    temp <- degs %>% select(ind,comp)
    overlap <- temp %>% mutate(s=rowSums(.)) %>% filter(s==2)
    if(ncol(temp)==1){
      val <- sum(temp[1])
    }
    else{
      val <- nrow(overlap)
    }
    r <- c(ind, comp, val)
    res <- rbind(res, r)
  }
}

colnames(res) <- c("ind","comp","overlap")

res <- res %>% mutate(merge_col=paste0(ind,"_",comp)) %>% select(merge_col,overlap)

final_tab <- merge(sim_mat_long,res, by="merge_col", all.x = TRUE)
final_tab <- final_tab %>% select(!merge_col) %>%
  mutate(overlap=ifelse(jaccard_similarity==17," ",overlap)) %>% 
  mutate(jaccard_similarity=ifelse(jaccard_similarity==17," ",jaccard_similarity)) %>% 
  pivot_wider(names_from = "comp", values_from = c("jaccard_similarity","overlap")) %>% 
  select(c(1,2,13,3,14,4,15,5,16,12,23,8,19,9,20,6,17,7,18,10,21,11,22)) %>% 
  dplyr::slice(match(rownames(sim_mat),ind))


write_csv(final_tab, "output/24h_similarity_lower_tri_FC_1_5.csv")
