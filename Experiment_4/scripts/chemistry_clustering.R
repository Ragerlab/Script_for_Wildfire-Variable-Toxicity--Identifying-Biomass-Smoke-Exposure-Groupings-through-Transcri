rm(list=ls())

#load libraries
library(janitor)
library(xlsx)
library(readxl)
library(reshape2)
library(tidyverse)
library(pheatmap)
library(vegan)
library(factoextra)
library(cluster)


setwd("Experiment_4")

#read in chemistry data and select the columns that reflect the chemical and the exposure conditions
chems <- read_xlsx("input/Chemistry_Calcs_012920.xlsx", sheet = "CleanedData")
chems <- chems %>% select(Chemical | contains("Flaming") | contains("Smoldering")) %>% column_to_rownames("Chemical")

#transpose so we can z-score center and standardize each of the chemicals
chems <- as.data.frame(t(scale(t(chems))))

#make a distance matrix of chems
D_chem <- vegdist(as.matrix(chems),method = "euclidean")
 
# sim_mat <- as.data.frame(as.matrix(D_samp))
# sim_mat <- sim_mat %>% mutate_all(function(x) 1-x)

#Determine the optimal number of chemical clusters. 6 looks pretty good.
wss_chem <- fviz_nbclust(chems, diss = D_chem, method = "wss", FUN=hcut, hc_func="diana", k.max=30)
wssplot_chem <- wss_chem$data
wssplot_chem <- wssplot_chem %>% ggplot(aes(x=clusters, y=y))+geom_point()+geom_line(group=1)+labs(x="number of clusters k", y="total within sum of squares")
ggsave("figures/optimal_clusters_chems_wss.png", height=10, width=16, unit="in",plot=wssplot_chem)

silho_chem <- fviz_nbclust(chems, diss = D_chem, method = "silhouette", FUN=hcut, hc_func="diana", k.max=30)
silhoplot_chem <- silho_chem$data
silhoplot_chem <- silhoplot_chem %>% ggplot(aes(x=clusters, y=y))+geom_point()+geom_line(group=1)+labs(x="number of clusters k", y="average silhouette width")
ggsave("figures/optimal_clusters_chems_silhouette.png", height=10, width=16, unit="in",plot=silhoplot_chem)



#prepare data for clustering on samples
samp_comp <- as.data.frame(t(chems))

#make a distance matrix for samples
D_samp <- vegdist(as.matrix(samp_comp),method = "euclidean")
D_df <- as.data.frame(as.matrix(D_samp))
write.csv(D_df, "output/euclidean_distance_matrix.csv")

#Determine the optimal number of sample clusters. 4 looks pretty good.
wss_samp <- fviz_nbclust(samp_comp, diss = D_samp, method = "wss", FUN=hcut, hc_func="diana", k.max=9)
wssplot_samp <- wss_samp$data
wssplot_samp <- wssplot_samp %>% ggplot(aes(x=clusters, y=y))+geom_point()+geom_line(group=1)+labs(x="number of clusters k", y="total within sum of squares")
ggsave("figures/optimal_clusters_samp_wss.png", height=10, width=16, unit="in",plot=wssplot_samp)

silho_samp <- fviz_nbclust(samp_comp, diss = D_samp, method = "silhouette", FUN=hcut, hc_func="diana", k.max=9)
silhoplot_samp <- silho_samp$data
silhoplot_samp <- silhoplot_samp %>% ggplot(aes(x=clusters, y=y))+geom_point()+geom_line(group=1)+labs(x="number of clusters k", y="average silhouette width")
ggsave("figures/optimal_clusters_samp_silhouette.png", height=10, width=16, unit="in",plot=silhoplot_samp)


#Organize data for heatmap
hm_data <- t(chems)

#Make heatmap with 4 clusters of samples and 6 clusters of chemicals
pheatmap(hm_data,
         cluster_rows=TRUE,
         cluster_cols = TRUE,
         clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean",
         show_rownames = TRUE,
         show_colnames = TRUE,
         treeheight_row = 35,
         treeheight_col = 35,
         cutree_rows = 4,
         cutree_cols = 6,
         legend = TRUE,
         cellwidth = 9,
         cellheight = 25.5,
         fontsize = 7,
         angle_col = 45,
         filename = "figures/clustered_chems.png",
         width = 18,
         height = 8)





########################################## Chem Organized to Match Omic Heatmap Clusters ####################################

#Reorganize burn conditions to match 4h transcriptomic clusters
org_4h <- as.data.frame(t(chems))
org_4h <- org_4h %>%
  rownames_to_column("Cond") %>%
  mutate(Cond=factor(Cond, levels = c("Peat_Flaming",
                                      "Eucalyptus_Flaming",
                                      "Eucalyptus_Smoldering",
                                      "Pine_Needles_Flaming",
                                      "Red_Oak_Flaming",
                                      "Pine_Flaming",
                                      "Pine_Needles_Smoldering",
                                      "Pine_Smoldering",
                                      "Peat_Smoldering",
                                      "Red_Oak_Smoldering"))) %>%
  arrange(Cond) %>% 
  column_to_rownames("Cond")

#Make heatmap 
pheatmap(org_4h,
         main = "4h match omics",
         cluster_rows=FALSE,
         cluster_cols = TRUE,
         clustering_distance_cols = "euclidean",
         show_rownames = TRUE,
         show_colnames = TRUE,
         treeheight_col = 35,
         cutree_cols = 6,
         gaps_row = c(1,5,6,8),
         legend = TRUE,
         cellwidth = 9,
         cellheight = 25.5,
         fontsize = 7,
         angle_col = 45,
         filename = "figures/4h_organized_chems.png",
         width = 18,
         height = 8)





#Reorganize burn conditions to match 24h transcriptomic clusters
org_24h <- as.data.frame(t(chems))
org_24h <- org_24h %>%
  rownames_to_column("Cond") %>%
  mutate(Cond=factor(Cond, levels = c("Eucalyptus_Flaming",
                                      "Eucalyptus_Smoldering",
                                      "Peat_Flaming",
                                      "Red_Oak_Smoldering",
                                      "Pine_Needles_Flaming",
                                      "Pine_Needles_Smoldering",
                                      "Peat_Smoldering",
                                      "Pine_Flaming",
                                      "Pine_Smoldering",
                                      "Red_Oak_Flaming"))) %>%
  arrange(Cond) %>% 
  column_to_rownames("Cond")

#Make heatmap 
pheatmap(org_24h,
         main = "24h match omics",
         cluster_rows=FALSE,
         cluster_cols = TRUE,
         clustering_distance_cols = "euclidean",
         show_rownames = TRUE,
         show_colnames = TRUE,
         treeheight_col = 35,
         cutree_cols = 6,
         gaps_row = c(3,4,6,9),
         legend = TRUE,
         cellwidth = 9,
         cellheight = 25.5,
         fontsize = 7,
         angle_col = 45,
         filename = "figures/24h_organized_chems.png",
         width = 18,
         height = 8)



