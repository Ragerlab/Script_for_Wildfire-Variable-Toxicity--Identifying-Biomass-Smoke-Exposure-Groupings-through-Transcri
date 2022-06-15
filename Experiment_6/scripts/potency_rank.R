rm(list=ls())

#load libraries
library(janitor)
library(xlsx)
library(readxl)
library(reshape2)
library(tidyverse)


setwd("Experiment_6")


make_rank <- function(df){
  u <- sort(unique(df[["n"]]))
  map <- data.frame(n=u)
  map$score <- 1:nrow(map)
  final <- merge(df, map, by="n")
  return(final)
}


deg_4h <- read.csv("input/DEGs_4h.csv")

sig_4h <- deg_4h %>% filter(expo_cond!="LPS")
sig_4h <- unique(sig_4h$gene)

deg_4h_counts <- deg_4h %>% filter(gene %in% all_of(sig_4h)) %>% count(expo_cond) %>% arrange(desc(n))

deg_4h_rank <- make_rank(deg_4h_counts)
deg_4h_rank <- deg_4h_rank %>% rename(deg_4h=n, deg_4h_score=score)


deg_24h <- read.csv("input/DEGs_24h.csv")

sig_24h <- deg_24h %>% filter(expo_cond!="LPS")
sig_24h <- unique(sig_24h$gene)

deg_24h_counts <- deg_24h %>% filter(gene %in% all_of(sig_24h)) %>% count(expo_cond) %>% arrange(desc(n))

deg_24h_rank <- make_rank(deg_24h_counts)
deg_24h_rank <- deg_24h_rank %>% rename(deg_24h=n, deg_24h_score=score)

tox <- read.csv("input/cardiopulm_endpoints.csv")
tox_4h <- tox %>%
  filter(Sample.Collection.Post.Exposure..h.==4) %>%
  select(!c("Sample.Collection.Post.Exposure..h.", "merge_col")) %>%
  column_to_rownames("Biomass.Burn")

tox_4h$n <- rowSums(tox_4h)
tox_4h <- tox_4h %>% rownames_to_column("expo_cond") %>% select(expo_cond,n)
tox_4h_rank <- make_rank(tox_4h)
tox_4h_rank <- tox_4h_rank %>% rename(tox_4h=n, tox_4h_score=score)

tox_24h <- tox %>%
  filter(Sample.Collection.Post.Exposure..h.==24) %>%
  select(!c("Sample.Collection.Post.Exposure..h.", "merge_col")) %>%
  column_to_rownames("Biomass.Burn")

tox_24h$n <- rowSums(tox_24h)
tox_24h <- tox_24h %>% rownames_to_column("expo_cond") %>% select(expo_cond,n)
tox_24h_rank <- make_rank(tox_24h)
tox_24h_rank <- tox_24h_rank %>% rename(tox_24h=n, tox_24h_score=score)

master <- merge(deg_4h_rank, deg_24h_rank, by="expo_cond")
master <- merge(master, tox_4h_rank, by="expo_cond")
master <- merge(master, tox_24h_rank, by="expo_cond")

master <- master %>%
  column_to_rownames("expo_cond") %>%
  mutate(final_score=deg_4h_score+deg_24h_score+tox_4h_score+tox_24h_score) %>% 
  rownames_to_column("expo_cond") %>%
  arrange(desc(final_score))

write_csv(master, "output/overall_biomass_burn_ranking.csv")




master_4h <- merge(deg_4h_rank, tox_4h_rank, by="expo_cond")

master_4h <- master_4h %>%
  column_to_rownames("expo_cond") %>%
  mutate(final_score=deg_4h_score+tox_4h_score) %>% 
  rownames_to_column("expo_cond") %>%
  arrange(desc(final_score))

write_csv(master_4h, "output/4h_biomass_burn_ranking.csv")




master_24h <- merge(deg_24h_rank, tox_24h_rank, by="expo_cond")

master_24h <- master_24h %>%
  column_to_rownames("expo_cond") %>%
  mutate(final_score=deg_24h_score+tox_24h_score) %>% 
  rownames_to_column("expo_cond") %>%
  arrange(desc(final_score))

write_csv(master_24h, "output/24h_biomass_burn_ranking.csv")



