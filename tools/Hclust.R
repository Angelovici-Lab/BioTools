#!/usr/bin/env Rscript

#######################################################################
## Load packages
#######################################################################
rm(list = ls())

library(jpeg)
library(dplyr)
library(tidyr)
library(tibble)
library(stringr)
library(purrr)
library(tidyselect)
library(scales)
library(argparse)
library(pheatmap)
library(data.table)
library(reshape2)
library(ggforce)
library(ComplexHeatmap)

set.seed(1)


#######################################################################
## Argument Parser
#######################################################################

parser <- argparse::ArgumentParser()

parser$add_argument("-i", type="character", help="Input file path", required=TRUE)
parser$add_argument("-o", type="character", help="Output folder path", required=TRUE)
parser$add_argument("-k", type="integer", default=5, help="Number of clusters  (optional; default value is 5)")

args <- parser$parse_args()

input <- args$i
output <- args$o
k <- args$k

if(!file.exists(input)){
  print("Invalid input file. Exiting ...")
  quit(status = 0)
}

if(!dir.exists(output)){
  dir.create(output, recursive = TRUE)
}

## Load data
if(endsWith(input, "csv")){
  df <- read.csv(file = input, header = TRUE)
} else if (endsWith(input, "txt")){
  df <- read.delim(file = input, header=TRUE)
} else{
  print("Invalid input file. Exiting ...")
  quit(status = 0)
}

cat("\n Raw data: \n")
print(head(df))


#######################################################################
## Hclust code starts here
#######################################################################

## making the first col as row.names (convert data to martix)
r_names <- df[,1]
mat_data <- data.matrix(df[,2:ncol(df)])
rownames(mat_data) <- r_names

## scale the data
scaledata <- t(scale(t(mat_data)))

cat("\n Scaled data: \n")
print(head(scaledata))

## finding optimum cluster
wss <- (nrow(scaledata)-1)*sum(apply(scaledata,2,var))
for (i in 2:20){
  wss[i] <- sum(kmeans(scaledata,centers=i)$withinss)
}

cat("\n Save image... \n")
jpeg(file.path(output, "num_of_clusters.jpg"))
plot(data.frame(1:length(wss), wss), type="b", xlab="Number of Clusters", ylab="Within groups sum of squares")
dev.off()
cat("\n")

## calculate distance and then do hierarchial cluster, euclidean method can be customized to other method if you like
distance <- dist(scaledata, method = "euclidean")


## HClust of genes using distance, again method can be customized
gene_hclust <- hclust(distance, method = "complete")


## We can use the cutree() function do this dendrogram "cutting".
## For example, if you want to cut it into 5 groups, you would simply do:
## for more info:  https://tavareshugo.github.io/data-carpentry-rnaseq/04b_rnaseq_clustering.html
gene_cluster <- cutree(gene_hclust, k = k) %>%
  enframe() %>%
  arrange_at(1) %>%
  as.data.frame()
colnames(gene_cluster)[1] <- colnames(df)[1]

cat("\n Gene clusters: \n")
print(head(gene_cluster))

## do inner join to join cluster number and data
trans_cts_cluster <- as.data.frame(scaledata) %>%
  rownames_to_column(var = colnames(df)[1]) %>%
  inner_join(gene_cluster, by = colnames(df)[1]) %>%
  arrange_at(1) %>%
  arrange(value) %>%
  as.data.frame()

cat("\n Scaled data and cluster value: \n")
print(head(trans_cts_cluster))

write.csv(trans_cts_cluster, file.path(output, "trans_cts_cluster.csv"), row.names=FALSE)

trans_cts_cluster_longer <- as.data.frame(scaledata) %>%
  rownames_to_column(var = colnames(df)[1]) %>%
  pivot_longer(cols = -c(colnames(df)[1]), names_to = "key", values_to = "val") %>%
  inner_join(gene_cluster, by = colnames(df)[1]) %>%
  arrange_at(1) %>%
  arrange(value) %>%
  as.data.frame()

cat("\n Scaled and pivoted data with cluster value: \n")
print(head(trans_cts_cluster_longer))

trans_cts_cluster_longer$key <- factor(
  trans_cts_cluster_longer$key,
  levels = str_sort(unique(trans_cts_cluster_longer$key), numeric = TRUE)
)

## Finally plot the cluster
p <- trans_cts_cluster_longer %>%
  ggplot(aes(key,val)) +
  geom_line(aes(group = Gene))+
  geom_line(stat = "summary", fun.y = "mean", colour = "brown", size = 1.5, aes(group = 1)) +
  facet_wrap( facets = 'value')

cat("\n Save image... \n")
ggsave(
  filename = "clusters.png",
  plot = p,
  path = output
)

cat("\n Save image... \n")
jpeg(file.path(output, "heatmap.jpg"), width = 800, height = 800)
Heatmap(scaledata, show_row_names = FALSE)
dev.off()
cat("\n")
