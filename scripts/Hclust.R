#!/usr/bin/env Rscript

#######################################################################
## Load packages
#######################################################################
rm(list = ls())

library(jpeg)
library(png)
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
library(pheatmap)

# Set seed
set.seed(1)


#######################################################################
## Input and output (Double check or make changes)
#######################################################################
# Specify input file (csv or tab delimited text file)
input <- file.path('/home/ycth8/data/projects/BioTools/data/Hclust.csv')

# Specify an output folder so that all the results can be stored into the folder
output <- file.path("/home/ycth8/data/projects/BioTools/output/07_24_2020")

# Number of cluster
k <- 5

# Total number of test cluster size
number_of_tests <- 20

clustering_distance <- "euclidean"

clustering_method <- "complete"


#######################################################################
## Prepare input and output
#######################################################################
# Check if the input file exists
# If it is not exists, the script will ends with a warning message
if(!file.exists(input)){
  print("Invalid input file. Exiting ...")
  quit(status = 0)
}

# If the output directory does not exists, it will be created.
if(!dir.exists(output)){
  dir.create(output, recursive = TRUE)
}

## Load data
# the input file extension should be csv or txt
if(endsWith(input, "csv")){
  df <- read.csv(file = input, header = TRUE)
} else if (endsWith(input, "txt")){
  df <- read.delim(file = input, header=TRUE)
} else{
  print("Invalid input file. Exiting ...")
  quit(status = 0)
}

# Display raw data
cat("\n Raw data: \n")
print(head(df))


#######################################################################
##
## Hclust code starts here
##
#######################################################################


#######################################################################
## Generate matrix and scaled matrix
#######################################################################

## making the first col as row.names (convert data to martix)
r_names <- df[,1]
mat_data <- data.matrix(df[,2:ncol(df)])
rownames(mat_data) <- r_names

# Show original data in matrix format
cat("\n Matrix data: \n")
print(head(mat_data))

## scale the data
scaledata <- t(scale(t(mat_data)))

# Show scaled data in matrix format
cat("\n Scaled data: \n")
print(head(scaledata))


#######################################################################
## Perform clustering
#######################################################################

## finding optimum cluster
wss <- (nrow(scaledata)-1)*sum(apply(scaledata,2,var))
if (number_of_tests > 1) {
  if (number_of_tests >= nrow(scaledata)) {
    number_of_tests <- nrow(scaledata) - 1
  }
  for (i in 2:number_of_tests){
    wss[i] <- sum(kmeans(scaledata,centers=i)$withinss)
  }
}


# Save the num_of_clusters image
cat("\n Save number of clusters plot... \n")
jpeg(file.path(output, "num_of_clusters.jpg"))
plot(data.frame(1:length(wss), wss), type="b", xlab="Number of Clusters", ylab="Within groups sum of squares")
dev.off()
cat("\n")

## calculate distance and then do hierarchial cluster, euclidean method can be customized to other method if you like
distance <- dist(scaledata, method = clustering_distance)

## HClust of genes using distance, again method can be customized
gene_hclust <- hclust(distance, method = clustering_method)

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


#######################################################################
## Add clusters to original data
#######################################################################

key_column <- colnames(df)[1]

original_data_df <- df %>%
  inner_join(gene_cluster, by = all_of(key_column)) %>%
  arrange_at(1) %>%
  arrange(value) %>%
  as.data.frame()

# Display scaled data and cluster value
cat("\n Original data and cluster value: \n")
print(head(original_data_df))

original_data_df_longer <- df %>%
  pivot_longer(cols = -all_of(key_column), names_to = "key", values_to = "val") %>%
  inner_join(gene_cluster, by = all_of(key_column)) %>%
  arrange_at(1) %>%
  arrange(value) %>%
  as.data.frame()

# Display scaled and pivoted data with cluster value
cat("\n Original and pivoted data with cluster value: \n")
print(head(original_data_df_longer))


#######################################################################
## Add cluster to scaled data
#######################################################################

key_column <- colnames(df)[1]

## do inner join to join cluster number and data
scaled_data_df <- as.data.frame(scaledata) %>%
  rownames_to_column(var = key_column) %>%
  inner_join(gene_cluster, by = all_of(key_column)) %>%
  arrange_at(1) %>%
  arrange(value) %>%
  as.data.frame()

# Display scaled data and cluster value
cat("\n Scaled data and cluster value: \n")
print(head(scaled_data_df))

# Save scaled data and cluster value
write.csv(scaled_data_df, file.path(output, "scaled_data_df.csv"), row.names=FALSE)

# Pivot table and join with gene information
scaled_data_df_longer <- as.data.frame(scaledata) %>%
  rownames_to_column(var = key_column) %>%
  pivot_longer(cols = -all_of(key_column), names_to = "key", values_to = "val") %>%
  inner_join(gene_cluster, by = all_of(key_column)) %>%
  arrange_at(1) %>%
  arrange(value) %>%
  as.data.frame()

# Display scaled and pivoted data with cluster value
cat("\n Scaled and pivoted data with cluster value: \n")
print(head(scaled_data_df_longer))


#######################################################################
## Plot cluster line plot
#######################################################################

# Add factor into the key column
scaled_data_df_longer$key <- factor(
  scaled_data_df_longer$key,
  levels = str_sort(unique(scaled_data_df_longer$key), numeric = TRUE)
)

## Finally plot the cluster
p <- scaled_data_df_longer %>%
  ggplot(aes(key,val)) +
  geom_line(aes(group = Gene))+
  geom_line(stat = "summary", colour = "brown", size = 1.5, aes(group = 1)) +
  facet_wrap( facets = 'value') +
  labs(x = "Mutant", y = "Values")

# Save cluster image
cat("\n Save facetted clusters image... \n")
ggsave(
  filename = "clusters_scaled_data.png",
  plot = p,
  path = output
)


for(i in sort(unique(scaled_data_df_longer$value))){
  p <- scaled_data_df_longer[scaled_data_df_longer$value == i, ] %>%
    ggplot(aes(key,val)) +
    #geom_line(aes(group = Gene))+
    geom_line(stat = "summary", colour = "brown", size = 1.5, aes(group = 1)) +
    labs(x = "Mutant", y = "Values")

  ggsave(
    filename = paste0(
      "clusters_scaled_data_",
      as.character(i),
      "__",
      as.character(nrow(scaled_data_df_longer[scaled_data_df_longer$value == i, ])),
      "_rows.png"
    ),
    plot = p,
    path = output
  )
}


#######################################################################
## Plot heatmaps
#######################################################################

# Save heatmap
cat("\n Save normal heatmap plot with scaled data... \n")
jpeg(file.path(output, "normal_heatmap_plot_with_scaled_data.jpg"), width = 800, height = 800)
Heatmap(scaledata, show_row_names = FALSE)
dev.off()
cat("\n")


for(i in sort(unique(gene_cluster$value))){
  tf <- (rownames(scaledata) %in% gene_cluster$Gene[gene_cluster$value==i])
  temp_scaledata <- data.matrix(scaledata[tf, ])
  pheatmap(
    temp_scaledata,
    fontsize=10,
    color=colorRampPalette(c("blue", "yellow", "red"))(299),
    cluster_rows=TRUE,
    cluster_cols=TRUE,
    fontsize_row=1,
    fontsize_col=10,
    scale="row",
    filename=file.path(
      output,
      paste0(
        "Heatmap_no_cluster_plot_with_scaled_data_",
        as.character(i),
        "__",
        as.character(nrow(temp_scaledata)),
        "_rows.jpg"
      )
    ),
    width = 20,
    height = 25
  )
}


# Plot and save heatmap using original data
#cat("\n Save heatmaps plot with original data... \n")
#pheatmap(
#  mat_data,
#  fontsize=10,
#  color=colorRampPalette(c("blue", "yellow", "red"))(299),
#  cluster_rows=TRUE,
#  cluster_cols=TRUE,
#  fontsize_row=1,
#  fontsize_col=10,
#  scale="row",
#  filename=file.path(output, "Heatmap_no_cluster_plot_with_original_data.jpg"),
#  width = 20,
#  height = 25
#)
#
#pheatmap(
#  mat_data,
#  fontsize=10,
#  color=colorRampPalette(c("blue", "yellow", "red"))(299),
#  cluster_rows=TRUE,
#  cluster_cols=TRUE,
#  cutree_rows = k,
#  clustering_distance_rows = clustering_distance,
#  clustering_method = clustering_method,
#  fontsize_row=1,
#  fontsize_col=10,
#  scale="row",
#  filename=file.path(output, "Heatmap_with_cluster_plot_with_original_data.jpg"),
#  width = 20,
#  height = 25
#)
#cat("\n")

# Plot and save heatmap using scaled data
cat("\n Save heatmaps plot with scaled data... \n")
pheatmap(
  scaledata,
  fontsize=10,
  color=colorRampPalette(c("blue", "yellow", "red"))(299),
  cluster_rows=TRUE,
  cluster_cols=TRUE,
  fontsize_row=1,
  fontsize_col=10,
  scale="row",
  filename=file.path(output, "Heatmap_no_cluster_plot_with_scaled_data.jpg"),
  width = 20,
  height = 25
)

pheatmap(
  scaledata,
  fontsize=10,
  color=colorRampPalette(c("blue", "yellow", "red"))(299),
  cluster_rows=TRUE,
  cluster_cols=TRUE,
  cutree_rows = k,
  clustering_distance_rows = clustering_distance,
  clustering_method = clustering_method,
  fontsize_row=1,
  fontsize_col=10,
  scale="row",
  filename=file.path(output, "Heatmap_with_cluster_plot_with_scaled_data.jpg"),
  width = 20,
  height = 25
)
cat("\n")
