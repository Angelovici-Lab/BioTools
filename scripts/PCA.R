#!/usr/bin/env Rscript

#######################################################################
## Load packages
#######################################################################
rm(list = ls())

library(jpeg)
library(dplyr)
library(tibble)
library(argparse)
library(GGally)

set.seed(1)


#######################################################################
## Argument Parser
#######################################################################

parser <- argparse::ArgumentParser()

parser$add_argument("-i", type="character", help="Input file path", required=TRUE)
parser$add_argument("-o", type="character", help="Output folder path", required=TRUE)
parser$add_argument("-s", type="integer", help="Start columns", default=1)

args <- parser$parse_args()

input <- args$i
output <- args$o
start_column <- args$s

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


#######################################################################
## PCA code starts here
#######################################################################

cat("\n Save image... \n")
jpeg(file.path(output, "scatterplot_matrix.jpg"), width = 1000, height = 800)
ggscatmat(df, columns = start_column:ncol(df))
dev.off()
cat("\n")


pca_info <- prcomp(df[, start_column:ncol(df)], center = TRUE, scale. = TRUE)

sdev <- data.frame(
  Trait = colnames(df)[start_column:ncol(df)],
  sdev = pca_info$sdev,
  stringsAsFactors = FALSE
)

rotation <- as.data.frame(pca_info$rotation, stringsAsFactors = FALSE)
rotation <- tibble::rownames_to_column(rotation, var = "Trait")

center <- as.data.frame(pca_info$center, stringsAsFactors = FALSE)
center <- tibble::rownames_to_column(center, var = "Trait")

scale <- as.data.frame(pca_info$scale, stringsAsFactors = FALSE)
scale <- tibble::rownames_to_column(scale, var = "Trait")

PC <- data.frame(pc1 = pca_info$x[, 1], pc2 = pca_info$x[, 2], stringsAsFactors = FALSE)

write.csv(sdev, file.path(output, "sdev.csv"), row.names = FALSE)
write.csv(rotation, file.path(output, "rotation.csv"), row.names = FALSE)
write.csv(center, file.path(output, "center.csv"), row.names = FALSE)
write.csv(scale, file.path(output, "scale.csv"), row.names = FALSE)
write.csv(PC, file.path(output, "PC.csv"), row.names = FALSE)
