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

# Set seed
set.seed(1)


#######################################################################
## Input and output (Double check or make changes)
#######################################################################
# Specify input file (csv or tab delimited text file)
input <- file.path("/home/ycth8/data/projects/HAPPI.GWAS/raw_data/Arabidopsis_1001/03_22_2020_Arabidopsis_1001_BAA.csv")

# Specify an output folder so that all the results can be stored into the folder
output <- file.path("/home/ycth8/data/projects/BioTools/output/08_06_2020")

# Start column in the dataset where the analysis should start
start_column <- 2


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
## PCA code starts here
#######################################################################

# Plot sactter plot matrix
cat("\n Save image... \n")
jpeg(file.path(output, "scatterplot_matrix.jpg"), width = 1000, height = 800)
ggscatmat(df, columns = start_column:ncol(df))
dev.off()
cat("\n")

#Performs a principal components analysis
pca_info <- prcomp(df[, start_column:ncol(df)], center = TRUE, scale. = TRUE)

# Get the standard deviations of the principal components
sdev <- data.frame(
  Trait = colnames(df)[start_column:ncol(df)],
  sdev = pca_info$sdev,
  stringsAsFactors = FALSE
)

# Get the matrix of variable loadings
rotation <- as.data.frame(pca_info$rotation, stringsAsFactors = FALSE)
rotation <- tibble::rownames_to_column(rotation, var = "Trait")

center <- as.data.frame(pca_info$center, stringsAsFactors = FALSE)
center <- tibble::rownames_to_column(center, var = "Trait")

scale <- as.data.frame(pca_info$scale, stringsAsFactors = FALSE)
scale <- tibble::rownames_to_column(scale, var = "Trait")

# Get PC
PC <- data.frame(pc1 = pca_info$x[, 1], pc2 = pca_info$x[, 2], stringsAsFactors = FALSE)

# Save all outputs to files
write.csv(sdev, file.path(output, "sdev.csv"), row.names = FALSE)
write.csv(rotation, file.path(output, "rotation.csv"), row.names = FALSE)
write.csv(center, file.path(output, "center.csv"), row.names = FALSE)
write.csv(scale, file.path(output, "scale.csv"), row.names = FALSE)
write.csv(PC, file.path(output, "PC.csv"), row.names = FALSE)
