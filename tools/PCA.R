#!/usr/bin/env Rscript

#######################################################################
## Load packages
#######################################################################
rm(list = ls())

library(jpeg)
library(dplyr)
library(tidyr)
library(tibble)
library(ggplot2)
library(argparse)
library(GGally)
library(ggfortify)

# Set seed
set.seed(1)


#######################################################################
## Argument Parser
#######################################################################

parser <- argparse::ArgumentParser()

parser$add_argument("-i", type="character", help="Input file path", required=TRUE)
parser$add_argument("-o", type="character", help="Output folder path", required=TRUE)
parser$add_argument("-s", type="integer", help="Start column index", required=TRUE)

args <- parser$parse_args()

input <- args$i
output <- args$o
start_column <- args$s


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
## PCA code starts here
##
#######################################################################

#######################################################################
## PCA do not allow NA, -Inf, and Inf
#######################################################################
# Drop all NAs
df <- df %>% drop_na() %>% as.data.frame(stringsAsFactors=FALSE)

# Drop rows with infinite sum
df <- df[is.finite(rowSums(df[, start_column:ncol(df)])),]


#######################################################################
## Plot sactter plot matrix
#######################################################################
cat("\n Save image... \n")
p <- ggscatmat(df, columns = start_column:ncol(df))
ggsave(
  filename = "scatterplot_matrix.jpg",
  plot = p,
  path = file.path(output),
  width = 10,
  height = 10,
  units = "in",
  dpi = 300
)
cat("\n")


#######################################################################
## Organize data for PCA
#######################################################################
pca_df <- df[, c(1, start_column:ncol(df))]

rownames(pca_df) <- pca_df[,1]
pca_df <- pca_df[,-1]

pca_df <- t(pca_df) %>%
  as.data.frame(stringsAsFactors=TRUE)

cat(rep("\n", 3))
if(ncol(pca_df)>10){
  print(head(pca_df[,1:10]))
} else{
  print(head(pca_df))
}


#######################################################################
## Performs a principal components analysis
#######################################################################
pca <- prcomp(pca_df, center = TRUE, scale. = TRUE)

pca_variance <- pca$sdev^2
pca_variance_percentage <- round(pca_variance/sum(pca_variance)*100, digits=2)

x <- pca$x %>% as.data.frame(stringsAsFactors=TRUE)


#######################################################################
## Plot PC1 vs PC2 with GGplot2
#######################################################################
p <- ggplot(x, aes(x=PC1, y=PC2)) +
  geom_text(label=rownames(x)) +
  labs(
    x=paste0("PC1 - ", pca_variance_percentage[1], "%"),
    y=paste0("PC2 - ", pca_variance_percentage[2], "%")
  )

ggsave(
  filename="PCA_Plot.jpg",
  plot=p,
  path=output
)


#######################################################################
## Plot PC1 vs PC2 with AutoPlot
#######################################################################
p <- autoplot(pca, data=pca_df, label=TRUE)

ggsave(
  filename="PCA_Autolot.jpg",
  plot=p,
  path=output
)


#######################################################################
## Get all PCA information
#######################################################################
sdev <- pca$sdev %>% as.data.frame(stringsAsFactors = FALSE)
rotation <- pca$rotation %>% as.data.frame(stringsAsFactors = FALSE)
center <- pca$center %>% as.data.frame(stringsAsFactors = FALSE)
scale <- pca$scale %>% as.data.frame(stringsAsFactors = FALSE)

# Get PC
PC <- pca$x %>% as.data.frame(stringsAsFactors = FALSE)

# # Save all outputs to files
write.csv(sdev, file.path(output, "sdev.csv"), row.names = FALSE)
write.csv(rotation, file.path(output, "rotation.csv"), row.names = TRUE)
write.csv(center, file.path(output, "center.csv"), row.names = FALSE)
write.csv(scale, file.path(output, "scale.csv"), row.names = FALSE)
write.csv(PC, file.path(output, "PC.csv"), row.names = TRUE)
