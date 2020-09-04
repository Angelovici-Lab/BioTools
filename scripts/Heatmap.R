#!/usr/bin/env Rscript

#######################################################################
## Load packages
#######################################################################
rm(list = ls())

library(pheatmap)

# Set seed
set.seed(1)


#######################################################################
## Input and output (Double check or make changes)
#######################################################################
# Specify input file (csv or tab delimited text file)
input <- file.path('/home/ycth8/data/projects/BioTools/data/Hclust.csv')

# Specify an output folder so that all the results can be stored into the folder
output <- file.path("/home/ycth8/data/projects/BioTools/output/08_14_2020")


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
print(dim(df))


#######################################################################
## Heatmap code starts here
#######################################################################

## making the first col as row.names (convert data to martix)
rnames_df <- df[,1]
matrix_df <- data.matrix(df[,2:ncol(df)])
rownames(matrix_df) <- rnames_df

# Display matrix
cat(rep("\n", 2))
print(head(matrix_df))
print(dim(matrix_df))

# Plot and save heatmap
cat("\n Save image... \n")
jpeg(file.path(output, "Heatmap.jpg"), width = 1000, height = 1200)
pheatmap(
  matrix_df,
  fontsize=10,
  color=colorRampPalette(c("blue", "yellow", "red"))(299),
  cluster_rows=TRUE,
  cluster_cols=TRUE,
  fontsize_row=1,
  fontsize_col=10,
  scale="row"
)
dev.off()
cat("\n")
