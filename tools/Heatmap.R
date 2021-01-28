#!/usr/bin/env Rscript

#######################################################################
## Load packages
#######################################################################
rm(list = ls())

library(pheatmap)


set.seed(1)


#######################################################################
## Argument Parser
#######################################################################

parser <- argparse::ArgumentParser()

parser$add_argument("-i", type="character", help="Input file path", required=TRUE)
parser$add_argument("-o", type="character", help="Output folder path", required=TRUE)

args <- parser$parse_args()

input <- args$i
output <- args$o


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
if(endsWith(input, "csv")){
  df <- read.csv(file = input, header = TRUE)
} else if (endsWith(input, "txt")){
  df <- read.delim(file = input, header=TRUE)
} else{
  print("Invalid input file. Exiting ...")
  quit(status = 0)
}


#######################################################################
##
## Heatmap code starts here
##
#######################################################################

# Show first 6 rows of the data
cat(rep("\n", 2))
print(dim(df))
print(head(df))

# Convert data into matrix format
rnames_df <- df[,1]
matrix_df <- data.matrix(df[,2:ncol(df)])
rownames(matrix_df) <- rnames_df

# Show the matrix
cat(rep("\n", 2))
print(dim(matrix_df))
print(head(matrix_df))

# Plot heatmap
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
