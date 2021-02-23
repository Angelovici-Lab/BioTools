#!/usr/bin/env Rscript

#######################################################################
## Load packages
#######################################################################
rm(list = ls())

library(dplyr)
library(tidyr)
library(tibble)
library(ggplot2)
library(lme4)
library(nlme)
library(argparse)

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
## Heritability code starts here
##
#######################################################################

#######################################################################
## Heritability model
#######################################################################
# Initialize empty arrays
trait_arr <- c()
heritability1_arr <- c()
heritability2_arr <- c()

for(i in start_column:ncol(df)){
  # Get trait
  trait_arr <- c(trait_arr, colnames(df)[i])

  # Subset columns
  selected_df <- df[,c(1,i)] %>%
    drop_na() %>%
    as.data.frame(stringsAsFactors = FALSE)

  # Convert column 1 to factor
  selected_df[,1] <- as.factor(selected_df[,1])

  # Fit a model
  frml <- paste0(colnames(selected_df)[2], " ~ ", "(1|", colnames(selected_df)[1], ")")
  mdl <- lmer(as.formula(frml), data=selected_df)

  # Extract variance and correlation components
  vars <- VarCorr(mdl)

  # Calculate heritability
  sigmas <- as.data.frame(vars, stringsAsFactors = FALSE)$vcov
  H2 <- sigmas[1] / sum(sigmas)

  # Get heritability1
  heritability1_arr <- c(heritability1_arr, H2)

  # Get geno and pheno
  geno <- selected_df[,1]
  pheno <- selected_df[,2]

  # Fit a model then extract variance and correlation components
  mdl <- lme(pheno ~ 1 , random= ~1 | geno)
  vars <- VarCorr(mdl)

  # Calculate heritability
  va<-as.numeric(vars[1])
  vw<-as.numeric(vars[2])
  H2 <- va/(va+vw)

  # Get heritability2
  heritability2_arr <- c(heritability2_arr, H2)
}


#######################################################################
## Output trait, heritability1, and heritability2 into a data frame
#######################################################################
heritability_df <- data.frame(
  "Trait" = trait_arr,
  "Heritability1" = heritability1_arr,
  "Heritability2" = heritability2_arr,
  stringsAsFactors = FALSE
)


#######################################################################
## Output heritability_df
#######################################################################
write.csv(
  x=heritability_df,
  file=file.path(output, "heritability_table.csv"),
  row.names=FALSE,
  na="",
  quote=FALSE
)
