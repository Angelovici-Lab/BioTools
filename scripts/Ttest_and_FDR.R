#!/usr/bin/env Rscript

#######################################################################
## Load packages
#######################################################################
rm(list = ls())

library(dplyr)
library(tidyr)
library(tibble)
library(stringr)

# Set seed
set.seed(1)


#######################################################################
## Input and output (Double check or make changes)
#######################################################################
# Specify input file (csv or tab delimited text file)
input <- file.path('/data/users/ycth8/projects/BioTools/data/Proteins_Drought.csv')

# Specify an output folder so that all the results can be stored into the folder
output <- file.path("/home/ycth8/data/projects/BioTools/output/01_15_2021")

# Output files could be in csv or txt format
output_filename <- "proteins_drought.txt"
filtered_output_filename <- "proteins_drought_filtered.txt"

# Check trt column in your input file
group1 <- 2
group2 <- 6

# all records below this threshold will be kept
fdr_filter_threshold <- 0.2


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
  if(!dir.exists(output)){
    quit(status=1)
  }
}

## Load data
# the input file extension should be csv or txt
if(endsWith(input, "csv")){
  dat <- read.csv(file = input, header = TRUE)
} else if (endsWith(input, "txt")){
  dat <- read.delim(file = input, header=TRUE)
} else{
  print("Invalid input file. Exiting ...")
  quit(status = 0)
}

# Display raw data
cat("\n Raw data: \n")
print(head(dat))


#######################################################################
## T-test and false discovery rate (FDR) correction
#######################################################################
df <- dat %>%
  group_by(Accession) %>%                                                 # group by individual protein id
  summarise_at(
    .vars = "measurement",
    .funs = list(
      Group1_mean = ~mean(.[trt == group1]),                              # mean for 10% trt
      Group2_mean = ~mean(.[trt == group2]),                              # mean for 100% trt (control)
      p_value = ~t.test(.[trt == group1], .[trt == group2])$p.value)      # p-value for t-test
  ) %>%
  mutate(
    p_fdr = p.adjust(p_value, method = "fdr", n = length(p_value)),       # corrected p-value for false discovery rate
    foldchange = Group1_mean/Group2_mean,                                 # fold change (trt/control)
    log2foldchange = log2(foldchange)                                     # log2 of fold change
  ) %>%
  left_join(., dat[,c(1,2)], by = "Accession") %>%
  distinct(Accession, .keep_all = TRUE) %>%
  select(Accession, Description, Group1_mean, Group2_mean, p_value, p_fdr, foldchange, log2foldchange) %>%
  arrange(Accession) %>%
  as.data.frame(stringsAsFactors=FALSE)


#######################################################################
## Write outputs
#######################################################################
if(endsWith(output_filename, "csv")){
  write.csv(
    x=df,
    file=file.path(output, output_filename),
    row.names=FALSE,
    quote=FALSE
  )
} else if(endsWith(output_filename, "txt")){
  write.table(
    x=df,
    file=file.path(output, output_filename),
    row.names=FALSE,
    quote=FALSE
  )
}

if(endsWith(filtered_output_filename, "csv")){
  write.csv(
    x=df[df$p_fdr < fdr_filter_threshold,],
    file=file.path(output, filtered_output_filename),
    row.names=FALSE,
    quote=FALSE
  )
} else if(endsWith(filtered_output_filename, "txt")){
  write.table(
    x=df[df$p_fdr < fdr_filter_threshold,],
    file=file.path(output, filtered_output_filename),
    row.names=FALSE,
    quote=FALSE
  )
}
