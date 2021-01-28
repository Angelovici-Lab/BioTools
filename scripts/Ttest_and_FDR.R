#!/usr/bin/env Rscript

#######################################################################
## Load packages
#######################################################################
rm(list = ls())

library(dplyr)
library(tidyr)
library(tibble)
library(stringr)

library(foreach)
library(iterators)
library(parallel)
library(doParallel)

# Set seed
set.seed(1)


#######################################################################
## Input and output (Double check or make changes)
#######################################################################
# Specify input file (csv or tab delimited text file)
input <- file.path('/data/users/ycth8/projects/BioTools/data/Proteins_Drought2_Abou_Jan2019.csv')

# Specify an output folder so that all the results can be stored into the folder
output <- file.path("/home/ycth8/data/projects/BioTools/output/01_27_2021")

# Number of processing cores
cores <- ifelse(detectCores()>1, detectCores()-1, 1)

fdr_threshold <- 0.2


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
## Ttest and FDR code starts here
##
#######################################################################

#######################################################################
## Register computng cores, so tasks can be ran in parallel
#######################################################################
registerDoParallel(cores = cores)

#######################################################################
## t-test for each sample
#######################################################################
result_table <- foreach(i=3:(ncol(df)-1), .combine = rbind) %dopar% {

  subject_1_list <- c()
  subject_2_list <- c()
  mean_1_list <- c()
  mean_2_list <- c()
  p_value_list <- c()

  for(j in (i+1):ncol(df)){

    # Calculate mean
    mean_1 <- mean(df[,i], na.rm=TRUE)
    mean_2 <- mean(df[,j], na.rm=TRUE)

    # Perform t-test
    t_test_result <- tryCatch({
      t.test(df[,i], df[,j])
    }, error = function(e) {
      return(NULL)
    })

    # Collect all results
    subject_1_list <- c(subject_1_list, colnames(df)[i])
    subject_2_list <- c(subject_2_list, colnames(df)[j])
    mean_1_list <- c(mean_1_list, mean_1)
    mean_2_list <- c(mean_2_list, mean_2)
    p_value_list <- c(p_value_list, ifelse(is.null(t_test_result), NA, t_test_result$p.value))

  }

  return(
    data.frame(
      "subject_1"=subject_1_list,
      "subject_2"=subject_2_list,
      "mean_1"=mean_1_list,
      "mean_2"=mean_2_list,
      "p_value"=p_value_list,
      stringsAsFactors=FALSE
    )
  )
}

# Drop all records that t-test cannot be performed
result_table <- result_table %>% drop_na()

# Calculate FDR values
result_table$FDR <- p.adjust(result_table$p_value, method = "fdr", n = length(result_table$p_value))

# Save the t-test result
write.csv(
  x=result_table,
  file=file.path(output, "t_test_for_each_sample.csv"),
  row.names=FALSE,
  na="",
  quote=FALSE
)


#######################################################################
## Organize table
#######################################################################
df <- df %>%
  pivot_longer(!c(Accession, Description), names_to = "Group", values_to = "Measurement") %>%
  separate(Group, c("Group", "Replicate"), sep="(\\.\\s*(?=[^\\.]+$))|(_\\s*(?=[^_]+$))", extra = "drop", fill = "right") %>%
  as.data.frame(stringsAsFactors=FALSE)

print(head(df))


#######################################################################
## t-test for each group
#######################################################################

# Get unique group
unique_groups <- unique(df$Group)

subject_1_list <- c()
subject_2_list <- c()
mean_1_list <- c()
mean_2_list <- c()
p_value_list <- c()

result_table <- foreach(i=1:(length(unique_groups)-1), .combine = rbind) %dopar% {
  for(j in (i+1):length(unique_groups)){

    #  Calculate mean
    mean_1 <- mean(df[df$Group == unique_groups[i],"Measurement"], na.rm=TRUE)
    mean_2 <- mean(df[df$Group == unique_groups[j],"Measurement"], na.rm=TRUE)

    # Perform t-test
    t_test_result <- tryCatch({
      t.test(
        df[df$Group == unique_groups[i],"Measurement"],
        df[df$Group == unique_groups[j],"Measurement"]
      )
    }, error = function(e) {
      return(NULL)
    })

    # Collect all results
    subject_1_list <- c(subject_1_list, unique_groups[i])
    subject_2_list <- c(subject_2_list, unique_groups[j])
    mean_1_list <- c(mean_1_list, mean_1)
    mean_2_list <- c(mean_2_list, mean_2)
    p_value_list <- c(p_value_list, ifelse(is.null(t_test_result), NA, t_test_result$p.value))

  }

  return(
    data.frame(
      "subject_1"=subject_1_list,
      "subject_2"=subject_2_list,
      "mean_1"=mean_1_list,
      "mean_2"=mean_2_list,
      "p_value"=p_value_list,
      stringsAsFactors=FALSE
    )
  )
}

# Drop all records that t-test cannot be performed
result_table <- result_table %>% drop_na()

# Calculate FDR values
result_table$FDR <- p.adjust(result_table$p_value, method = "fdr", n = length(result_table$p_value))

# Save the t-test result
write.csv(
  x=result_table,
  file=file.path(output, "t_test_for_each_group.csv"),
  row.names=FALSE,
  na="",
  quote=FALSE
)


#######################################################################
## t-test for each accession
#######################################################################

# Get unique accession
unique_accession <- unique(df$Accession)

# For each accession, do t-test between group
result_table <- foreach(i=1:length(unique_accession), .combine = rbind) %dopar% {
  dat <- df[df$Accession == unique_accession[i],]

  # Get unique group
  unique_groups <- unique(dat$Group)

  subject_1_list <- c()
  subject_2_list <- c()
  mean_1_list <- c()
  mean_2_list <- c()
  p_value_list <- c()

  if(nrow(dat) > 0){

    for(j in 1:(length(unique_groups)-1)){
      for(k in (j+1):length(unique_groups)){

        #  Calculate mean
        mean_1 <- mean(dat[dat$Group == unique_groups[j],"Measurement"], na.rm=TRUE)
        mean_2 <- mean(dat[dat$Group == unique_groups[k],"Measurement"], na.rm=TRUE)

        # Perform t-test (if all replicates have same value, t-test cannot be performed, so return NULL instead)
        t_test_result <- tryCatch({
          t.test(
            dat[dat$Group == unique_groups[j],"Measurement"],
            dat[dat$Group == unique_groups[k],"Measurement"]
          )
        }, error = function(e) {
          return(NULL)
        })

        # Collect all results
        subject_1_list <- c(subject_1_list, unique_groups[j])
        subject_2_list <- c(subject_2_list, unique_groups[k])
        mean_1_list <- c(mean_1_list, mean_1)
        mean_2_list <- c(mean_2_list, mean_2)
        p_value_list <- c(p_value_list, ifelse(is.null(t_test_result), NA, t_test_result$p.value))
      }
    }

    result_df <- data.frame(
      "Accession"=unique_accession[i],
      "Description"=dat$Description[1],
      "subject_1"=subject_1_list,
      "subject_2"=subject_2_list,
      "mean_1"=mean_1_list,
      "mean_2"=mean_2_list,
      "p_value"=p_value_list,
      stringsAsFactors=FALSE
    )
  } else{
    result_df <- data.frame(
      "Accession"=unique_accession[i],
      "Description"="",
      "subject_1"="",
      "subject_2"="",
      "mean_1"=NA,
      "mean_2"=NA,
      "p_value"=NA,
      stringsAsFactors=FALSE
    )
  }

  # Return the result
  return(result_df)
}

# Drop all records that t-test cannot be performed
result_table <- result_table %>% drop_na()

# Calculate FDR values
result_table$FDR <- p.adjust(result_table$p_value, method = "fdr", n = length(result_table$p_value))

# Save the t-test result
write.csv(
  x=result_table,
  file=file.path(output, "t_test_for_each_accession.csv"),
  row.names=FALSE,
  na="",
  quote=FALSE
)

# Save the t-test result (significant only)
write.csv(
  x=result_table[result_table$FDR <= fdr_threshold,],
  file=file.path(output, "t_test_for_each_accession_significant_only.csv"),
  row.names=FALSE,
  na="",
  quote=FALSE
)


#######################################################################
## free the processing cores
#######################################################################
stopImplicitCluster()
