# PLease ensure you have all the following packages installed before starting your analysis

# Library list
library(monocle3)
library(SeuratWrappers)
library(Seurat)
library(ggplot2) #misc
library(tidyverse) #misc
library(data.table) #misc
library(readr) #misc
library(aws.s3) #for EC2&S3
library(aws.ec2metadata) #for EC2&S3

# 1. Monocle3 Installation | adapted from: Monocle 3 Installation tutorial:--------------------
# https://cole-trapnell-lab.github.io/monocle3/docs/installation/
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.15") # for latest version of R
# previous iteration can utilize 3.14. This is the bare minimum for monocle3

BiocManager::install(c('BiocGenerics', 'DelayedArray', 'DelayedMatrixStats',
                       'limma', 'lme4', 'S4Vectors', 'SingleCellExperiment',
                       'SummarizedExperiment', 'batchelor', 'Matrix.utils',
                       'HDF5Array', 'terra', 'ggrastr'))

install.packages("devtools")
devtools::install_github('cole-trapnell-lab/monocle3')
#install_github("lme4/lme4",dependencies=TRUE)


# Note to reader. When downloading Seurat and SeuratWrapper you may encounter a message asking
# if you want to download files that need additional compiling. Always say yes
# if these packages have a non-zero exit status, utilize devtools 
# and download straight from repository 
# once you do that restart the installation and it should work fine


# 2. Seurat Installation-----------------------------------------------------------------------
install.packages('Seurat')


#3. SeuratWrapper Installation-----------------------------------------------------------------
remotes::install_github('satijalab/seurat-wrappers')
# OR
devtools::install_github('satijalab/seurat-wrappers')
# I prefer devtools but theoretically both commands work
# remotes may result in the message "API rate limit exceeded"

#4. AWS Interaction----------------------------------------------------------------------------
# Access ID: 
# Secret access key: 
# Use your access ID and Secret access key, as well as the region/server information 
# to allow this Rstudio instance access to the files in the S3 Bucket
Sys.setenv("AWS_ACCESS_KEY_ID" = "...",
           "AWS_SECRET_ACCESS_KEY" = "...",
           "AWS_DEFAULT_REGION" = "ap-southeast-1")
# No access ID and access key would be needed if the bucket can be accessed publicly


# below is an example of how to read the files for the Seurat/Monocle3 Workflow via RStudio Image 
data <- aws.s3::s3read_using(read.csv, object = "s3://bucketname/GSE202695_metadata.csv")
