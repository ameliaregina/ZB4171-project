setwd("D:/Monocle3/GSE181919") # set working directory and make sure it has all the files.

set.seed(1234) # lock seed to this value to ensure reproducibility

library(monocle3)
library(SeuratWrappers)
library(Seurat)
library(ggplot2)
library(tidyverse)
library(data.table)
library(readr)


# read in data [To be changed depending on our own dataset]
# markers <- read.delim('ABC_Marker.txt', header = T) # gene metadata. [This seems to be missing from our dataset]
#metadata.trial1 <- read.delim('GSE202695_metadata.csv', header = T, sep = ',') # cell metadata [This exists as csv and with minimal information]
#expr1 <- read.delim('GSE202695_counts_afterQC.csv', header = T, sep = ',') # expression matrix 
# for expression matrix which of the three do I utilize?

#actual code run
expr <- read_delim('GSE181919_UMI_counts.txt', col_names = T) # gene metadata
metadata <- read.delim('GSE181919_Barcode_metadata.txt', header = T)
expr <- fread('GSE181919_UMI_counts.txt', header = T)
expr1 <- fread('test.csv')