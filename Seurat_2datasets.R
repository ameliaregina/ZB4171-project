library(monocle3)
library(SeuratWrappers)
library(Seurat)
library(ggplot2)
library(tidyverse)
library(aws.s3)
library(aws.ec2metadata)
library(dplyr)

# 2. Seurat Installation-----------------------------------------------------------------------
install.packages('Seurat')


#3. SeuratWrapper Installation-----------------------------------------------------------------
remotes::install_github('satijalab/seurat-wrappers')
# OR
devtools::install_github('satijalab/seurat-wrappers')
# I prefer devtools but theoretically both commands work
# remotes may result in the message "API rate limit exceeded"

#4. AWS Interaction----------------------------------------------------------------------------
# Use the above ID and Secret access key, as well as the region/server information 
# to allow this Rstudio instance access to the files in the S3 Bucket
Sys.setenv("AWS_DEFAULT_REGION" = "ap-southeast-1")
# Note that for our case, Amelia has taken the liberty of providing public access to the bucket. 
# While not preferable, due to time constraints we are forced to cut corners at the moment

# Important Addresses for whoever is using this EC2 Instance
# GSE202695 Bucket: s3://zb4171-acww/GSE202695/(filename).filetype
# GSE200218 Bucket: s3://zb4171-acww/GSE200218/(filename).filetype
# GSE200218 Individual Datasets Bucket: s3://zb4171-acww/GSE200218/GSE200218_IndivData/(filename).filetype
# GSE181919 Bucket: s3://zb4171-acww/GSE181919/(filename).filetype

####################### GSE181919 #############################
expr <- aws.s3::s3read_using(read.delim, object = "s3://zb4171-acww/GSE181919/GSE181919_UMI_counts.txt") # gene metadata
metadata <- aws.s3::s3read_using(read.delim, object = "s3://zb4171-acww/GSE181919/GSE181919_Barcode_metadata.txt")

# make expr col names same as metadata row names
# use test instead of expr 
test <- expr
rownames <- rownames(metadata)
colnames(test) <- rownames

# make metadata variable row names an actual column
metadata <- tibble::rownames_to_column(metadata, "temp")

# create seurat object ---------------
# we need rows as genes and columns as cell_id
seu.obj1 <- CreateSeuratObject(counts = test)
View(seu.obj1@meta.data) #checking the metadata
colnames(metadata)[1] <- "cell_id" # fill-in Null ID
seu.obj1@meta.data <- merge(seu.obj1@meta.data, metadata, by.x = 'row.names', by.y = 'cell_id') #name in axes depends on data
View(seu.obj1@meta.data)
seu.obj1@meta.data <- seu.obj1@meta.data %>% 
  column_to_rownames(var = 'Row.names')

# percentage of reads that map to mitochondrial genes
seu.obj1$mitopercent <- PercentageFeatureSet(seu.obj1, pattern = '^MT-') 

# Filter out cells with less than 800 reads, less than 500 genes 
# and greater than 10% reads map to mito.genes
seu.obj1.filtered <- subset(seu.obj1, subset = nCount_RNA > 800 &
                              nFeature_RNA > 500 &
                              mitopercent < 10)


# subset my seurat object - (B cells) - change this to cell of interest
# this step can be skipped for our experiment
# would it be possible to skip all the steps up till MONOCLE3 Workflow?
# I assume yes OR we perform individual analyses for each subtype

# check the populations present in filtered dataset
unique(seu.obj1.filtered@meta.data$population) # returns NULL
# would "model" be better suited for this

# Check the identity (for tut this shows other identities such as BNK CD4T)
Idents(seu.obj1.filtered) # i think this returns NULL as well

# Change the identity to what is the population label (put cells into pop. group)
Idents(seu.obj1.filtered) <- seu.obj1.filtered$model
#we utilized the different model titles as a basis for their identity

# see what type of cells are present in the data
unique(seu.obj1.filtered@meta.data$model) # What else other than redefined_cluster
# check if this would be an ideal dataset to see whether there are the start, intermediate and end points
# i dont think i should have used model here, not sure what would be more appropriate

# pre-processing using seurat (this is only done for the subset of cells) can this be run for whole population?
# time-stamp 21:00 ish
seu.obj1.filtered <- NormalizeData(seu.obj1.filtered)
seu.obj1.filtered <- FindVariableFeatures(seu.obj1.filtered)
seu.obj1.filtered<- ScaleData(seu.obj1.filtered)
seu.obj1.filtered <- RunPCA(seu.obj1.filtered)
seu.obj1.filtered <- FindNeighbors(seu.obj1.filtered, dims = 1:30)
seu.obj1.filtered <- FindClusters(seu.obj1.filtered, resolution = 0.7)
# optimal resolution to avoid accidental clustering (how do you obtain optimal resolution?)
seu.obj1.filtered  <- RunUMAP(seu.obj1.filtered , dims = 1:30, n.neighbors = 50)

# How to find optimal resolution
# Iterative process, check using various resolutions and visualize as UMAP
# Can also take a look at the feature plot, expression of markers
# Repeat until you find the optimal separation of cells into distinct cell types

# Data visualization (this can be a checkpoint to see if resolutions are optimal)
a1 <- DimPlot(seu.obj1.filtered, reduction = 'umap', group.by = 'tissue.type', label = T)
a2 <- DimPlot(seu.obj1.filtered, reduction = 'umap', group.by = 'cell.type', label = T)
a1|a2


####################### GSE200218 #############################
expr_218 <- aws.s3::s3read_using(read.csv, object = "s3://zb4171-acww/GSE200218/GSE200218_sc_sn_integrated_data.csv") # gene metadata
metadata_218 <- aws.s3::s3read_using(read.csv, object = "s3://zb4171-acww/GSE200218/GSE200218_sc_sn_metadata.csv")
expr_218.with.rownames <- data.frame(expr_218[,-1], row.names=expr_218[,1])

# create seurat object ---------------
# we need rows as genes and columns as cell_id
seu.obj1 <- CreateSeuratObject(counts = expr_218.with.rownames)
View(seu.obj1@meta.data) #checking the metadata
colnames(metadata_218)[1] <- "cell_id" # fill-in Null ID
seu.obj1@meta.data <- merge(seu.obj1@meta.data, metadata_218, by.x = 'row.names', by.y = 'cell_id') #name in axes depends on data
View(seu.obj1@meta.data)
seu.obj1@meta.data <- seu.obj1@meta.data %>% 
  column_to_rownames(var = 'Row.names')

# percentage of reads that map to mitochondrial genes
seu.obj1$mitopercent <- PercentageFeatureSet(seu.obj1, pattern = '^MT-') 

# Filter out cells with less than 800 reads, less than 500 genes 
# and greater than 10% reads map to mito.genes
seu.obj1.filtered <- subset(seu.obj1, subset = nCount_RNA > 800 &
                              nFeature_RNA > 500 &
                              mitopercent < 10)


# subset my seurat object - (B cells) - change this to cell of interest
# this step can be skipped for our experiment
# would it be possible to skip all the steps up till MONOCLE3 Workflow?
# I assume yes OR we perform individual analyses for each subtype

# check the populations present in filtered dataset
unique(seu.obj1.filtered@meta.data$population) # returns NULL
# would "model" be better suited for this

# Check the identity (for tut this shows other identities such as BNK CD4T)
Idents(seu.obj1.filtered) # i think this returns NULL as well

# Change the identity to what is the population label (put cells into pop. group)
Idents(seu.obj1.filtered) <- seu.obj1.filtered$model
#we utilized the different model titles as a basis for their identity

# see what type of cells are present in the data
unique(seu.obj1.filtered@meta.data$model) # What else other than redefined_cluster
# check if this would be an ideal dataset to see whether there are the start, intermediate and end points
# i dont think i should have used model here, not sure what would be more appropriate

# pre-processing using seurat (this is only done for the subset of cells) can this be run for whole population?
# time-stamp 21:00 ish
seu.obj1.filtered <- NormalizeData(seu.obj1.filtered)
seu.obj1.filtered <- FindVariableFeatures(seu.obj1.filtered)
seu.obj1.filtered<- ScaleData(seu.obj1.filtered)
seu.obj1.filtered <- RunPCA(seu.obj1.filtered)
seu.obj1.filtered <- FindNeighbors(seu.obj1.filtered, dims = 1:30)
seu.obj1.filtered <- FindClusters(seu.obj1.filtered, resolution = 0.7)
# optimal resolution to avoid accidental clustering (how do you obtain optimal resolution?)
seu.obj1.filtered  <- RunUMAP(seu.obj1.filtered , dims = 1:30, n.neighbors = 50)

# How to find optimal resolution
# Iterative process, check using various resolutions and visualize as UMAP
# Can also take a look at the feature plot, expression of markers
# Repeat until you find the optimal separation of cells into distinct cell types

# Data visualization (this can be a checkpoint to see if resolutions are optimal)
a1 <- DimPlot(seu.obj1.filtered, reduction = 'umap', group.by = 'tissue.type', label = T)
a2 <- DimPlot(seu.obj1.filtered, reduction = 'umap', group.by = 'cell.type', label = T)
a1|a2