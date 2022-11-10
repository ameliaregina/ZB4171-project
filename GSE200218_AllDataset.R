# script to perform trajectory analysis on GSE202695

setwd("D:/Monocle3/GSE200218") # set working directory and make sure it has all the files.

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
expr1 <- read_csv('GSE200218_sc_sn_integrated_data.csv', col_names = TRUE)
metadata <- read.csv('GSE200218_sc_sn_metadata.csv', header= TRUE)
markers <- read.csv('GSE200218_sc_sn_gene_names.csv', row.names = 1, header= TRUE)
counts <- read.matrix('GSE200218_sc_sn_counts.mtx', )
expr1 = fread("D:/Monocle3/GSE200218/GSE200218_sc_sn_integrated_data.csv", header = TRUE)



# create seurat object ---------------
# we need rows as genes and columns as cell_id
seu.obj1 <- CreateSeuratObject(counts = expr1)
View(seu.obj1@meta.data) #checking the metadata
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
a1 <- DimPlot(seu.obj1.filtered, reduction = 'umap', group.by = 'model', label = T)
a2 <- DimPlot(seu.obj1.filtered, reduction = 'umap', group.by = 'type', label = T)
a1|a2

# Check whether our data set has been annotated or not




# MONOCLE3 WORKFLOW ---------------------
# monocle3 requires cell_data_set object
# convert seurat object to cell_data_set object for monocle3




# ...1 Convert to cell_data_set object ------------------------

# convert seurat object into object of cell dataset class (SeuratWrapper)
cds1 <- as.cell_data_set(seu.obj1.filtered)
cds1

# to get cell metadata
colData(cds1)
# to gene/feature metadata
fData(cds1)
rownames(fData(cds1))[1:10] #rownames are the genes

# since it misses the gene_short_name column, let's add it
fData(cds1)$gene_short_name <- rownames(fData(cds1))
# i think this dataset already utilize the short naming convention for genes

# to get counts
counts(cds1)
# obtain sparse matrix of genes(Rows), Columns(Cell_ID) and numbers as the counts



# ...2. Cluster cells (using clustering info from seurat's UMAP)---------------------------
# let's use the clustering information have

# assign partitions
# is the number 1 important here
recreate.partition <- c(rep(1,length(cds1@colData@rownames)))
names(recreate.partition) <- cds1@colData@rownames
recreate.partition <- as.factor(recreate.partition)
# forms a named list which contains all the cells with value 1


cds1@clusters$UMAP$partitions <- recreate.partition

# Assign the cluster info 
# We can obtain cluster info from seurat object shown below
list_cluster <- seu.obj1.filtered@active.ident 
cds1@clusters$UMAP$clusters <- list_cluster


# Assign UMAP coordinate - cell embeddings
# UMAP is used because UMAp was what was done in seurat
cds1@int_colData@listData$reducedDims$UMAP <- seu.obj1.filtered@reductions$umap@cell.embeddings



# plot

#data visualization before trajectory
cluster.before.trajectory <- plot_cells(cds1,
                                        color_cells_by = 'model',
                                        label_groups_by_cluster = FALSE,
                                        group_label_size = 5) +
  theme(legend.position = "right")

cluster.names <- plot_cells(cds1,
                            color_cells_by = 'model',
                            label_groups_by_cluster = FALSE,
                            group_label_size = 5) +
  scale_color_manual(values = c('red', 'blue', 'green', 'maroon', 'yellow', 'grey', 'cyan')) +
  theme(legend.position = "right")

cluster.before.trajectory | cluster.names
# NOTE: for the colors we need to see how many groups are present in the dataset and adjust accordingly
# This step might be unnecessary as I can't get it to run

# ...3. Learn trajectory graph ------------------------
cds1 <- learn_graph(cds1, use_partition = FALSE)

plot_cells(cds1,
           color_cells_by = 'model',
           label_groups_by_cluster = FALSE,
           label_branch_points = FALSE,
           label_roots = FALSE,
           label_leaves = FALSE,
           group_label_size = 5)

plot_cells(cds1, 
           color_cells_by = "partition", 
           label_groups_by_cluster = FALSE,
           label_branch_points = FALSE,
           label_roots = FALSE,
           label_leaves = FALSE,)




# ...4. Order the cells in pseudotime -------------------

cds1 <- order_cells(cds1, reduction_method = 'UMAP', root_cells = colnames(cds1[,clusters(cds1) == 5]))
# Change parameter of root_cells to which cells are at the beginning of the trajectory based on
# prior knowledge (which in our case should be the original cancer line i can't remember the name)
# give the parameters all the cellID of the original/root group of cells.

plot_cells(cds1,
           color_cells_by = 'pseudotime',
           label_groups_by_cluster = FALSE,
           label_branch_points = FALSE,
           label_roots = FALSE,
           label_leaves = FALSE)

# cells ordered by monocle3 pseudotime

pseudotime(cds1)
cds1$monocle3_pseudotime <- pseudotime(cds1)
data.pseudo.GSE202695 <- as.data.frame(colData(cds1))

ggplot(data.pseudo.GSE202695, aes(monocle3_pseudotime, reorder(model, monocle3_pseudotime, median), fill = model)) +
  geom_boxplot()



# ...5. Highlighting Differentially expressed genes ------------------------


# ...6. Finding genes that change as a function of pseudotime --------------------

deg_bcells <- graph_test(cds1, neighbor_graph = 'principal_graph', cores = 4)

# graph_test test the genes on differential expression based on the low dimension embeddings
# uses test statistic to determine the different gene expression levels in cells at nearby positions on the trajectory 
# We would be interested in differentially expressed genes using the q-value column

deg_bcells %>% 
  arrange(q_value) %>% 
  filter(status == 'OK') %>% 
  head()

# This step requires us to analyze the genes 
FeaturePlot(seu.obj1.filtered, features = c('ARPC4', 'ARPC1B', 'LINC02067'))
#gene list incomplete, need to analyze data more


# visualizing pseudotime in seurat as a feature plot

seu.obj1.filtered$pseudotime <- pseudotime(cds1)
Idents(seu.obj1.filtered) <- seu.obj1.filtered$model
FeaturePlot(seu.obj1.filtered, features = "pseudotime", label = T)

# Within each cluster/over all the groups identify most differentially expressed gene and lit review.



