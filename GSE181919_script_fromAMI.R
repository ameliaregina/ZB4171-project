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
BiocManager::install(version = "3.15") #for latest version of R
# previous iteration can use 3.14. This is the bare minimum for monocle3

BiocManager::install(c('BiocGenerics', 'DelayedArray', 'DelayedMatrixStats',
                       'limma', 'lme4', 'S4Vectors', 'SingleCellExperiment',
                       'SummarizedExperiment', 'batchelor', 'Matrix.utils',
                       'HDF5Array', 'terra', 'ggrastr'))

install.packages('devtools') #fuck u devtools
devtools::install_github('cole-trapnell-lab/monocle3')
#install_github("lme4/lme4",dependencies=TRUE)

#wei jing, [11/8/2022 3:36 PM]
#> library("devtools")
#> install_github("lme4/lme4",dependencies=TRUE)

#wei jing, [11/8/2022 3:36 PM]
#> packageurl<-"https://cran.r-project.org/src/contrib/Archive/nloptr/nloptr_1.2.1.tar.gz"
#> 
#> install.packages(packageurl, repos=NULL, type="source")
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


#4. AWS S3 Installation------------------------------------------------------------------------
# for aws.s3
install.packages("aws.s3", repos = c("cloudyr" = "http://cloudyr.github.io/drat"))

# for aws.ec2metadata
install.packages("aws.ec2metadata", repos = c(cloudyr = "http://cloudyr.github.io/drat", getOption("repos")))

#5. AWS Interaction----------------------------------------------------------------------------
# Access ID: AKIAR4S22DSAKPBWY5ZA 
# Secret access key: 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 
# Use the above ID and Secret access key, as well as the region/server information 
# to allow this Rstudio instance access to the files in the S3 Bucket
Sys.setenv("AWS_ACCESS_KEY_ID" = "AKIAR4S22DSAKPBWY5ZA ",
           "AWS_SECRET_ACCESS_KEY" = "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",
           "AWS_DEFAULT_REGION" = "ap-southeast-1")
# Note that for our case, Amelia has taken the liberty of providing public access to the bucket. 
# While not preferable, due to time constraints we are forced to cut corners at the moment

# Important Addresses for whoever is using this EC2 Instance
# GSE202695 Bucket: s3://zb4171-acww/GSE202695/(filename).filetype
# GSE200218 Bucket: s3://zb4171-acww/GSE200218/(filename).filetype
# GSE200218 Individual Datasets Bucket: s3://zb4171-acww/GSE200218/GSE200218_IndivData/(filename).filetype
# GSE181919 Bucket: s3://zb4171-acww/GSE181919/(filename).filetype

# below is an example of how to read the files for the Seurat/Monocle3 Workflow via RStudio Image 
data <- aws.s3::s3read_using(read.csv, object = "s3://zb4171-acww/GSE202695_metadata.csv")



####################### GSE181919 #############################
expr <- aws.s3::s3read_using(read.delim, object = "s3://zb4171-acww/GSE181919/GSE181919_UMI_counts.txt") # gene metadata
metadata <- aws.s3::s3read_using(read.delim, object = "s3://zb4171-acww/GSE181919/GSE181919_Barcode_metadata.txt")

# make expr col names same as metadata row names
# use test instead of expr 
test <- expr
rownames <- rownames(metadata)
colnames(test) <- rownames

## save removed rows from metadata
metadata_new <- metadata
removed <- rownames(metadata_new[grepl("NL", metadata_new$tissue.type),])
removed_LP <- rownames(metadata_new[grepl("LP", metadata_new$tissue.type),])
removed_final <- append (removed, removed_LP)

## remove metadata rows with "sn"
metadata_remove_nl <- metadata_new[!grepl("NL", metadata_new$tissue.type),]
metadata_remove_nl_lp <- metadata_remove_nl[!grepl("LP", metadata_remove_nl$tissue.type),]

# make metadata variable row names an actual column
metadata_remove_nl_lp <- tibble::rownames_to_column(metadata_remove_nl_lp, "temp")

## remove same rows in expr data
test <- test[,!(names(test) %in% removed_final)]

# create seurat object ---------------
# we need rows as genes and columns as cell_id
seu.obj1 <- CreateSeuratObject(counts = test)
View(seu.obj1@meta.data) #checking the metadata
colnames(metadata_remove_nl_lp)[1] <- "cell_id" # fill-in Null ID
seu.obj1@meta.data <- merge(seu.obj1@meta.data, metadata_remove_nl_lp, by.x = 'row.names', by.y = 'cell_id') #name in axes depends on data
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

# MONOCLE3 WORKFLOW ---------------------
# monocle3 requires cell_data_set object
# convert seurat object to cell_data_set object for monocle3




# ...1 Convert to cell_data_set object ------------------------

# convert seurat object into object of cell dataset class (SeuratWrapper)
cds <- as.cell_data_set(seu.obj1.filtered)
cds

# to get cell metadata
colData(cds)
# to gene/feature metadata
fData(cds)
rownames(fData(cds))[1:10] #rownames are the genes

# since it misses the gene_short_name column, let's add it
fData(cds)$gene_short_name <- rownames(fData(cds))
# i think this dataset already utilize the short naming convention for genes

# to get counts
counts(cds)
# obtain sparse matrix of genes(Rows), Columns(Cell_ID) and numbers as the counts



# ...2. Cluster cells (using clustering info from seurat's UMAP)---------------------------
# let's use the clustering information have

# assign partitions
# is the number 1 important here
recreate.partition <- c(rep(1,length(cds@colData@rownames)))
names(recreate.partition) <- cds@colData@rownames
recreate.partition <- as.factor(recreate.partition)
# forms a named list which contains all the cells with value 1


cds@clusters$UMAP$partitions <- recreate.partition

# Assign the cluster info 
# We can obtain cluster info from seurat object shown below
list_cluster <- seu.obj1.filtered@active.ident 
cds@clusters$UMAP$clusters <- list_cluster


# Assign UMAP coordinate - cell embeddings
# UMAP is used because UMAp was what was done in seurat
cds@int_colData@listData$reducedDims$UMAP <- seu.obj1.filtered@reductions$umap@cell.embeddings



# plot

#data visualization before trajectory
cluster.before.trajectory <- plot_cells(cds,
                                        color_cells_by = 'tissue.type',
                                        label_groups_by_cluster = FALSE,
                                        group_label_size = 5,show_trajectory_graph = FALSE) +
  theme(legend.position = "right")

cluster.names <- plot_cells(cds,
                            color_cells_by = 'tissue.type',
                            label_groups_by_cluster = FALSE,
                            group_label_size = 5, show_trajectory_graph = FALSE) +
  scale_color_manual(values = c('pink', 'blue', 'yellow', 'cyan', 'purple', 'grey', 'purple')) +
  theme(legend.position = "right")

cluster.before.trajectory | cluster.names
# NOTE: for the colors we need to see how many groups are present in the dataset and adjust accordingly
# This step might be unnecessary as I can't get it to run

# ...4. Highlighting Differentially expressed genes ------------------------

# This contains different metrics for how specifically expressed each gene is in each partition
marker_test_res <- top_markers(cds, group_cells_by="tissue.type", 
                               reference_cells=1000, cores=8)
marker_test_res <- top_markers(cds, group_cells_by="cell.type", 
                               reference_cells=1000, cores=8)



# Pseudo_R2 is a metric: other metrics include: (to add)
top_specific_markers <- marker_test_res %>%
  filter(fraction_expressing >= 0.10) %>%
  group_by(cell_group) %>%
  top_n(1, pseudo_R2)

top_specific_marker_ids <- unique(top_specific_markers %>% pull(gene_id))


# Plot expression and fraction of cells that express each marker in each group
plot_genes_by_group(cds,
                    top_specific_marker_ids,
                    group_cells_by="tissue.type",
                    ordering_type="maximal_on_diag",
                    max.size=10)


# Vary size and visualize (top_n(x, Metric) where x is the size)
top_specific_markers <- marker_test_res %>%
  filter(fraction_expressing >= 0.10) %>%
  group_by(cell_group) %>%
  top_n(3, pseudo_R2)

top_specific_marker_ids <- unique(top_specific_markers %>% pull(gene_id))

plot_genes_by_group(cds,
                    top_specific_marker_ids,
                    group_cells_by="tissue.type",
                    ordering_type="cluster_row_col",
                    max.size=8)

# Monocle 3 allows for isolation of groups and clusters observed----------
cds_subset <- choose_cells(cds) # will lead to a pop-up


## for each subset of cell type subsetted, run code below to get dot plot
top_specific_markers <- marker_test_res %>%
  filter(fraction_expressing >= 0.10) %>%
  group_by(cell_group) %>%
  top_n(11, pseudo_R2)

top_specific_marker_ids <- unique(top_specific_markers %>% pull(gene_id))

plot_genes_by_group(cds_subset,
                    top_specific_marker_ids,
                    group_cells_by="tissue.type",
                    ordering_type="cluster_row_col",
                    max.size=10)


################# error in below code; use the one on top instead #######################
pr_graph_test_res <- graph_test(cds_subset, neighbor_graph="knn", cores=8)
pr_deg_ids <- row.names(subset(pr_graph_test_res, morans_I > 0.01 & q_value < 0.05))

gene_module_df <- find_gene_modules(cds_subset[pr_deg_ids,], resolution=1e-3)

plot_cells(cds_subset, genes=gene_module_df, 
           show_trajectory_graph=FALSE, 
           label_cell_groups=FALSE)
##########################################################################################


# While by-right we do not have a clear idea of which identity each cluster is, 
# i think it would be safe for us to assume based on gene expression
cds_subset <- cluster_cells(cds_subset, resolution=1e-2)
plot_cells(cds_subset, color_cells_by="cluster")


# ...3. Learn trajectory graph ------------------------
cds1 <- learn_graph(cds, use_partition = FALSE)

plot_cells(cds1,
           color_cells_by = 'tissue.type',
           label_groups_by_cluster = FALSE,
           label_branch_points = FALSE,
           label_roots = FALSE,
           label_leaves = FALSE,
           group_label_size = 5)

plot_cells(cds1, 
           color_cells_by = "tissue.type", 
           label_groups_by_cluster = FALSE,
           label_branch_points = FALSE,
           label_roots = FALSE,
           label_leaves = FALSE,)




# ...5. Order the cells in pseudotime -------------------

cds1 <- order_cells(cds1, reduction_method = 'UMAP', root_cells = colnames(cds1[,clusters(cds1) == 7]))
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
data.pseudo.GSE181919 <- as.data.frame(colData(cds1))

ggplot(data.pseudo.GSE181919, aes(monocle3_pseudotime, reorder(cell.type, monocle3_pseudotime, median), fill = cell.type)) +
  geom_boxplot()


# Another method to subset cells based on the branching you see is
cds_sub <- choose_graph_segments(cds)
# However this would be done after trajectory generation
# And I am not sure if this can be done for the 3D trajectory

# ...6. Working with 3D trajectories---------------------------------------------
# God I fucking hope this works
cds_3d <- reduce_dimension(cds, max_components = 3)
cds_3d <- cluster_cells(cds_3d)
cds_3d <- learn_graph(cds_3d)

## got error for this line
## just "cds_3d <- order_cells(cds_3d)" and choose root node yourself
cds_3d <- order_cells(cds_3d, root_pr_nodes=get_earliest_principal_node(cds))

cds_3d_plot_obj <- plot_cells_3d(cds_3d, color_cells_by="tissue.type")


# ...7. Finding genes that change as a function of pseudotime --------------------

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

