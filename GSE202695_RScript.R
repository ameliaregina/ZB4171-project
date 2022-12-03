# script to perform trajectory analysis on GSE202695 and its subsets

setwd("...") # set working directory and make sure it has all the files.

set.seed(1234) # lock seed to this value to ensure reproducibility

# Some features in this code will require the use of additional packages
# Please install them when prompted
library(monocle3)
library(SeuratWrappers)
library(Seurat)
library(ggplot2)
library(tidyverse)


# read in data [Varies depending on the data format]

expr1 <- read.csv('GSE202695_counts_afterQC.csv', row.names = 1, header= TRUE) # expression matrix
metadata <- read.csv('GSE202695_metadata.csv', header= TRUE) # cell metadata

# no gene metadata for this particular dataset

# for raw data of expression matrix[if available]
# if you plan on performing the QC with different parameters use this code below
expr2 <- read.csv('GSE202695_counts_raw.csv', row.names = 1, header= TRUE)

# verify how content is organized for each data and modify code accordingly

# create seurat object ---------------
# we need rows as genes and columns as cell_id
seu.obj1 <- CreateSeuratObject(counts = expr1)
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


# subset my seurat object -----------------------------------

# Check whether our data set has been annotated or not

# check the populations present in filtered dataset
unique(seu.obj1.filtered@meta.data$model)
unique(seu.obj1.filtered@meta.data$type)
# change "model" to the title present in your dataset

# Check the identity (This would shows other identities such as BNK CD4T)
Idents(seu.obj1.filtered)

# Change the identity to what is the population label (put cells into pop. group)
Idents(seu.obj1.filtered) <- seu.obj1.filtered$model
# GSE202695 Identities: HBRX2353 MDAMB231 HBRX3078 HBRX1921 HBRX2344

# You can opt to skip this step and perform the pipeline using the whole population
# Some of the names used downstream will need to be adjusted accordingly

# subset to only members in specific population
MDAMB231.seu <- subset(seu.obj1.filtered, idents = "MDAMB231")
MDAMB231.seu

# see what type of cells are present in the data
unique(MDAMB231.seu@meta.data$type)
# check if this would be an ideal dataset to see whether there are the start, intermediate and end points

# pre-processing using seurat (this is only done for the subset of cells) can this be run for whole population?
# time-stamp 21:00 ish
MDAMB231.seu <- NormalizeData(MDAMB231.seu)
MDAMB231.seu <- FindVariableFeatures(MDAMB231.seu)
MDAMB231.seu <- ScaleData(MDAMB231.seu)
MDAMB231.seu <- RunPCA(MDAMB231.seu)
MDAMB231.seu <- FindNeighbors(MDAMB231.seu, dims = 1:30)
MDAMB231.seu <- FindClusters(MDAMB231.seu, resolution = 0.9)
# optimal resolution to avoid accidental clustering (how do you obtain optimal resolution?)
MDAMB231.seu <- RunUMAP(MDAMB231.seu, dims = 1:30, n.neighbors = 30)

# How to find optimal resolution
# Iterative process, check using various resolutions and visualize as UMAP
# Can also look at the feature plot, expression of markers
# Repeat until you find the optimal separation of cells into distinct cell types


# Data visualization (this can be a checkpoint to see if resolutions are optimal)
a1 <- DimPlot(MDAMB231.seu, reduction = 'umap', group.by = 'model', label = T)
a2 <- DimPlot(MDAMB231.seu, reduction = 'umap', group.by = 'type', label = T)
a1|a2



# MONOCLE3 WORKFLOW ---------------------
# monocle3 requires cell_data_set object
# convert seurat object to cell_data_set object for monocle3

# ...1 Convert to cell_data_set object ------------------------

# convert seurat object into object of cell dataset class (SeuratWrapper)
cds1 <- as.cell_data_set(MDAMB231.seu)
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
# is the number 1 important here?
recreate.partition <- c(rep(1,length(cds1@colData@rownames)))
names(recreate.partition) <- cds1@colData@rownames
recreate.partition <- as.factor(recreate.partition)
# forms a named list which contains all the cells with value 1


cds1@clusters$UMAP$partitions <- recreate.partition

# Assign the cluster info 
# We can obtain cluster info from seurat object shown below
list_cluster <- MDAMB231.seu@active.ident 
cds1@clusters$UMAP$clusters <- list_cluster


# Assign UMAP coordinate - cell embeddings
# UMAP is used because UMAP was what was done in seurat
# and many downstream Monocle 3 functionality requires UMAP
cds1@int_colData@listData$reducedDims$UMAP <- MDAMB231.seu@reductions$umap@cell.embeddings



# plot 
#data visualization before trajectory
cluster.before.trajectory <- plot_cells(cds1,
                                        color_cells_by = 'cluster',
                                        label_groups_by_cluster = TRUE,
                                        group_label_size = 5,
                                        cell_size = 1.5) +
  theme(legend.position = "right")

cluster.names <- plot_cells(cds1,
                            color_cells_by = 'type',
                            label_groups_by_cluster = FALSE,
                            group_label_size = 5,
                            cell_size = 1.5) +
  scale_color_manual(values = c('red', 'blue', 'green', 'maroon', 'yellow', 'grey', 'cyan')) +
  theme(legend.position = "right")

cluster.before.trajectory | cluster.names
# NOTE: for the colors we need to see how many groups are present in the dataset and adjust accordingly
# you can adjust and create other ways to cluster (e.g. by cell type). Use this as a checkpoint of sorts


# ...3. Highlighting marker genes for each cluster ------------------------

# This contains different metrics for how specifically expressed each gene is in each partition
marker_test_res <- top_markers(cds1, group_cells_by="type", 
                               reference_cells=1000, cores=8)
marker_test_res <- top_markers(cds1, group_cells_by="partition", 
                               reference_cells=1000, cores=8)
marker_test_res <- top_markers(cds1, group_cells_by="cluster", 
                               reference_cells=1000, cores=8)



# Pseudo_R2 is a metric. Other metrics are described in the monocle3 documentation
top_specific_markers <- marker_test_res %>%
  filter(fraction_expressing >= 0.10) %>%
  group_by(cell_group) %>%
  top_n(1, pseudo_R2)

top_specific_marker_ids <- unique(top_specific_markers %>% pull(gene_id))


# Plot expression and fraction of cells that express each marker in each group
plot_genes_by_group(cds1,
                    top_specific_marker_ids,
                    group_cells_by="type",
                    ordering_type="maximal_on_diag",
                    max.size=10)


# Vary size and visualize (top_n(x, Metric) where x is the size)
top_specific_markers <- marker_test_res %>%
  filter(fraction_expressing >= 0.10) %>%
  group_by(cell_group) %>%
  top_n(7, pseudo_R2)

top_specific_marker_ids <- unique(top_specific_markers %>% pull(gene_id))

plot_genes_by_group(cds1,
                    top_specific_marker_ids,
                    group_cells_by="type",
                    ordering_type="cluster_row_col",
                    max.size=10)


# Monocle 3 allows for isolation of groups and clusters observed----------
# Unless necessary, additional clustering is not needed for Subset analysis script
cds_subset <- choose_cells(cds1) # will lead to a pop-up

# Start of error - I can't seem to run this ------------------------------
# is it due to the missing cell metadata file? 
# Note to Chris: Try running a separate script with only Monocle 3
pr_graph_test_res <- graph_test(cds_subset, neighbor_graph="knn", cores=8)
pr_deg_ids <- row.names(subset(pr_graph_test_res, morans_I > 0.01 & q_value < 0.05))

gene_module_df <- find_gene_modules(cds_subset, resolution=1e-3)

plot_cells(cds_subset, genes=gene_module_df, 
           show_trajectory_graph=FALSE, 
           label_cell_groups=FALSE)
# End of Error -----------------------------------------------------------

# While by-right we do not have a clear idea of which identity each cluster is, 
# I think it would be safe for us to assume based on gene expression
cds_subset <- cluster_cells(cds_subset, resolution=1e-2)
c1.1 <- plot_cells(cds_subset, color_cells_by="cluster", group_label_size = 5)
c1.2 <- plot_cells(cds_subset, color_cells_by="type", group_label_size = 5)

c1.1|c1.2


# ...4.1.1. Learn trajectory graph ------------------------
# for whole population [or for the individual subset analysis]
cds1 <- learn_graph(cds1, use_partition = FALSE)
cds1 <- learn_graph(cds1, use_partition = TRUE)

plot_cells(cds1,
           color_cells_by = 'model',
           label_groups_by_cluster = FALSE,
           label_branch_points = TRUE,
           label_roots = TRUE,
           label_leaves = TRUE,
           group_label_size = 5,
           cell_size = 1.5)+
  theme(legend.position = "right")

plot_cells(cds1,
           color_cells_by = 'type',
           label_groups_by_cluster = FALSE,
           label_branch_points = TRUE,
           label_roots = TRUE,
           label_leaves = TRUE,
           group_label_size = 5,
           cell_size = 1.5)+
  theme(legend.position = "right")

plot_cells(cds1, 
           color_cells_by = "cluster", 
           label_groups_by_cluster = TRUE,
           label_branch_points = TRUE,
           label_roots = TRUE,
           label_leaves = TRUE,
           group_label_size = 5,
           cell_size = 1.5)+
  theme(legend.position = "right")

plot_cells(cds1, 
           color_cells_by = "partition", 
           label_groups_by_cluster = FALSE,
           label_branch_points = TRUE,
           label_roots = TRUE,
           label_leaves = TRUE,
           group_label_size = 5,
           cell_size = 1.5)+
  theme(legend.position = "right")

# ...4.1.2. Learn trajectory graph ------------------------
# for subset if performed in previous step
cds_subset <- learn_graph(cds_subset, use_partition = FALSE)

c1.1.plot <- plot_cells(cds_subset,
                        color_cells_by = 'cluster',
                        label_groups_by_cluster = TRUE,
                        label_branch_points = TRUE,
                        label_roots = FALSE,
                        label_leaves = FALSE,
                        group_label_size = 5)+
  theme(legend.position = "right")

c1.2.plot <- plot_cells(cds_subset, 
                        color_cells_by = "type", 
                        label_groups_by_cluster = FALSE,
                        label_branch_points = TRUE,
                        label_roots = FALSE,
                        label_leaves = FALSE,)+
  theme(legend.position = "right")

# to check against whole population
c.plot <- plot_cells(cds1,
                     color_cells_by = 'cluster',
                     label_groups_by_cluster = TRUE,
                     label_branch_points = TRUE,
                     label_roots = FALSE,
                     label_leaves = FALSE,
                     group_label_size = 5)+
  theme(legend.position = "right")

c1.1.plot|c1.2.plot
c1.1.plot|c.plot


# ...5. Order the cells in pseudotime -------------------

cds1 <- order_cells(cds1, reduction_method = 'UMAP')
# Change parameter of root_cells to which cells are at the beginning of the trajectory based on
# prior knowledge (which in our case should be the original cancer line i can't remember the name)
# give the parameters all the cellID of the original/root group of cells.
# if no root nodes are provided, a graphical user interface will be launched 
# This is for selecting one or more root nodes. We can pick as many as we want

plot_cells(cds1,
           color_cells_by = 'pseudotime',
           label_groups_by_cluster = FALSE,
           label_branch_points = TRUE,
           label_roots = TRUE,
           label_leaves = TRUE,
           group_label_size = 5,
           cell_size = 1.5)

# cells ordered by monocle3 pseudotime

pseudotime(cds1)
cds1$monocle3_pseudotime <- pseudotime(cds1)
data.pseudo.GSE202695 <- as.data.frame(colData(cds1))

ggplot(data.pseudo.GSE202695, aes(monocle3_pseudotime, reorder(model, monocle3_pseudotime, median), fill = type)) +
  geom_boxplot()

# Another method to subset cells based on the branching you see is
cds_sub <- choose_graph_segments(cds)
# However this would be done after trajectory generation
# And I am not sure if this can be done for the 3D trajectory


# Start of error - I could generate the helper function
# But applying the function for analysis generates an error instead
# Could be due to lack of data or simply due to my inexperience

# Programmatic specification of root node
# a helper function to identify the root principal points:
get_earliest_principal_node <- function(cds1, time_bin="130-170"){
  cell_ids <- which(colData(cds1)[, "embryo.time.bin"] == time_bin) # is embryo.time.bin universal?
  
  closest_vertex <-
    cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
  closest_vertex <- as.matrix(closest_vertex[colnames(cds1), ])
  root_pr_nodes <-
    igraph::V(principal_graph(cds1)[["UMAP"]])$name[as.numeric(names
                                                               (which.max(table(closest_vertex[cell_ids,]))))]
  
  root_pr_nodes
}
cds1 <- order_cells(cds1, root_pr_nodes = get_earliest_principal_node(cds1))
# utilizing this function would be useful for downstream 3D trajectory
# ONLY if we opt to use above helper function downstream
# End of error - Note: verify integrity with other datasets if possible


# ...6. Working with 3D trajectories---------------------------------------------
cds_3d <- reduce_dimension(cds1, reduction_method = 'UMAP', max_components = 3)
cds_3d <- cluster_cells(cds_3d)
cds_3d <- learn_graph(cds_3d)

# Check clustering - Useful to verify position of root node if unclear
cds_3d_plot_obj1 <- plot_cells_3d(cds_3d, color_cells_by="type")
cds_3d_plot_obj2 <- plot_cells_3d(cds_3d, color_cells_by="model")

cds_3d <- order_cells(cds_3d)
# root_pr_nodes= get_earliest_principal_node(cds1) was removed
# if function used upstream, use code below
cds_3d <- order_cells(cds_3d, root_pr_nodes=get_earliest_principal_node(cds1))

cds_3d_plot_obj3 <- plot_cells_3d(cds_3d, color_cells_by="pseudotime")

cds_3d_plot_obj1
cds_3d_plot_obj2
cds_3d_plot_obj3

# RStudio could visualize the 3D model for us
# If you wish to publish the results on an HTML page: Download the necessary packages.

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
FeaturePlot(MDAMB231.seu, features = c('TMSB10', 'TOP2A', 'CCNB1', 'TPX2', 'ANP32E', 'RPL41'))
# gene list incomplete, need to analyze data more
# Alternatively, come back to this after downstream analysis
FeaturePlot(MDAMB231.seu, features = c('TMSB10', 'TOP2A', 'CCNB1', 'TPX2', 'ANP32E', 'RPL41'))|cluster.names

# visualizing pseudotime in seurat as a feature plot

MDAMB231.seu$pseudotime <- pseudotime(cds1)
Idents(MDAMB231.seu) <- MDAMB231.seu$type
FeaturePlot(MDAMB231.seu, features = "pseudotime", label = T)

# Within each cluster/over all the groups identify most differentially expressed gene and lit review.


# Note: additional steps can be found at
# https://cole-trapnell-lab.github.io/monocle3/docs/differential/

