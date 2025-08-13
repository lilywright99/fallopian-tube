
library(monocle3)
library(anndata)
library(igraph)

ann<- anndata::read_h5ad('Users/user/Documents/Fallopiantube/newanalysis/finalepi_weigert.h5ad')

expression_matrix<-ann$X
cell_metadata<-ann$obs
gene_metadata<-ann$var

head(expression_matrix, c(5,5))

#set expression matrix as a matrix, transpose so genes=row 
expression_matrix <- as.matrix(expression_matrix)
expression_matrix<-t(expression_matrix)

head(gene_metadata)
colnames(gene_metadata)[7] <-'gene_short_name'

##extracting UMAP coordinates from scanpy analysis
umapcoords <- as.data.frame(ann$obsm[["X_umap"]])
rownames(umapcoords) <- rownames(ann$obs)
colnames(umapcoords)[1:2] <- c("UMAP1", "UMAP2")

##generate cell_data_set object, used by monocle to hold single cell expression data 
cds<- new_cell_data_set(expression_matrix, cell_metadata, gene_metadata)

##add umap coordinates to cds 
reducedDims(cds)$UMAP <- umapcoords
plot_cells(cds, reduction_method='UMAP')

plot_cells(cds, color_cells_by = 'lily_cell_type')
plot_cells(cds, genes=c("OVGP1", "FOXJ1"))

##add clustering and partition data, use 1 partition as all epithelial cells
scanpy_clusters<-ann$obs['leiden']
scanpy_clusters<-as.factor(scanpy_clusters)

scanpy_clusters <- as.factor(ann$obs$leiden)
names(scanpy_clusters) <- rownames(ann$obs)
all(colnames(cds) %in% names(scanpy_clusters)) 

cds@clusters[["UMAP"]]$clusters<- cds@colData$lily_cell_type
cds@colData$clusters <- cds@colData$lily_cell_type
cds@clusters[["UMAP"]]$partitions <- factor(rep("1", ncol(cds)))
cds@colData$partitions<- cds@clusters[["UMAP"]]$partitions
names(cds@clusters$UMAP$partitions) <- colnames(cds)


cds@clusters$UMAP$clusters <- scanpy_clusters
cds@clusters$UMAP$partitions <- factor(rep("1", ncol(cds)), names=colnames(cds))

head(cds@clusters$UMAP$clusters)

##final bidirectional graph 

graph_controls <- list(minimal_branch_len = 4,
                       orthogonal_proj_tip = F,
                       nn.k = 40,
                       prune_graph = T,
                       L1.sigma = 0.05
                      )

cds <- learn_graph(cds, learn_graph_control = graph_controls, close_loop = T,)

plot_cells(cds, color_cells_by = "clusters", cell_size=1)

plot_cells(cds, label_principal_points = T)

get_earliest_principal_node <- function(cds, cell_type="Progenitor"){
  cell_ids <- which(colData(cds)[, "lily_cell_type"] == cell_type)
  
  closest_vertex <-
    cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
  closest_vertex <- as.matrix(closest_vertex[colnames(cds), ])
  root_pr_nodes <-
    igraph::V(principal_graph(cds)[["UMAP"]])$name[as.numeric(names
                                                              (which.max(table(closest_vertex[cell_ids,]))))]
  
  root_pr_nodes
}
cds <- order_cells(cds, root_pr_nodes=get_earliest_principal_node(cds))

cds <- order_cells(cds)
                  
plot_cells(cds,
           color_cells_by = "pseudotime",
           label_cell_groups=FALSE,
           label_leaves=FALSE,
           label_branch_points=TRUE,
           graph_label_size=1.5,
           cell_size = 1,
           trajectory_graph_segment_size = 1.5,)
          
secretory <- choose_cells(cds)

secretory_genes <- c("OVGP1", "PAX8", "MSLN")
secretorylineage <- secretory[rowData(secretory)$gene_short_name %in% secretory_genes]
secretorylineage <- order_cells(secretory)

plot_genes_in_pseudotime(secretorylineage,
                         color_cells_by="lily_cell_type",
                         min_expr=0.5)

cdssubset <- choose_cells(cds)

subset_pr_test_res <- graph_test(cdssubset, neighbor_graph="UMAP", cores=4)
pr_deg_ids <- row.names(subset(subset_pr_test_res, q_value < 0.05))

### meta analysis data

ann<- anndata::read_h5ad('Users/user/Documents/Fallopiantube/newanalysis/oviductepitheliumforconcat.h5ad')

expression_matrix<-ann$X
cell_metadata<-ann$obs
gene_metadata<-ann$var

head(expression_matrix, c(5,5))

#set expression matrix as a matrix, transpose so genes=row 
expression_matrix <- as.matrix(expression_matrix)
expression_matrix<-t(expression_matrix)

head(gene_metadata)
colnames(gene_metadata)[7] <-'gene_short_name'

##extracting UMAP coordinates from scanpy analysissum(is.na(cell_metadata))
umapcoords <- as.data.frame(ann$obsm[["X_umap"]])
rownames(umapcoords) <- rownames(ann$obs)
colnames(umapcoords)[1:2] <- c("UMAP1", "UMAP2")
sum(is.na(expression_matrix))
expression_matrix[is.na(expression_matrix)] <- 0


##generate cell_data_set object, used by monocle to hold single cell expression data 
cds<- new_cell_data_set(expression_matrix, cell_metadata, gene_metadata)

##add umap coordinates to cds 
reducedDims(cds)$UMAP <- umapcoords
cds@reducedDims$UMAP <- umapcoords

plot_cells(cds, reduction_method='UMAP')

plot_cells(cds, color_cells_by = 'epithelial_celltypes')
plot_cells(cds, genes=c("OVGP1", "FOXJ1"))

##add clustering and partition data, use 1 partition as all epithelial cells
scanpy_clusters<-ann$obs['epithelial_celltypes']
scanpy_clusters<-as.factor(scanpy_clusters)

cds@clusters[["UMAP"]]$clusters<- cds@colData$lily_cell_type
cds@colData$clusters <- cds@colData$lily_cell_type
cds@clusters[["UMAP"]]$partitions <- factor(rep("1", ncol(cds)))
cds@colData$partitions<- cds@clusters[["UMAP"]]$partitions
names(cds@clusters$UMAP$partitions) <- colnames(cds)


cds@clusters$UMAP$clusters <- scanpy_clusters
cds@clusters$UMAP$partitions <- factor(rep("1", ncol(cds)), names=colnames(cds))

cds@clusters$UMAP$partitions <- factor(rep("1", ncol(cds)))  # Create a factor
names(cds@clusters$UMAP$partitions) <- colnames(cds)  # Assign names separate

head(cds@clusters$UMAP$clusters)

##final bidirectional graph 

graph_controls <- list(minimal_branch_len = 4,
                       orthogonal_proj_tip = F,
                       nn.k = 40,
                       prune_graph = T,
                       L1.sigma = 0.05
)

cds <- learn_graph(cds, learn_graph_control = graph_controls, close_loop = T,)

plot_cells(cds, color_cells_by = "clusters")

plot_cells(cds, label_principal_points = T)

get_earliest_principal_node <- function(cds, cell_type="Progenitor"){
  cell_ids <- which(colData(cds)[, "lily_cell_type"] == cell_type)
  
  closest_vertex <-
    cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
  closest_vertex <- as.matrix(closest_vertex[colnames(cds), ])
  root_pr_nodes <-
    igraph::V(principal_graph(cds)[["UMAP"]])$name[as.numeric(names
                                                              (which.max(table(closest_vertex[cell_ids,]))))]
  
  root_pr_nodes
}
cds <- order_cells(cds, root_pr_nodes=get_earliest_principal_node(cds))

cds <- order_cells(cds)

plot_cells(cds,
           color_cells_by = "pseudotime",
           label_cell_groups=FALSE,
           label_leaves=FALSE,
           label_branch_points=TRUE,
           graph_label_size=1.5,
           cell_size = 1,
           trajectory_graph_segment_size = 1.5,)

