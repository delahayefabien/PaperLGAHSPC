### Project Setup ==================================================================================
library(here)
out<- here("outputs", "12-Pseudotime")
dir.create(out, recursive = TRUE, showWarnings = FALSE, mode = "0775")


### Load Packages ==================================================================================
# renv::install("bioc::batchelor")
# renv::install('cole-trapnell-lab/leidenbase')
# renv::install("cole-trapnell-lab/monocle3")
# renv::install("satijalab/seurat-wrappers")

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratWrappers)
  library(monocle3)
  library(Matrix)
  library(ggplot2)
  library(patchwork)
  set.seed(1234)

  
})


### Tables and Figures Theme =======================================================================
# theme_set(theme_light())


### Functions ======================================================================================
fp<-function(...)file.path(...)


# a helper function to identify the root principal points:
get_earliest_principal_node <- function(cbp.cds, time_bin="LT-HSC"){
  cell_ids <- which(colData(cbp.cds)[, "cell_type"] == time_bin & colData(cbp.cds)[, "group_hto"] == "ctrlFALSE")
  
  closest_vertex <-
  cbp.cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
  closest_vertex <- as.matrix(closest_vertex[colnames(cbp.cds), ])
  root_pr_nodes <-
  igraph::V(principal_graph(cbp.cds)[["UMAP"]])$name[as.numeric(names
  (which.max(table(closest_vertex[cell_ids,]))))]
  
  root_pr_nodes
}


### Analysis =======================================================================================

cbps<-readRDS("outputs/06-integr_singlecell_cbps/cbps_filtered.rds")

DefaultAssay(cbps) <- "RNA"


cbps.cds <- as.cell_data_set(cbps)
cbps.cds <- cluster_cells(cds = cbps.cds, reduction_method = "UMAP")
cbps.cds <- learn_graph(cbps.cds, use_partition = TRUE)

cbps.cds <- order_cells(cbps.cds, root_pr_nodes=get_earliest_principal_node(cbp.cds))
cbps <- AddMetaData(
  object = cbps,
  metadata = cbps.cds@principal_graph_aux@listData$UMAP$pseudotime,
  col.name = "Pseudotime_ALL_LTHSC_CTRLonly"
)

saveRDS(cbps,file=fp(out,"cbps_RNA_Pseudotime_ComputRoot.rds"))

### Complete =======================================================================================
message("Success!", appendLF = TRUE)
