### Project Setup ==================================================================================
library(here)
out<- here("outputs", "12A-Pseudotime_integrated")
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
  library(data.table)
  set.seed(1234)

  
})


### Tables and Figures Theme =======================================================================
# theme_set(theme_light())


### Functions ======================================================================================
fp<-function(...)file.path(...)


# a helper function to identify the root principal points:
get_earliest_principal_node <- function(cbp.cds, time_bin="LT-HSC"){
  cell_ids <- which(colData(cbp.cds)[, "cell_type_hmap"] == time_bin & colData(cbp.cds)[, "group_hto"] == "ctrlFALSE")
  
  closest_vertex <-
  cbp.cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
  closest_vertex <- as.matrix(closest_vertex[colnames(cbp.cds), ])
  root_pr_nodes <-
  igraph::V(principal_graph(cbp.cds)[["UMAP"]])$name[as.numeric(names
  (which.max(table(closest_vertex[cell_ids,]))))]
  
  root_pr_nodes
}


### Analysis =======================================================================================

cbps<-readRDS("outputs/10A-classical_integr/cbps0-8_clean.rds")

mtd<-data.table(cbps@meta.data,keep.rownames = "bc")
mtd[,n.sample:=.N,"sample_hto"]
mtd[,pct.lin:=.N/n.sample,c("sample_hto","lineage_hmap")]

table(mtd$group)
table(mtd$hto)
table(mtd$sample_hto)



DefaultAssay(cbps) <- "RNA"

cbps[["pca"]]<-cbps@reductions$integrated.pca
cbps[["umap"]]<-cbps@reductions$integrated.umap

cbps.cds <- as.cell_data_set(cbps)
cbps.cds <- cluster_cells(cds = cbps.cds, reduction_method = "UMAP")
cbps.cds <- learn_graph(cbps.cds, use_partition = TRUE)

cbps.cds <- order_cells(cbps.cds, root_pr_nodes=get_earliest_principal_node(cbps.cds))

cbps <- AddMetaData(
  object = cbps,
  metadata = cbps.cds@principal_graph_aux@listData$UMAP$pseudotime,
  col.name = "Pseudotime_ALL_LTHSC_CTRLonly"
)

plot_cells(
  cds = cbps.cds,
  color_cells_by = "pseudotime",
  show_trajectory_graph = TRUE,
  label_branch_points= TRUE,
  label_leaves=FALSE
)

cbps_inf<-subset(cbps,Pseudotime_ALL_LTHSC_CTRLonly!=Inf)
Idents(cbps_inf)<-"cell_type_hmap"
FeaturePlot(cbps_inf,"Pseudotime_ALL_LTHSC_CTRLonly",reduction = 'umap',label=T)


p2<-FeaturePlot(cbps_inf,c("Pseudotime_ALL_LTHSC_CTRLonly"),reduction = 'umap',label=T)

p1<-DimPlot(cbps_inf,group.by = "seurat_clusters",label = T)
p1+p2

p3<-DimPlot(cbps_inf,group.by = "lineage_hmap",label = T)
p1+p3

p2|p3
m11<-FindMarkers

mtd<-data.table(cbps@meta.data,keep.rownames = "bc")

fwrite(mtd,fp(out,"metadata_RNA_Pseudotime_ComputRoot.csv"))
mtd<-fread("outputs/12A-Pseudotime_integrated")


ggplot(mtd)+geom_boxplot(aes(y=Pseudotime_ALL_LTHSC_CTRLonly,x=lineage_hmap),alpha=0.7)+theme_minimal()


unique(mtd$batch)
ggplot(mtd)+geom_density(aes(x=Pseudotime_ALL_LTHSC_CTRLonly,fill=group,col=group),alpha=0.7)+
  facet_grid(hto~)+theme_minimal()

mtd[,n.sample:=.N,"sample_hto"]

mtd[,pct.lin:=.N/n.sample,c("sample_hto","lineage_hmap")]

ggplot(unique(mtd,by=c("sample_hto","lineage_hmap")))+geom_boxplot(aes(x=hto,y=pct.lin,fill=group))+facet_wrap("lineage_hmap")+theme_minimal()

table(unique(mtd,by=c("sample_hto"))$group_hto)
# ctrlFALSE  ctrlTRUE  lgaFALSE   lgaTRUE 
#         7         6         6         4 



ggplot(mtd)+geom_bar(aes(x=group,fill=lineage_hmap),position = "fill")+facet_wrap("hto")+theme_minimal()
ggplot(unique(mtd,by=c("sample_hto","lineage_hmap"))[differentiated==F])+geom_boxplot(aes(x=hto,y=pct.lin,fill=group))+facet_wrap("lineage_hmap")+theme_minimal()


saveRDS(cbps,file=fp(out,"cbps_RNA_Pseudotime_ComputRoot.rds"))

### Complete =======================================================================================
message("Success!", appendLF = TRUE)
