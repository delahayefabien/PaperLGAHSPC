#HSC enrichment after hto stim due to what ?
#due to hto mark mieux HSC que les autres ? 

renv::install("bioc::glmGamPoi")
hmap<-readRDS("../singlecell/outputs/02-hematopo_datasets_integration/hematomap_ctrls_sans_stress/hematomap_ctrls_sans_stress.rds")
cbp4<-readRDS("../singlecell/outputs/01-Analyses_Individuelles/CBP4_HiDepth/cbp4_all_cells.rds")

cbp4<-SCTransform(cbp4,method = "glmGamPoi")
cbp4 <- CellCycleScoring(cbp4,s.features = cc.genes$s.genes,
                          g2m.features = cc.genes$g2m.genes,
                          set.ident = TRUE,
                          search=TRUE)


  
cbp4$CC.Difference <- cbp4$S.Score - cbp4$G2M.Score
 

cbp4<-SCTransform(cbp4,vars.to.regress=c("percent.mt","CC.Difference"),
                  return.only.var.genes=F, 
                  method = "glmGamPoi")

DefaultAssay(hmap)<-"integrated"
hmap[["pca.annoy.neighbors"]] <- LoadAnnoyIndex(object = hmap[["pca.annoy.neighbors"]], file = "../singlecell/outputs/02-hematopo_datasets_integration/cbps0_8_Map_on_hematomap/reftmp.idx")

anchors <- FindTransferAnchors(
    reference = hmap,
    query = cbp4,
    k.filter = NA,
    reference.reduction = "pca", 
    reference.neighbors = "pca.annoy.neighbors", 
    dims = 1:50
  )

cbp4 <- MapQuery(
    anchorset = anchors, 
    query = cbp4,
    reference = hmap, 
    refdata = list(
      cell_type = "cell_type", 
      lineage = "lineage"),
    reference.reduction = "pca",
    reduction.model = "ref.umap"
  )

cbp4$sample<-cbp4$new.ID

mtd_cbp4<-data.table(cbp4@meta.data,keep.rownames = "bc")
mtd_cbp4$predicted.lineage<-factor(mtd_cbp4$predicted.lineage,levels = c("LT-HSC","HSC","MPP/LMPP","Lymphoid","B cell","T cell","Erythro-Mas","Mk/Er","Myeloid","DC"))

ggplot(mtd_cbp4)+geom_bar(aes(x=sample,fill=predicted.lineage),position = "fill")

#primed / MPP CD34+ cells more in dying cells than HSC after hto stress ?

samples.umis<- Read10X("~/RUN/Run_539_10x_standard/Output/cellranger_count_cbp4_tri/single_cell_barcode_539_HTO_cbp4b/outs/filtered_feature_bc_matrix/")$`Gene Expression`

cbp4_all <- CreateSeuratObject(counts = samples.umis,project = "cbp4")
cbp4_all[["percent.mt"]] <- PercentageFeatureSet(object = cbp4_all, pattern = "^MT-")

VlnPlot(object = cbp4_all, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

cbp4_all<-SCTransform(cbp4_all,method = "glmGamPoi")
cbp4_all <- CellCycleScoring(cbp4_all,s.features = cc.genes$s.genes,
                          g2m.features = cc.genes$g2m.genes,
                          set.ident = TRUE,
                          search=TRUE)


  
cbp4_all$CC.Difference <- cbp4_all$S.Score - cbp4_all$G2M.Score
 

cbp4_all<-SCTransform(cbp4_all,vars.to.regress=c("percent.mt","CC.Difference"),
                  return.only.var.genes=F, 
                  method = "glmGamPoi")

anchors <- FindTransferAnchors(
    reference = hmap,
    query = cbp4_all,
    k.filter = NA,
    reference.reduction = "pca", 
    reference.neighbors = "pca.annoy.neighbors", 
    dims = 1:50
  )

cbp4_all <- MapQuery(
    anchorset = anchors, 
    query = cbp4_all,
    reference = hmap, 
    refdata = list(
      cell_type = "cell_type", 
      lineage = "lineage"),
    reference.reduction = "pca",
    reduction.model = "ref.umap"
  )


mtd_cbp4_all<-data.table(cbp4_all@meta.data,keep.rownames = "bc")
mtd_cbp4_all_s<-merge(mtd_cbp4_all,mtd_cbp4[,.(bc,sample)],all.x=T)
mtd_cbp4_all_s$predicted.lineage<-factor(mtd_cbp4_all_s$predicted.lineage,levels = c("LT-HSC","HSC","MPP/LMPP","Lymphoid","B cell","T cell","Erythro-Mas","Mk/Er","Myeloid","DC"))

ggplot(mtd_cbp4_all_s)+
  geom_bar(aes(x=sample,fill=predicted.lineage),
           position = "fill") 
ggplot(mtd_cbp4_all_s)+geom_bar(aes(x=sample,fill=predicted.cell_type),position = "fill")
#dying cells ++ MPP and EMP ! 
DimPlot(cbp4_all,group.by = "predicted.cell_type",label=T,reduction = "ref.umap")
DimPlot(cbp4_all,group.by = "predicted.lineage",label=T,reduction = "ref.umap")

DimPlot(cbp4_all,group.by = "predicted.cell_type",label=T,reduction = "ref.umap",split.by='dying_cand')
DimPlot(cbp4_all,group.by = "predicted.lineage",label=T,reduction = "ref.umap",split.by='dying_cand')

FeaturePlot(cbp4_all, "nFeature_RNA",reduction = "ref.umap")

#parce que pas bcp de genes/marqueurs ?

cbp4_all$sample<-data.frame(mtd_cbp4_all_s,row.names = "bc")[colnames(cbp4_all),"sample"]
DimPlot(cbp4_all,group.by = "sample",reduction = "ref.umap")

cbp4_all$dying_cand<-is.na(cbp4_all$sample)

table(cbp4_all$dying_cand) #2.5k/7.5 33% of dying cells
#have genes markers of MPP / EMP
FeaturePlot(cbp4_all,c("CDK6","AVP","MLLT3"), #MPP markers
            split.by = "dying_cand",
            reduction = "ref.umap")
#"dying" cells have expression of MPP markers

FeaturePlot(cbp4_all,c("GATA2","GATA1"), #EMP markers
            split.by = "dying_cand",
            reduction = "ref.umap")

VlnPlot(cbp4_all,c("GATA2","GATA1"), #EMP markers
        group.by = "predicted.cell_type",
            split.by = "dying_cand")

VlnPlot(cbp4_all,c("percent.mt"), 
        group.by = "predicted.cell_type",
            split.by = "dying_cand")

VlnPlot(cbp4_all,c("percent.mt","nFeature_RNA","nCount_RNA"), 
        group.by = "orig.ident",
            split.by = "dying_cand")

#so est ce que les cellules filtrés sont reellement mourrantes ?
#> enrichment for apoptosis/cell death pathway ?

#MPP
Idents(cbp4_all)<-"predicted.cell_type"
MPP_dying_cand_markers<-FindMarkers(cbp4_all,subset.ident = "MPP",only.pos = T,group.by = "dying_cand",ident.1 = "TRUE")
head(MPP_dying_cand_markers,100)
MPP_dying_cand_markers<-data.table(MPP_dying_cand_markers,keep.rownames = "gene")
MPP_dying_cand_markers[p_val_adj<0.01]$gene

#
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
res_or_kegg<-enrichKEGG(bitr(MPP_dying_cand_markers[p_val_adj<0.01]$gene,fromType = "SYMBOL",toType = "ENTREZID",OrgDb = org.Hs.eg.db)$ENTREZID,
                 organism = "hsa",pvalueCutoff = 0.05)

#24% failed to map
genes<-bitr(MPP_dying_cand_markers[p_val_adj<0.01]$gene,fromType = "SYMBOL",toType = "ENTREZID",OrgDb = org.Hs.eg.db)
genes #no MTgene

genes<-TransSymboltoEnsembl(MPP_dying_cand_markers[p_val_adj<0.01]$gene)
genes2<-bitr(genes$ensembl_gene_id,fromType = "ENSEMBL",toType = "ENTREZID",OrgDb = org.Hs.eg.db)
res_or_kegg<-enrichKEGG(genes2$ENTREZID,
                 organism = "hsa",pvalueCutoff = 0.05)
as.data.frame(res_or_kegg)
emapplot(pairwise_termsim(res_or_kegg)) #there is apoptosis pathway
or_kegg<-data.table(as.data.frame(res_or_kegg))
or_kegg[,genes:=paste(tr(geneID,tradEntrezInSymbol = T),collapse = "|"),by="ID"]
or_kegg[Description=="Apoptosis"]

res_or_kegg<-enrichKEGG(bitr(MPP_dying_cand_markers[p_val_adj<0.01]$gene,fromType = "SYMBOL",toType = "ENTREZID",OrgDb = org.Hs.eg.db)$ENTREZID,
                 organism = "hsa",pvalueCutoff = 0.05)



res_or_go<-enrichGO(genes2$ENTREZID,ont = "BP",
                        OrgDb = org.Hs.eg.db,
                      pvalueCutoff = 0.05)
as.data.frame(res_or_go)
emapplot(pairwise_termsim(res_or_go,showCategory = 50),showCategory = 50) #there is apoptosis pathway

or_go<-data.table(as.data.frame(res_or_go))
or_go[,genes:=paste(tr(geneID,tradEntrezInSymbol = T),collapse = "|"),by="ID"]
or_go[Description=="response to oxidative stress"]

#EMP
EMP_dying_cand_markers<-FindMarkers(cbp4_all,subset.ident = "EMP",only.pos = T,group.by = "dying_cand",ident.1 = "TRUE")
head(EMP_dying_cand_markers,100)
EMP_dying_cand_markers<-data.table(EMP_dying_cand_markers,keep.rownames = "gene")
EMP_dying_cand_markers[p_val_adj<0.01]$gene

genes<-TransSymboltoEnsembl(EMP_dying_cand_markers[p_val_adj<0.01]$gene)
genes2<-bitr(genes$ensembl_gene_id,fromType = "ENSEMBL",toType = "ENTREZID",OrgDb = org.Hs.eg.db)
res_or_kegg<-enrichKEGG(genes2$ENTREZID,
                 organism = "hsa",pvalueCutoff = 0.05)
as.data.frame(res_or_kegg)
emapplot(pairwise_termsim(res_or_kegg)) #there is apoptosis pathway
or_kegg<-data.table(as.data.frame(res_or_kegg))
or_kegg[,genes:=paste(tr(geneID,tradEntrezInSymbol = T),collapse = "|"),by="ID"]
or_kegg[Description=="Signaling pathways regulating pluripotency of stem cells"]


res_or_go<-enrichGO(genes2$ENTREZID,ont = "BP",
                        OrgDb = org.Hs.eg.db,
                      pvalueCutoff = 0.05)
as.data.frame(res_or_go)
emapplot(pairwise_termsim(res_or_go,showCategory = 50),showCategory = 50) #there is apoptosis pathway

or_go<-data.table(as.data.frame(res_or_go))
or_go[,genes:=paste(tr(geneID,tradEntrezInSymbol = T),collapse = "|"),by="ID"]
or_go[Description=="response to oxidative stress"]

#CCL : 
#- les cellules filtrées a haut poucentage mito sont enrichis en MPP et EMP, et en pathway de reponse au stimulis/ activation de la differentiation
#=> il faudra les inclures dans l'analyse en prenant soin de distinguer cellules morte / mourante et activé
#- vaut mieux basé son annotation sur predicted.cell_type que predicted.lineage car plus accurate
