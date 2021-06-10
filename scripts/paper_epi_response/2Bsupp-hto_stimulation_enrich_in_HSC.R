#HSC enrichment after hto stim due to what ?
#due to hto mark mieux HSC que les autres ? 
source("scripts/utils/new_utils.R")
library(Seurat)
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



mtd_cbp4<-data.table(readRDS('../singlecell/outputs/01-Analyses_Individuelles/CBP4_HiDepth/cbp4_all_cells.rds')@meta.data,keep.rownames = "bc")
mtd_cbp4[,sample:=new.ID]

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

Idents(cbp4_all)<-"predicted.cell_type"
FeaturePlot(cbp4_all,c("GATA2","GATA1"), #EMP markers
            split.by = "dying_cand",
            reduction = "ref.umap",label=T)

VlnPlot(cbp4_all,c("GATA2","GATA1","HBD"), #EMP markers
        group.by = "dying_cand",idents = "EMP") #pas d'expr de ces markers d'EMP, d'autres markeurs que MT ?
cbp4_all$ct_himt<-paste(cbp4_all$predicted.cell_type,cbp4_all$dying_cand,sep="_")
others_lin<-paste0(unique(cbp4_all$predicted.cell_type),"_FALSE")
others_lin<-others_lin[!str_detect(others_lin,"EMP|Er")]
EMP_himt_markers<-FindMarkers(cbp4_all,
                              group.by = "ct_himt",
                              ident.1 = "EMP_TRUE",
                              ident.2 = others_lin,
                              only.pos = T)
EMP_himt_markers #no erythroid markers except MT


VlnPlot(cbp4_all,c("percent.mt"), 
        group.by = "predicted.cell_type",
            split.by = "dying_cand")

VlnPlot(cbp4_all,c("percent.mt","nFeature_RNA","nCount_RNA"), 
        group.by = "dying_cand",idents = "EMP")

VlnPlot(cbp4_all,c("percent.mt","nFeature_RNA","nCount_RNA"), 
        group.by = "dying_cand",idents = "MPP")
VlnPlot(cbp4_all,c("nCount_RNA"), 
        group.by = "dying_cand",
        idents = "MPP",log = T)

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
emapplot(pairwise_termsim(res_or_kegg)) 
or_kegg<-data.table(as.data.frame(res_or_kegg))
or_kegg[,genes:=paste(tr(geneID,tradEntrezInSymbol = T),collapse = "|"),by="ID"]
or_kegg[Description=="Signaling pathways regulating pluripotency of stem cells"]


res_or_go<-enrichGO(genes2$ENTREZID,ont = "BP",
                        OrgDb = org.Hs.eg.db,
                      pvalueCutoff = 0.05)
as.data.frame(res_or_go)
emapplot(pairwise_termsim(res_or_go,showCategory = 50),showCategory = 50) 

or_go<-data.table(as.data.frame(res_or_go))
or_go[,genes:=paste(tr(geneID,tradEntrezInSymbol = T),collapse = "|"),by="ID"]
or_go[Description=="response to oxidative stress"]
#EMP hi mt have activation response markers but doesnt seems to be EMP specific
#seems to have cluster in EMP for bad reason : due to high MT merkers as in EMP
FeaturePlot(hmap,"percent.mt",reduction = "ref.umap",max.cutoff = "q95")
FeaturePlot(hmap,c("MT-CO2","MT-CO3"),reduction = "ref.umap",label = T)

#check if EMP hi mt have lo precition score comapred to MPP / other prediction
VlnPlot(cbp4_all,c("predicted.cell_type.score"), 
        group.by = "predicted.cell_type",
            split.by = "dying_cand")

#there is clearly a 2 souspop in dyning cand : les mortes des activés


#CCL : 
#- les cellules filtrées a haut poucentage mito sont enrichis en MPP et EMP, et en pathway de reponse au stimulis/ activation de la differentiation
#=> il faudra les inclures dans l'analyse en prenant soin de distinguer cellules morte / mourante et activé
#- vaut mieux basé son annotation sur predicted.cell_type que predicted.lineage car plus accurate

#test redefine percent mt
cbp4_all<-readRDS("../singlecell/outputs/01-Analyses_Individuelles/CBP4_HiDepth/cbp4_all_cells_recov_mt_activ.rds")
VlnPlot(cbp4_all,c("percent.mt"), 
        group.by = "predicted.lineage")
Idents(cbp4_all)<-"predicted.lineage"
wrap_plots(lapply(levels(cbp4_all),function(lin)FeatureScatter(cbp4_all,"percent.mt","nCount_RNA",cells =WhichCells(cbp4_all,idents = lin ))+ggtitle(lin)),)


cbp4_all$status_cells<-paste(ifelse(cbp4_all$dying_cand,"filtered","kept"),ifelse(cbp4_all$percent.mt<ifelse(str_detect(cbp4_all$predicted.cell_type,"MPP"),75,50),"activ","dead"),sep="_")
VlnPlot(cbp4_all,c("predicted.cell_type.score"), 
        group.by = "predicted.cell_type",
            split.by = "status_cells")
VlnPlot(cbp4_all,"CDK6",group.by = "predicted.cell_type",
            split.by = "status_cells") #OK


cbp4_all$to_recover<-cbp4_all$status_cells=="filtered_activ"
sum(cbp4_all$to_recover) #1792
cbp4_recov<-subset(cbp4_all,dying_cand==F|to_recover==T)


#reassign samples

cbp4_recov.htos<-as.matrix(Read10X("~/RUN/Run_539_10x_standard/Output/cellranger_count_cbp4_tri/single_cell_barcode_539_HTO_cbp4b/outs/filtered_feature_bc_matrix/")$`Antibody Capture`)
rownames(cbp4_recov.htos)<-c("ctrlM555",
                          "ctrlM518",
                          "ctrlM537",
                          "lgaF551",
                          "lgaF543")

cbp4_recov[["HTO"]] <- CreateAssayObject(counts = cbp4_recov.htos[,colnames(cbp4_recov)])


# Normalize HTO data, here we use centered log-ratio (CLR) transformation
cbp4_recov <- NormalizeData(cbp4_recov, assay = "HTO", normalization.method = "CLR")
cbp4_recov <- HTODemux(cbp4_recov, assay = "HTO",positive.quantile = 0.90)
table(cbp4_recov$HTO_classification.global)

 # Doublet Negative  Singlet 
 #    1555     2273     3008 
#instead of
# Doublet Negative  Singlet 
#     1173     1636     2235

 

#sex based recovery
source("../singlecell/scripts/utils/HTO_utils.R")
cbp4_recov<-checkHTOSex(cbp4_recov,gene_male="RPS4Y1",gene_female="XIST")
# calculating pct of real singlet male/female cells expressing sex marker ( save in 'misc' of 'HTO' assay):
# for male :  90 % > 73 %  express the male gene
# for female :  98 % > 96% express the female gene
#  
# Based on expression of this sex biomarkers :
#Based on expression of this sex biomarkers :
 # - 622 doublet ( 9 %), 
 #      Doublet Sex          
 #           FALSE TRUE
 #  Doublet   1329  226
 #  Negative  2115  158
 #  Singlet   2770  238
 # 
 # - 1628 cells with 'HTO_maxID' not good ( 24 % )
 # -in which  565 singlets badly assign ( 19 % of the 3008 singlets)
 #      Bad sex assign          
 #           FALSE TRUE
 #  Doublet   1067  488
 #  Negative  1698  575
 #  Singlet   2443  565
#   Singlet   1771  464
cbp4_recov<-sexBasedHTOAssign(cbp4_recov)
#  => save 'signal_stat' and 'background_stat'(min,median.. without outliers of positive HTO) for each sample save in 'misc' of HTO assay
# create metadata 'new.ID' based on 'hash.ID'
#  flag true doublet 
#  162 / 905 doublet male_female, are really doublets => can recover 743 cells
# recovering 743 doublets ...
#  534 recover from HTO_maxID
#  209 recover from HTO_secondID
# clean / reassign singlets ...
# for that, creating the metadata 'second[HTO]_is_sex_diff' 327 singlets badly assign
# on which,  217 are recoverable because second HTO is sex different
# try to recover the x cells, if this one have i) a sex diff as third hto signal or ii) a Hto signal marge important enough
# for i) creating the metadata 'HTO_thirdID' 
# for ii) creating the metadata 'HTO_2_3_marge' 
# need to have a signal differences between 2nd and 3rd HTO > quantile 0.95  to reassign the singlet 
#  197 / 217 singlet recoverable can be recover with this margin130 singlets non recoverable and non sex_doublet, identify as 'Bad_HTO_assign' 
#  finally, trying to recover the 1698 negative recoverable cells...
# 1) readjust threshold for negative based on new minimum HTO counts accept for singlets.
# new cutoff :           sample cutoff
# ctrlM555 ctrlM555    186
# lgaF551   lgaF551    101
# lgaF543   lgaF543     48
# ctrlM537 ctrlM537    101
# ctrlM518 ctrlM518    110
# 
# new first assignation of negative cells :
# Bad_HTO_assign       ctrlM518       ctrlM537       ctrlM555        Doublet        lgaF543        lgaF551       Negative 
#             66            126            114             41           1064             92            107            663 
# try to recover the new doublets with the same process than before..
# 1) create metadata 'new.HTO_classif' based on 'new.ID', with the new doublet pairs ('sample1_sample2') or 'Multiplet' if >2 HTO positives 
# new HTO classif for previously negative cells :
#    Bad_HTO_assign          ctrlM518          ctrlM537 ctrlM537_ctrlM518          ctrlM555 ctrlM555_ctrlM518 ctrlM555_ctrlM537 
#                66               126               114                99                41                31                37 
#  ctrlM555_lgaF543  ctrlM555_lgaF551           lgaF543  lgaF543_ctrlM518  lgaF543_ctrlM537           lgaF551  lgaF551_ctrlM518 
#                13                31                92                40                57               107                83 
#  lgaF551_ctrlM537   lgaF551_lgaF543         Multiplet          Negative 
#                82                66               367               246 
# 
#  332 are doublet m_f, in which 307 are concordant with seurat assigation of second and maxID. 
#  so 307  are recoverable !
# 2) recovery... 
# 24 / 307 doublet male_female, are really doublets  
#  => can recover 283 cells
# recovering 283 doublets ...
# 283 recover from HTO_maxID
# 
#  283 recover from HTO_maxID
#  0 recover from HTO_secondID
# sample distribution before sex based reassign : 
#  Doublet ctrlM555 Negative  lgaF551  lgaF543 ctrlM537 ctrlM518 
#     1555      847     2273      526      760      365      510 
# 
# sample distribution after sex based reassign : 
# Bad_HTO_assign       ctrlM518       ctrlM537       ctrlM555        Doublet        lgaF543        lgaF551       Negative 
#            196            779            562            683           1831           1228            894            663 
# 
# create metadata 'new.HTO_classif.global'
#  
# distribution of the old global classification : 
#  Doublet Negative  Singlet 
#     1555     2273     3008 
# 
# distribution of the new global classification : 
# Bad_HTO_assign        Doublet       Negative        Singlet 
#            196           1831            663           4146 
# 
# save new sex assignation in 'sex' meta.data: 

#  2995 => 4146 singlets
#=> +1151 singlet..
saveRDS(cbp4_recov,"../singlecell/outputs/01-Analyses_Individuelles/CBP4_HiDepth/cbp4_all_cells_recov_mt_activ.rds")

#diff of HTO by activ ?
VlnPlot(cbp4_recov,"nCount_HTO",group.by="status_cells",log=T) #nop

#by ct ?
VlnPlot(cbp4_recov,"nCount_HTO",group.by="predicted.cell_type",log=T,pt.size = 0) #nop

#lot of change with ancienne assign ?
cbp4_recov_s<-subset(cbp4_recov,new.HTO_classif.global=="Singlet")
sum(cbp4_recov_s$sample == cbp4_recov_s$new.ID,na.rm = T) #2816 / 2824
#NOP DONC VALIDEYY
cbp4_recov_s$sample<-cbp4_recov_s$new.ID
saveRDS(cbp4_recov_s,"../singlecell/outputs/01-Analyses_Individuelles/CBP4_HiDepth/cbp4_singlet_recov_mt_activ.rds")

