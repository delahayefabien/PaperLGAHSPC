library(Seurat)
library(Matrix)
library(Matrix.utils)
library(DESeq2)
source("scripts/utils/new_utils.R")
out<-"outputs/08-HTO_signature"
dir.create(out)

cbps<-readRDS("outputs/06-integr_singlecell_cbps/cbps_filtered.rds")

#0) signature : same sample hto vs not
#a) signature all cbps
library(Seurat)
library(Matrix)
library(Matrix.utils)
library(DESeq2)
source("../methyl/scripts/utils/new_utils.R")


cbps_dup<-subset(cbps,sample%in%c("ctrlM555","ctrlM518","ctrlM537"))
table(cbps_dup$hto,cbps_dup$sample)
#get mtd of interest
mtd<-data.table(cbps_dup@meta.data,keep.rownames = "bc")

mts<-unique(mtd,by=c("sample","orig.ident"))


#get counts and filter genes lowly express
counts<-as.matrix(cbps_dup@assays$RNA@counts)
dim(counts) 

counts <- counts[rowSums(counts > 0) >= 100, ] 
dim(counts) #13203 genes

# Aggregate across cluster-sample groups
sample_counts <- t(aggregate.Matrix(t(counts[,mtd$bc]), 
                       groupings = mtd$sample_hto, fun = "sum")) 


#DEseq2_analysis
dds <- DESeqDataSetFromMatrix(sample_counts, 
                               colData = data.frame(mts,row.names="sample_hto")[colnames(sample_counts),], 
                               design = ~ hto )

dds <- DESeq(dds)
resultsNames(dds)
res <- results(dds,name =  "htoTRUE",
                alpha = 0.05)


res_dt<-data.table(as.data.frame(res),keep.rownames="gene")

res_dt[is.na(padj),padj:=1]
res_dt[padj<0.05&abs(log2FoldChange)>0.6] 
fwrite(res_dt,fp(out,"res_pseudobulk_DESeq2_3replicates.csv"),sep=";")


#pathway enrichment of the signature
library(clusterProfiler)
library(enrichplot)
library(org.Hs.eg.db)

res_dt<-fread("outputs/08-HTO_signature/res_pseudobulk_DESeq2_3replicates.csv")

# -kegg
res_kegg<-enrichKEGG(bitr(res_dt[padj<0.05&abs(log2FoldChange)>0.6]$gene,fromType = "SYMBOL",
                             toType = "ENTREZID",OrgDb = org.Hs.eg.db)$ENTREZID,
                 organism = "hsa",pvalueCutoff = 0.05)
saveRDS(res_kegg,fp(out,"res_hto_signature_kegg.rds"))
res_kegg_dt<-data.table(as.data.frame(res_kegg))
res_kegg_dt#yes !
fwrite(res_kegg_dt,fp(out,"res_hto_signature_kegg.csv"))


# -go
res_go_mf<-enrichGO(bitr(res_dt[padj<0.05&abs(log2FoldChange)>0.6]$gene,fromType = "SYMBOL",
                             toType = "ENTREZID",OrgDb = org.Hs.eg.db)$ENTREZID,ont = "MF",
                 OrgDb = org.Hs.eg.db,pvalueCutoff = 0.05)

saveRDS(res_go_mf,fp(out,"res_hto_signature_go_mf.rds"))
res_go_mf_dt<-data.table(as.data.frame(res_go_mf))
fwrite(res_go_mf_dt,fp(out,"res_hto_signature_go_mf.csv"))

res_go_bp<-enrichGO(bitr(res_dt[padj<0.05&abs(log2FoldChange)>0.6]$gene,fromType = "SYMBOL",
                             toType = "ENTREZID",OrgDb = org.Hs.eg.db)$ENTREZID,ont = "BP",
                 OrgDb = org.Hs.eg.db,pvalueCutoff = 0.05)

saveRDS(res_go_bp,fp(out,"res_hto_signature_go_bp.rds"))
res_go_bp_dt<-data.table(as.data.frame(res_go_bp))
fwrite(res_go_bp_dt,fp(out,"res_hto_signature_go_bp.csv"))

#b) signature by lineage
out1<-fp(out,"by_lineage")
dir.create(out1)
Idents(cbps_dup)<-"lineage_hmap"
res_lin_dup<-Reduce(rbind,lapply(levels(cbps_dup),function(lin){
  print(lin)
  cbps_sub<-subset(cbps_dup,lineage_hmap==lin)
  #get mtd of interest
  mtd<-data.table(cbps_sub@meta.data,keep.rownames = "bc")
  mts<-unique(mtd,by=c("sample","orig.ident"))
  #get counts and filter genes lowly express
  counts<-as.matrix(cbps_sub@assays$RNA@counts)
  dim(counts) 

  counts <- counts[rowSums(counts > 0) >= 100|rowSums(counts > 0)>=ncol(counts)*0.1, ] 
  message(nrow(counts)," genes kept after filtering") 
  if(nrow(counts)>0){
    # Aggregate across cluster-sample groups
  sample_counts <- t(aggregate.Matrix(t(counts[,mtd$bc]), 
                       groupings = mtd$sample_hto, fun = "sum"))
  #DEseq2_analysis
  dds <- DESeqDataSetFromMatrix(sample_counts, 
                                 colData = data.frame(mts,row.names="sample_hto")[colnames(sample_counts),], 
                                 design = ~ hto)
  
  dds <- DESeq(dds)
  
  res <- results(dds,name = "htoTRUE",alpha = 0.05)
  
  return(data.table(as.data.frame(res),keep.rownames="gene")[,lineage:=lin])
  }else{
    return(data.table())
      }
  

  }))

fwrite(res_lin_dup,fp(out1,"res_pseudobulk_DESeq2_3replicates.csv.gz"))



#sign enrich by lineage
res_lin_dup<-fread("outputs/08-HTO_signature/by_lineage/res_pseudobulk_DESeq2_3replicates.csv.gz")

# -kegg
res_kegg_lin<-lapply(unique(res_lin_dup$lineage), function(lin)enrichKEGG(bitr(res_lin_dup[padj<0.05&abs(log2FoldChange)>0.6&lineage==lin]$gene,fromType = "SYMBOL",
                             toType = "ENTREZID",OrgDb = org.Hs.eg.db)$ENTREZID,
                 organism = "hsa",pvalueCutoff = 0.05))
names(res_kegg_lin)<-unique(res_lin_dup$lineage)
saveRDS(res_kegg_lin,fp(out1,"res_hto_signature_kegg_by_lineage.rds"))

res_kegg_lin_dt<-Reduce(rbind,lapply(names(res_kegg_lin)[!sapply(res_kegg_lin, is.null)],function(lin)data.table(as.data.frame(res_kegg_lin[[lin]]))[,lineage:=lin]))
res_kegg_lin_dt#yes !
fwrite(res_kegg_lin_dt,fp(out1,"res_hto_signature_kegg_by_lineage.csv"))


# -go
res_go_mf_lin<-lapply(unique(res_lin_dup$lineage), function(lin)enrichGO(bitr(res_lin_dup[padj<0.05&abs(log2FoldChange)>0.6&lineage==lin]$gene,fromType = "SYMBOL",
                             toType = "ENTREZID",OrgDb = org.Hs.eg.db)$ENTREZID,ont = "MF",
                 OrgDb = org.Hs.eg.db,pvalueCutoff = 0.05))
names(res_go_mf_lin)<-unique(res_lin_dup$lineage)

saveRDS(res_go_mf_lin,fp(out1,"res_hto_signature_go_mf_by_lineage.rds"))
res_go_mf_lin_dt<-Reduce(rbind,lapply(names(res_go_mf_lin)[!sapply(res_go_mf_lin, is.null)],function(lin)data.table(as.data.frame(res_go_mf_lin[[lin]]))[,lineage:=lin]))
res_go_mf_lin_dt
table(res_go_mf_lin_dt$lineage)
  # B cell          DC Erythro-Mas         HSC      LT-HSC    Lymphoid    MPP/LMPP     Myeloid 
  #        12          25           6          61          13          10          15          14 
fwrite(res_go_mf_lin_dt,fp(out1,"res_hto_signature_go_mf_by_lineage.csv"))


res_go_bp_lin<-lapply(unique(res_lin_dup$lineage), function(lin)enrichGO(bitr(res_lin_dup[padj<0.05&abs(log2FoldChange)>0.6&lineage==lin]$gene,fromType = "SYMBOL",
                             toType = "ENTREZID",OrgDb = org.Hs.eg.db)$ENTREZID,ont = "BP",
                 OrgDb = org.Hs.eg.db,pvalueCutoff = 0.05))
names(res_go_bp_lin)<-unique(res_lin_dup$lineage)

saveRDS(res_go_bp_lin,fp(out1,"res_hto_signature_go_bp_by_lineage.rds"))
res_go_bp_lin_dt<-Reduce(rbind,lapply(names(res_go_bp_lin)[!sapply(res_go_bp_lin, is.null)],function(lin)data.table(as.data.frame(res_go_bp_lin[[lin]]))[,lineage:=lin]))
res_go_bp_lin_dt
table(res_go_bp_lin_dt$lineage)
   # B cell          DC Erythro-Mas         HSC    Lymphoid    MPP/LMPP     Myeloid 
   #       57         210          98         443         129          46         155 
fwrite(res_go_bp_lin_dt,fp(out1,"res_hto_signature_go_bp_by_lineage.csv"))


#I) ctrl vs ctrl hto
source("scripts/utils/new_utils.R")
cbps<-readRDS("outputs/06-integr_singlecell_cbps/cbps_filtered.rds")
cbps_sub<-subset(cbps,group=="ctrl")


out<-"outputs/08-HTO_signature/pseudobulk_DESeq2_ctrl_hto"
dir.create(out,recursive=T)

#get mtd of interest
mtd<-data.table(cbps_sub@meta.data,keep.rownames = "bc")
table(mtd$hto)
# FALSE  TRUE 
# 18520  3985 

mts<-unique(mtd,by=c("sample","hto"))
table(mts$hto)
# FALSE  TRUE 
#     7     6 

#get counts and filter genes lowly express
counts<-as.matrix(cbps_sub@assays$RNA@counts)
dim(counts) #34889 22505

counts <- counts[rowSums(counts > 0) >= 100|rowSums(counts > 0)>=ncol(counts)*0.1, ] 
nrow(counts) # 14593 genes

  # Aggregate across cluster-sample groups
sample_counts <- t(aggregate.Matrix(t(counts[,mtd$bc]), 
                     groupings = mtd$sample_hto, fun = "sum"))

#DEseq2_analysis
dds <- DESeqDataSetFromMatrix(sample_counts, 
                               colData = data.frame(mts,row.names="sample_hto")[colnames(sample_counts),], 
                               design = ~ orig.ident+sex)

dds <- DESeq(dds)
resultsNames(dds)
# get the model matrix
mod_mat <- model.matrix(design(dds), colData(dds))
mod_mat
# calculate the vector of coefficient weights in hto
hto <- colMeans(mod_mat[dds$hto == T, ])
hto

# calculate the vector of coefficient weights in non hto
basal <- colMeans(mod_mat[dds$hto == F, ])
basal

# The contrast we are interested in is the difference between hto and not
hto - basal

# get the results for this contrast
res <- results(dds, contrast = hto - basal)
res_dt<-data.table(as.data.frame(res),keep.rownames="gene")
res_dt[padj<0.05&log2FoldChange>0.6]
res_dt[,lineage:="all_cbps"]

fwrite(res_dt,fp(out,"res_all_cbps_de_analysis.csv"),sep=";")


#II) same in lga
cbps_sub<-subset(cbps,group=="lga")

out<-"outputs/08-HTO_signature/pseudobulk_DESeq2_lga_hto"
dir.create(out,recursive=T)

#get mtd of interest

mtd<-data.table(cbps_sub@meta.data,keep.rownames = "bc")
table(mtd$hto)
# FALSE  TRUE 
# 16791  5097  
mts<-unique(mtd,by=c("sample","orig.ident"))
table(mts$hto)
# FALSE  TRUE 
#     6     4 
#get counts and filter genes lowly express
counts<-as.matrix(cbps_sub@assays$RNA@counts)
dim(counts) #34889 21888

counts <- counts[rowSums(counts > 0) >= 100|rowSums(counts > 0)>=ncol(counts)*0.1, ] 
nrow(counts) # 14659 genes

  # Aggregate across cluster-sample groups
sample_counts <- t(aggregate.Matrix(t(counts[,mtd$bc]), 
                     groupings = mtd$sample_hto, fun = "sum"))

#DEseq2_analysis
dds <- DESeqDataSetFromMatrix(sample_counts, 
                               colData = data.frame(mts,row.names="sample_hto")[colnames(sample_counts),], 
                               design = ~ orig.ident+sex)

dds <- DESeq(dds)
resultsNames(dds)
# get the model matrix
mod_mat <- model.matrix(design(dds), colData(dds))
mod_mat
# calculate the vector of coefficient weights in hto
hto <- colMeans(mod_mat[dds$hto == T, ])
hto

# calculate the vector of coefficient weights in non hto
basal <- colMeans(mod_mat[dds$hto == F, ])
basal

# The contrast we are interested in is the difference between hto and not
hto - basal

# get the results for this contrast
res <- results(dds, contrast = hto - basal)

res_dt<-data.table(as.data.frame(res),keep.rownames="gene")
res_dt[padj<0.05] #1369
res_dt[,lineage:="all_cbps"]
fwrite(res_dt,fp(out,"res_all_cbps_de_analysis.csv"),sep=";")

#IV) ctrl lga hto
cbps_sub<-subset(cbps,hto==T&group%in%c("ctrl","lga")&ambigous==F&orig.ident!="cd34_hto1_0C1I1L")

out<-"outputs/08-HTO_signature/pseudobulk_DESeq2_lga_vs_ctrl_hto"
dir.create(out,recursive=T)

#get mtd of interest

mtd<-data.table(cbps_sub@meta.data,keep.rownames = "bc")
table(mtd$hto)
#   TRUE 
#  9082
mts<-unique(mtd,by=c("sample","orig.ident"))
table(mts$hto)
# TRUE 
#  10 
#get counts and filter genes lowly express
counts<-as.matrix(cbps_sub@assays$RNA@counts)

counts <- counts[rowSums(counts > 0) >= 100|rowSums(counts > 0)>=ncol(counts)*0.1, ] 
nrow(counts) # 13613 genes

  # Aggregate across cluster-sample groups
sample_counts <- t(aggregate.Matrix(t(counts[,mtd$bc]), 
                     groupings = mtd$sample_hto, fun = "sum"))

#DEseq2_analysis
dds <- DESeqDataSetFromMatrix(sample_counts, 
                               colData = data.frame(mts,row.names="sample_hto")[colnames(sample_counts),], 
                               design = ~ group+orig.ident+sex)

dds <- DESeq(dds)
resultsNames(dds)

# get the model matrix
mod_mat <- model.matrix(design(dds), colData(dds))
mod_mat
# calculate the vector of coefficient weights in hto
lga <- colMeans(mod_mat[dds$group == "lga", ])
lga

# calculate the vector of coefficient weights in non hto
ctrl <- colMeans(mod_mat[dds$group == "ctrl", ])
ctrl

# The contrast we are interested in is the difference between hto and not
lga - ctrl

# get the results for this contrast
res <- results(dds, contrast = lga - ctrl)

res_dt<-data.table(as.data.frame(res),keep.rownames="gene")
res_dt[padj<0.05] #304
res_dt[,lineage:="all_cbps"]
fwrite(res_dt,fp(out,"res_all_cbps_de_analysis.csv"),sep=";")



#now that we have validated hsc response to stimulation, we can evaluate the differential response of lga HSC compared to control

#see 09-
