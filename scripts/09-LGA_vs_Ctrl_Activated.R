library(Seurat)
library(Matrix)
library(Matrix.utils)
library(DESeq2)
source("scripts/utils/new_utils.R")
out<-"outputs/09-LGA_vs_Ctrl_Activated"
dir.create(out)
cbps<-readRDS("outputs/06-integr_singlecell_cbps/cbps_filtered.rds")

cbps_h<-subset(cbps,hto==T)

#pseudobulk all cbps
#get mtd of interest
mtd<-data.table(cbps_h@meta.data,keep.rownames = "bc")
mts<-unique(mtd,by=c("sample"))
#get counts and filter genes lowly express
counts<-as.matrix(cbps_h@assays$RNA@counts)
dim(counts) 

counts <- counts[rowSums(counts > 0) >= 100|rowSums(counts > 0)>=ncol(counts)*0.1, ] 
message(nrow(counts)," genes kept after filtering") 

  # Aggregate across cluster-sample groups
sample_counts <- t(aggregate.Matrix(t(counts[,mtd$bc]), 
                     groupings = mtd$sample, fun = "sum"))
#DEseq2_analysis
dds <- DESeqDataSetFromMatrix(sample_counts, 
                               colData = data.frame(mts,row.names="sample")[colnames(sample_counts),], 
                               design = ~ group+orig.ident+sex)

dds <- DESeq(dds)

res <- results(dds,contrast = c("group","lga","ctrl"),alpha = 0.05)

res<-data.table(as.data.frame(res),keep.rownames="gene")

res[padj<0.05]$gene
fwrite(res,fp(out,"res_pseudobulkDESeq2_all_cbps.csv.gz"))


#Distribution change
cbps_h$lineage_hmap<-factor(cbps_h$lineage_hmap,levels = c("LT-HSC","HSC","MPP/LMPP","Lymphoid","B cell","T cell","Erythro-Mas","Mk/Er","Myeloid","DC"))


mtd<-data.table(cbps_h@meta.data,keep.rownames="bc")
mtd[,n.sample:=.N,by="sample_hto"]
mtd[,pct.lin:=.N/n.sample,by=c("sample_hto","lineage_hmap")]

ggplot(unique(mtd,by=c("sample","lineage_hmap")))+
  geom_boxplot(aes(x=lineage_hmap,y=pct.lin,fill=group))
ggsave("outputs/figures_epi_response/figure2/2d-distribution_lineage_control_lga_basal.pdf")

#expr change in lineage
#pseudobulk
Idents(cbps_h)<-"lineage_hmap"
res_lin<-Reduce(rbind,lapply(levels(cbps_h),function(lin){
  print(lin)
  cbps_sub<-subset(cbps_h,lineage_hmap==lin)
  #get mtd of interest
  mtd<-data.table(cbps_sub@meta.data,keep.rownames = "bc")
  mts<-unique(mtd,by=c("sample"))
  #get counts and filter genes lowly express
  counts<-as.matrix(cbps_sub@assays$RNA@counts)
  dim(counts) 

  counts <- counts[rowSums(counts > 0) >= 100|rowSums(counts > 0)>=ncol(counts)*0.1, ] 
  message(nrow(counts)," genes kept after filtering") 
  if(nrow(counts)>0){
    # Aggregate across cluster-sample groups
  sample_counts <- t(aggregate.Matrix(t(counts[,mtd$bc]), 
                       groupings = mtd$sample, fun = "sum"))
  #DEseq2_analysis
  dds <- DESeqDataSetFromMatrix(sample_counts, 
                                 colData = data.frame(mts,row.names="sample")[colnames(sample_counts),], 
                                 design = ~ group+orig.ident+sex)
  
  dds <- DESeq(dds)
  
  res <- results(dds,contrast = c("group","lga","ctrl"),alpha = 0.05)
  
  return(data.table(as.data.frame(res),keep.rownames="gene")[,lineage:=lin])
  }else{
    return(data.table())
      }
  

  }))

fwrite(res_lin,fp(out,"res_pseudobulkDESeq2_by_lineage.csv.gz"))

#volcano by lin
res_lin<-fread(fp(out,"res_pseudobulkDESeq2_by_lineage.csv.gz"))

genes_of_interest<-c("SOCS3","HES1","JUN","FOS","JUNB","ZFP36","EGR1",
                      "DUSP2","DUSP1","FOSB","SOCS1","KLF2","KLF4",
                       "PLK2","PLK3","ID1","MYC","","ID2","IDS","RGCC","PIK3R1","MT-ND3")

ggplot(res_lin[lineage%in%c("LT-HSC","HSC","MPP/LMPP","Erythro-Mas","Myeloid","Lymphoid")],aes(x=log2FoldChange,y=-log10(padj),col=padj<0.11&abs(log2FoldChange)>0.6))+
  geom_point()+ 
  geom_label_repel(aes(label = ifelse(padj<0.1&
                                        abs(log2FoldChange)>0.6&gene%in%genes_of_interest,gene,"")),
                   max.overlaps = 5000,
                   box.padding   = 0.35,
                   point.padding = 0.5,
                   segment.color = 'grey50')+
  facet_wrap("lineage")+
  scale_color_manual(values = c("grey","red")) +
  theme_minimal() +
  theme(legend.position = "bottom")
ggsave("outputs/figures_epi_response/figure2/2D-pseudo_bulk_deseq2_by_lineage_lga_vs_ctrl_activated.pdf")

#sc 
#run 09A
res_lin_act_sc<-fread("outputs/09-LGA_vs_Ctrl_Activated/res_scEdgeR_by_lineage.csv.gz")
table(res_lin_act_sc[p_val_adj<0.001&abs(avg_logFC)>0.6]$lineage_hmap)

res_lin_act_sc[p_val_adj<0.001&abs(avg_logFC)>0.6&lineage_hmap=="HSC"]

ggplot(res_lin_act_sc[lineage_hmap%in%c("LT-HSC","HSC","MPP/LMPP","Erythro-Mas","Myeloid","Lymphoid")],aes(x=avg_logFC,y=-log10(p_val_adj),col=p_val_adj<0.001&abs(avg_logFC)>0.6))+
  geom_point()+ 
  geom_label_repel(aes(label = ifelse(p_val_adj<0.001&
                                        abs(avg_logFC)>0.6&gene%in%genes_of_interest,gene,"")),
                   max.overlaps = 5000,
                   box.padding   = 0.35,
                   point.padding = 0.5,
                   segment.color = 'grey50')+
  facet_wrap("lineage_hmap")+
  scale_color_manual(values = c("grey","red")) +
  theme_minimal() +
  theme(legend.position = "bottom")

ggsave("outputs/figures_epi_response/figure2/2D-sc_edger_by_lineage_lga_vs_ctrl_activated.pdf")


#pathways/biological process up and dn by lineage [TO UPDATE] ####

library(clusterProfiler)
library(enrichplot)
library(org.Hs.eg.db)
res_dt<-fread("")
# -kegg
res_kegg<-enrichKEGG(bitr(res_dt[padj<0.05&abs(log2FoldChange)>0.6]$gene,fromType = "SYMBOL",
                             toType = "ENTREZID",OrgDb = org.Hs.eg.db)$ENTREZID,
                 organism = "hsa",pvalueCutoff = 1,qvalueCutoff = 1)
res_kegg_dt<-data.table(as.data.frame(res_kegg))
res_kegg_dt[p.adjust<0.5]#1
saveRDS(res_kegg,fp(out,"res_kegg.rds"))
fwrite(res_kegg_dt,fp(out,"res_kegg.csv"))

#lower thr
res_kegg2<-enrichKEGG(bitr(res_dt[padj<0.2&abs(log2FoldChange)>0.6]$gene,fromType = "SYMBOL",
                             toType = "ENTREZID",OrgDb = org.Hs.eg.db)$ENTREZID,
                 organism = "hsa",pvalueCutoff = 1,qvalueCutoff = 1)
res_kegg_dt<-data.table(as.data.frame(res_kegg2))
res_kegg_dt[p.adjust<0.2]#15
saveRDS(res_kegg2,fp(out,"res_kegg_padj0.2.rds"))
fwrite(res_kegg_dt,fp(out,"res_kegg_padj0.2.csv"))

#up
res_kegg_up<-enrichKEGG(bitr(res_dt[padj<0.2&log2FoldChange>0.6]$gene,fromType = "SYMBOL",
                             toType = "ENTREZID",OrgDb = org.Hs.eg.db)$ENTREZID,
                 organism = "hsa",pvalueCutoff = 1,qvalueCutoff = 1)
res_kegg_dt<-data.table(as.data.frame(res_kegg_up))
res_kegg_dt[p.adjust<0.2]#5
saveRDS(res_kegg_up,fp(out,"res_kegg_up_padj0.2.rds"))
fwrite(res_kegg_dt,fp(out,"res_kegg_up_padj0.2.csv"))

#dn
res_kegg_dn<-enrichKEGG(bitr(res_dt[padj<0.2&log2FoldChange<(-0.6)]$gene,fromType = "SYMBOL",
                             toType = "ENTREZID",OrgDb = org.Hs.eg.db)$ENTREZID,
                 organism = "hsa",pvalueCutoff = 1,qvalueCutoff = 1)
res_kegg_dt<-data.table(as.data.frame(res_kegg_dn))
res_kegg_dt[p.adjust<0.2]#44
saveRDS(res_kegg_dn,fp(out,"res_kegg_dn_padj0.2.rds"))
fwrite(res_kegg_dt,fp(out,"res_kegg_dn_padj0.2.csv"))


# -go

#molecular function
res_go_mf<-enrichGO(bitr(res_dt[padj<0.05&abs(log2FoldChange)>0.6]$gene,fromType = "SYMBOL",
                             toType = "ENTREZID",OrgDb = org.Hs.eg.db)$ENTREZID,
                 OrgDb = org.Hs.eg.db,pvalueCutoff = 1,qvalueCutoff = 1)

res_go_dt<-data.table(as.data.frame(res_go_mf))
res_go_dt[p.adjust<0.2]#47
saveRDS(res_go_mf,fp(out,"res_go_mf.rds"))
fwrite(res_go_dt,fp(out,"res_go_mf.csv.gz"))

#biological process
res_go_bp<-enrichGO(bitr(res_dt[padj<0.05&abs(log2FoldChange)>0.6]$gene,fromType = "SYMBOL",
                             toType = "ENTREZID",OrgDb = org.Hs.eg.db)$ENTREZID,ont = "BP",
                 OrgDb = org.Hs.eg.db,pvalueCutoff = 1,qvalueCutoff = 1)

res_go_dt<-data.table(as.data.frame(res_go_bp))
res_go_dt[p.adjust<0.1]#48
saveRDS(res_go_bp,fp(out,"res_go_bp.rds"))
fwrite(res_go_dt,fp(out,"res_go_bp.csv.gz"))

#up
res_go_mf_up<-enrichGO(bitr(res_dt[padj<0.05&log2FoldChange>0.6]$gene,fromType = "SYMBOL",
                             toType = "ENTREZID",OrgDb = org.Hs.eg.db)$ENTREZID,
                 OrgDb = org.Hs.eg.db,pvalueCutoff = 1,qvalueCutoff = 1)

res_go_dt<-data.table(as.data.frame(res_go_mf_up))
res_go_dt[p.adjust<0.2]#54, wiuth methyl transferase acctivity
saveRDS(res_go_mf_up,fp(out,"res_go_mf_up.rds"))
fwrite(res_go_dt,fp(out,"res_go_mf_up.csv.gz"))

res_go_bp_up<-enrichGO(bitr(res_dt[padj<0.05&log2FoldChange>0.6]$gene,fromType = "SYMBOL",
                             toType = "ENTREZID",OrgDb = org.Hs.eg.db)$ENTREZID,ont = "BP",
                 OrgDb = org.Hs.eg.db,pvalueCutoff = 1,qvalueCutoff = 1)

res_go_dt<-data.table(as.data.frame(res_go_bp_up))
res_go_dt[p.adjust<0.2]#0
saveRDS(res_go_bp_up,fp(out,"res_go_bp_up.rds"))
fwrite(res_go_dt,fp(out,"res_go_bp_up.csv.gz"))

#dn
res_go_mf_dn<-enrichGO(bitr(res_dt[padj<0.05&log2FoldChange<(-0.6)]$gene,fromType = "SYMBOL",
                             toType = "ENTREZID",OrgDb = org.Hs.eg.db)$ENTREZID,
                 OrgDb = org.Hs.eg.db,pvalueCutoff = 1,qvalueCutoff = 1)

res_go_dt<-data.table(as.data.frame(res_go_mf_dn))
res_go_dt[p.adjust<0.2]#21
saveRDS(res_go_mf_dn,fp(out,"res_go_mf_dn.rds"))
fwrite(res_go_dt,fp(out,"res_go_mf_dn.csv.gz"))

res_go_bp_dn<-enrichGO(bitr(res_dt[padj<0.05&log2FoldChange<(-0.6)]$gene,fromType = "SYMBOL",
                             toType = "ENTREZID",OrgDb = org.Hs.eg.db)$ENTREZID,ont = "BP",
                 OrgDb = org.Hs.eg.db,pvalueCutoff = 1,qvalueCutoff = 1)

res_go_dt<-data.table(as.data.frame(res_go_bp_dn))
res_go_dt[p.adjust<0.1]#108!
saveRDS(res_go_bp_dn,fp(out,"res_go_bp_dn.rds"))
fwrite(res_go_dt,fp(out,"res_go_bp_dn.csv.gz"))




