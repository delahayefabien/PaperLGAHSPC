library(Seurat)
library(Matrix)
library(Matrix.utils)
library(DESeq2)
source("scripts/utils/new_utils.R")
out<-"outputs/09-LGA_vs_Ctrl_Activated"
dir.create(out)
cbps<-readRDS("outputs/06-integr_singlecell_cbps/cbps_filtered.rds")

cbps_h<-subset(cbps,hto==T)

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
genes_of_interest<-c("SOCS3","HES1","JUN","FOS","JUNB","ZFP36","EGR1",
                      "DUSP2","DUSP1","FOSB","SOCS1","KLF2","KLF4",
                       "PLK2","PLK3","ID1","MYC","","ID2","IDS","RGCC")

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
ggsave("outputs/figures_epi_response/figure2/2D-pseudo_bulk_deseq2_by_lineage_lga_vs_ctrl_basal.pdf")

#sc 
#run 09A


