
source("scripts/utils/new_utils.R")
out0<-"outputs/figures_epi_response"
out<-fp(out0,"figure1")
dir.create(out,recursive = T)

#EPIGENETIC PROGRAMMING
#figure 1
#1A : volcanoplot
res<-fread("outputs/01-lga_vs_ctrl_limma_DMCs_analysis/res_limma.tsv.gz")
res[,meth.change:=logFC][,p_val_adj:=adj.P.Val][,p_val:=P.Value]

ggplot(res)+
  geom_point(aes(x=meth.change,y=-log10(p_val),col=p_val_adj<0.1&abs(meth.change)>25),size=0.1)+
  scale_color_manual(values = c("grey","red"))+theme_minimal()
ggsave(fp(out,"1A-volcano_plot_hypermet_LGA.png"))
ggsave(fp(out,"1A-volcano_plot_hypermet_LGA.pdf"))

#1B : pathway GSEA
library(enrichplot)

#kegg
res_kegg<-readRDS("outputs/03-pathway_analysis/res_gsea_kegg.rds")

pdf(fp(out,"1B-emapplot_gsea_kegg.pdf"),width = 14,height = 8)
emapplot(pairwise_termsim(res_kegg,showCategory = 30),showCategory = 30)
dev.off()

#go
res_go<-readRDS("outputs/03-pathway_analysis/res_gsea_go.rds")
pdf(fp(out,"1B-emapplot_gsea_go.pdf"),width = 14,height = 8)
emapplot(pairwise_termsim(res_go,showCategory = 40),showCategory = 40)
dev.off()

gsea_go<-fread("outputs/03-pathway_analysis/res_gsea_go.csv")

pdf(fp(out,"1B-dotplot_gsea_go_x_avg_genescore.pdf"),width = 8,height = 10)
dotplot(res_go,x=gsea_go[order(p.adjust)]$gene_score.avg[1:40],showCategory=40)
dev.off()

# both :
df_kegg_go<-merge(res_kegg,res_go,all=T) #merge in a dataframe
head(df_kegg_go)
#need transfoom in enrichment results :
# renv::install("bioc::ComplexHeatmap") #need this dependensies 
# renv::install("jmw86069/multienrichjam")
library(multienrichjam)
res_kegg_go<-enrichDF2enrichResult(df_kegg_go,keyColname = "ID",pvalueColname = "p.adjust",geneColname = "leading_edge")
data.table(as.data.frame(res_kegg_go))[!str_detect(ID,"GO")]

emapplot(pairwise_termsim(res_kegg_go,showCategory = 80),showCategory = 80)#not well informative


#gwas :
res_gwas<-readRDS("outputs/03-pathway_analysis/res_gsea_gwas.rds")
pdf(fp(out,"1B-emapplot_gsea_gwas.pdf"),width = 14,height = 8)
emapplot(pairwise_termsim(res_gwas,showCategory = 86),showCategory = 86)
dev.off()

pdf(fp(out,"1B-dotplot_gsea_gwas.pdf"),width = 8,height = 12)
dotplot(res_gwas,showCategory = 86)
dev.off()

#1C : TF motif enrichment
#with own background
res<-fread("outputs/03B-motif_analysis/knownResults.txt",
           select = c(1,2,3,5,6,7,8,9),
           col.names = c("motif","consensus","pval","padj","n_dmc_with_motif","pct_dmc_with_motif","n_background_with_motif","pct_background_with_motif"))
res[,pct_dmc_with_motif:=as.numeric(str_remove(pct_dmc_with_motif,"%"))]
res[,pct_background_with_motif:=as.numeric(str_remove(pct_background_with_motif,"%"))]
res[,motif:=str_remove(motif,"/Homer")]
res[,fold.enrichment:=pct_dmc_with_motif/pct_background_with_motif]
ggplot(res[padj<=0.2&  n_dmc_with_motif>30][order(pval)])+geom_point(aes(x=motif,col=-log10(pval),size=fold.enrichment,y=n_dmc_with_motif))+
  scale_x_discrete(limits=res[padj<=0.2&  n_dmc_with_motif>30][order(pval)]$motif)+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 70, hjust=1))

ggsave("outputs/figures_epi_response/figure1/1C-motif_enrichment_homer_own_background.pdf")

#with auto background
res<-fread("outputs/03B-motif_analysis/with_auto_background/knownResults.txt",
           select = c(1,2,3,5,6,7,8,9),
           col.names = c("motif","consensus","pval","padj","n_dmc_with_motif","pct_dmc_with_motif","n_background_with_motif","pct_background_with_motif"))
res[,pct_dmc_with_motif:=as.numeric(str_remove(pct_dmc_with_motif,"%"))]
res[,pct_background_with_motif:=as.numeric(str_remove(pct_background_with_motif,"%"))]
res[,motif:=str_remove(motif,"/Homer")]
res[,fold.enrichment:=pct_dmc_with_motif/pct_background_with_motif]
ggplot(res[padj<=0.2&  n_dmc_with_motif>30][order(pval)])+geom_point(aes(x=motif,col=-log10(pval),size=fold.enrichment,y=n_dmc_with_motif))+
  scale_x_discrete(limits=res[padj<=0.2&  n_dmc_with_motif>30][order(pval)]$motif)+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 70, hjust=1))

ggsave("outputs/figures_epi_response/figure1/1C-motif_enrichment_homer_auto_background.pdf",width = 16)



#1supp : Network [to optimize]

networks<-readRDS("outputs/04-comethylation_analysis/megena_spearman_pfn_mca_outputs.rds")
networks

networks.sum<-readRDS("outputs/04-comethylation_analysis/summary_megena_outputs.rds")
head(networks.sum$module.table,100)

library(ggplot2)
library(ggraph)

g <- readRDS("outputs/04-comethylation_analysis/graph.rds")
pnet.obj <- plot_module(output.summary = networks.sum,PFN = g,subset.module = "c1_29",
	layout = "kamada.kawai",label.hubs.only = FALSE,
	gene.set = NULL,color.code =  "grey",
	output.plot = FALSE,out.dir = "modulePlot",col.names = c("magenta","green","cyan"),label.scaleFactor = 20,
	hubLabel.col = "black",hubLabel.sizeProp = 1,show.topn.hubs = Inf,show.legend = TRUE)
print(pnet.obj[[1]]) #not informative so next



#INFLUENCE ON HSPC
#figure 2 : 
out<-fp(out0,"figure2")
dir.create(out)
source("scripts/utils/new_utils.R")
library(Seurat)
hmap<-readRDS("../singlecell/outputs/02-hematopo_datasets_integration/hematomap_ctrls_sans_stress/hematomap_ctrls_sans_stress.rds")

#2A : UMAP
#rm cluster 18
hmap<-subset(hmap,cell_type!="18")

#reorder celltype and lineage levels
Idents(hmap)<-"cell_type"
levels(hmap)<-c("LT-HSC",
                "HSC-1",
                "HSC-2",
                "HSC-3",
                "HSC-4",
                "MPP",
                "LMPP",
                "CLP",
                "proB",
                "B cell",
                "T cell",
               "MPP-Ery",
               "EMP",
               "EMP-cycle",
               "ErP-cycle",
               "Mk/Er",
               "GMP",
               "GMP-cycle",
               "DC")
hmap[["cell_type"]]<-Idents(hmap)

Idents(hmap)<-"lineage"

levels(hmap)<-c("LT-HSC",
                "HSC",
                "MPP/LMPP",
                "Lymphoid",
                "B cell",
                "T cell",
                "Erythro-Mas",
                "Mk/Er",
                "Myeloid",
                "DC")
hmap[["lineage"]]<-Idents(hmap)
DimPlot(hmap,group.by = c("cell_type","lineage"),label = T)
ggsave(fp(out,"2A-umap.pdf"))


#2supp : key genes/ct :
DefaultAssay(hmap)<-"SCT"
key_genes_ct<-c("ID1","ID2","DUSP2", #LT-HSC
           "EGR1","AVP", #HSC-1
         "CXCL8","ZFP36", #HSC-2
         "IRF1","STAT1", #HSC-3/MkP
          "KLF2","TSC22D3", #HSC-4/proT
           "MLLT3","CDK6", #MPP
           "SELL","CD99", #LMPP
           "LTB", #CLP
         "VPREB1","IGLL1", # proB
           "IGHM","CD37", #B cell
           "TNFAIP3","CD7", #T cell
           "GATA2", #MPP-Ery
           "GATA1", #EMP
           "BIRC5","MKI67","TOP2A", #ErP-cycle
           "HDC", #Mast
           "PLEK","HBD", #Mk/Er
           "MPO","AZU1", #GMP
         "CST3","CD83") #DC
DotPlot(hmap,features = key_genes_ct,group.by = "cell_type",)
ggsave(fp(out,"2supp-key_genes_by_celltype.pdf"),width = 23)


#2supp : key genes by lineage
key_genes_lin<-c("ID1","ID2","DUSP2", #LT-HSC
           "EGR1","AVP", #HSC
           "MLLT3","CDK6", #MPP
           "SELL","CD99", #LMPP
           "LTB", #CLP
         "VPREB1","IGLL1", # proB
           "IGHM","CD37", #B cell
           "TNFAIP3","CD7", #T cell
           "GATA2", #MPP-Ery
           "GATA1", #EMP
           "MKI67","TOP2A", #ErP-cycle
           "PLEK","HBD", #Mk/Er
           "MPO","CEBPA","CTSG","AZU1", #GMP
         "CST3","CD83") #DC
DotPlot(hmap,features = key_genes_lin,group.by = "lineage")
ggsave(fp(out,"2supp-key_genes_by_lin.pdf"),width = 12,height = 12)

ps<-FeaturePlot(hmap,features = c("ID1","EGR1","AVP","GATA1","HDC","MPO","CST3","LTB","VPREB1"),combine = F)
ps<-lapply(ps, function(x)x+NoAxes()+NoLegend())
wrap_plots(ps)
ggsave(fp(out,"2supp-umap_key_genes_by_lin.pdf"),width = 12,height = 12)


#2supp : CTRL LGA basale
#all 
res_ctrl_lga<-fread("../singlecell/outputs/08-DEGs_LGA_no_stress/pseudobulk_deseq2_all_cbps/res_de_analysis_all_genes.csv")

ggplot(res_ctrl_lga,aes(x=log2FoldChange,y=-log10(padj),col=padj<0.05&abs(log2FoldChange)>0.6))+
  geom_point()+ 
  geom_label_repel(aes(label = ifelse(padj<0.05&
                                        abs(log2FoldChange)>0.6
                                      ,gene,"")),
                   max.overlaps = 3000,
                   box.padding   = 0.35,
                   point.padding = 0.5,
                   segment.color = 'grey50')+
  scale_color_manual(values = c("grey","red")) +
 # scale_x_continuous(limits = c(-10,10))+
  theme_minimal() +
  theme(legend.position = "bottom")
ggsave(fp(out,"2supp-volcano_lga_vs_ctrl_basal.pdf"))

#by lineage
res_ctrl_lga_lin<-fread("../singlecell/outputs/04-HSC_stimulation_response_lga_vs_ctrl/pseudobulk_deseq2/comparison_basal/res_all_lineages_all_genes.csv")

ggplot(res_ctrl_lga_lin,aes(x=log2FoldChange,y=-log10(padj),col=padj<0.05&abs(log2FoldChange)>0.6))+
  geom_point()+ 
  geom_label_repel(aes(label = ifelse(padj<0.05&
                                        abs(log2FoldChange)>0.6
                                      ,gene,"")),
                   max.overlaps = 3000,
                   box.padding   = 0.35,
                   point.padding = 0.5,
                   segment.color = 'grey50')+
  scale_color_manual(values = c("grey","red")) +
  facet_wrap("lineage")+
 # scale_x_continuous(limits = c(-10,10))+
  theme_minimal() +
  theme(legend.position = "bottom")
ggsave(fp(out,"2supp-volcano_lga_vs_ctrl_basal_by_lineage.pdf"))

#2B : Ab activation signature
res_hto_dup<-fread("../singlecell/outputs/03-HTOs_stimulation/duplicates/res_pseudobulk_analysis_all_genes.csv")
res_hto_dup[padj<0.05&abs(log2FoldChange)>0.6] #1417
res_hto_dup[padj<0.05&abs(log2FoldChange)>0.6&gene%in%c("SOCS3","HES1","JUN","EGR1")]
genes_to_highlight<-c("SOCS3","HES1","JUN","FOSB","EGR1","ID1","ID2","DUSP2","SOCS1","JUNB","DUSP2","DUSP1","PLK2","PLK3")

#volcano
ggplot(res_hto_dup,aes(x=log2FoldChange,y=-log10(padj),col=padj<0.05&abs(log2FoldChange)>0.6))+
  geom_point()+
  geom_label_repel(aes(label=ifelse(padj<0.05&abs(log2FoldChange)>0.6&gene%in%genes_to_highlight,gene,"")),
                   max.overlaps = 3000)+
  scale_color_manual(values = c("grey","red")) +
  theme_minimal() +
  theme(legend.position = "bottom")
ggsave(fp(out,"2B-volcano_antibody_activation_signature.pdf"))

#2B-supp pathway activation of diff/prolif signaling ?
library(clusterProfiler)
library(enrichplot)

res_kegg<-readRDS("../singlecell/outputs/03-HTOs_stimulation/duplicates/res_kegg.rds")

pdf(fp(out,"2Bsupp-kegg_enrichment_antibody_activation_signature.pdf"),width = 12)
emapplot(pairwise_termsim(res_kegg,showCategory = 36),showCategory = 36)
dev.off()

res_go<-readRDS("../singlecell/outputs/03-HTOs_stimulation/duplicates/res_go.rds")

pdf(fp(out,"2Bsupp-go_enrichment_antibody_activation_signature.pdf"),width = 12)
emapplot(pairwise_termsim(res_go,showCategory = 41),showCategory = 41)
dev.off()


#2B-supp : heatmap HTO by lineage
res_hto_lin<-fread("../singlecell/outputs/03-HTOs_stimulation/duplicates/by_lineage/res_de_analysis_all_genes.csv")
res_hto_all<-rbind(res_hto_lin,fread("../singlecell/outputs/03-HTOs_stimulation/duplicates/res_pseudobulk_analysis_all_genes.csv")[,lineage:="all_cbps"])

res_hto_mat<-dcast(res_hto_all[gene%in%hto_signature&lineage%in%c("all_cbps","LT-HSC","HSC","MPP/LMPP","Lymphoid","Myeloid","Erythro-Mas")],gene~lineage,value.var ="log2FoldChange")
head(res_hto_mat)
res_hto_mat<-as.matrix(data.frame(res_hto_mat,row.names = "gene"))
head(res_hto_mat)
res_hto_mat[is.na(res_hto_mat)]<-0
res_hto_mat_scaled<-t(scale(t(res_hto_mat),center = F))

pdf(fp(out,"2Bsupp_heatmap_antibody_activation_signature_by_lineage.pdf"))
pheatmap(res_hto_mat_scaled[which(res_hto_mat_scaled[,"all_cbps"]<quantile(res_hto_mat_scaled[,"all_cbps"],0.90)),c(7,2:6,1)],cluster_cols = F,show_rownames = F,na_col = "grey")
dev.off()

#2Ca-heatmap signature
library(pheatmap)
hto_signature<-res_hto_dup[padj<0.05&abs(log2FoldChange)>0.6]$gene
res_hto<-Reduce(rbind,list(fread("../singlecell/outputs/03-HTOs_stimulation/pseudobulk_DESeq2_ctrl_hto/res_all_cbps_de_analysis.csv")[,compa:="ctrl_hto"],
                           fread("../singlecell/outputs/03-HTOs_stimulation/pseudobulk_DESeq2_lga_hto/res_all_cbps_de_analysis.csv")[,compa:="lga_hto"],
                           fread("../singlecell/outputs/03-HTOs_stimulation/pseudobulk_DESeq2_lga_vs_ctrl_hto/res_all_cbps_de_analysis.csv")[,compa:="lga_vs_ctrl_hto"]
))

res_hto_mat<-dcast(res_hto[gene%in%hto_signature],gene~compa,value.var ="log2FoldChange")
head(res_hto_mat)
res_hto_mat<-as.matrix(data.frame(res_hto_mat,row.names = "gene"))
head(res_hto_mat)

pdf(fp(out,"2Ca-heatmap_antibody_activation_signature_by_group.pdf"),width = 3)
pheatmap(res_hto_mat,cluster_cols = F,show_rownames = F)
dev.off()

#2Cb-LGA CTRL HTO
res_lga_ctrl_hto<-fread("../singlecell/outputs/03-HTOs_stimulation/pseudobulk_DESeq2_lga_vs_ctrl_hto/res_all_cbps_de_analysis.csv")
ggplot(res_lga_ctrl_hto,aes(x=log2FoldChange,y=-log10(padj),col=padj<0.05&abs(log2FoldChange)>0.6))+
  geom_point()+
  scale_color_manual(values = c("grey","red")) +
    scale_alpha_manual(values = c(0.5,1))+
  theme_minimal() +
  theme(legend.position = "bottom")
ggsave(fp(out,"2Cb-volcano_LGA_vs_Ctrl_activation_response.pdf"))




res_lga_ctrl_hto[padj<0.05&abs(log2FoldChange)>0.6,deg_sig:="deg"]

res_lga_ctrl_hto[padj<0.05&abs(log2FoldChange)>0.6&gene %in% hto_signature,deg_sig:="deg_hto"]
res_lga_ctrl_hto[is.na(deg_sig),deg_sig:="other_gene"]


ggplot(res_lga_ctrl_hto,aes(x=log2FoldChange,y=-log10(padj),col=deg_sig))+
  geom_point()+
   geom_label_repel(aes(label=ifelse(padj<0.05&abs(log2FoldChange)>0.6&gene%in%genes_to_highlight,gene,"")),
                   max.overlaps = 3000)+
  scale_color_manual(values = c("red","gold","grey")) +
  theme_minimal() +
  theme(legend.position = "bottom")

ggsave(fp(out,"2Cb-volcano_LGA_vs_Ctrl_activation_response_signature_highlight.pdf"))

#2D : HSC activation_response difference
res_lga_ctrl_hto_hsc<-fread("../singlecell/outputs/04-HSC_stimulation_response_lga_vs_ctrl/pseudobulk_deseq2/res_HSC_de_analysis.csv")

res_lga_ctrl_hto_hsc[padj<0.05&abs(log2FoldChange)>0.6,deg_sig:="deg"]
res_lga_ctrl_hto_hsc[padj<0.05&abs(log2FoldChange)>0.6&gene %in% hto_signature,deg_sig:="deg_hto"]
res_lga_ctrl_hto_hsc[is.na(deg_sig),deg_sig:="other_gene"]


ggplot(res_lga_ctrl_hto_hsc,aes(x=log2FoldChange,y=-log10(padj),col=deg_sig))+
  geom_point()+
   geom_label_repel(aes(label=ifelse(padj<0.05&abs(log2FoldChange)>0.6&gene%in%genes_to_highlight,gene,"")),
                   max.overlaps = 3000)+
  scale_color_manual(values = c("red","gold","grey")) +
  theme_minimal() +
  theme(legend.position = "bottom")

res_lga_ctrl_hto_hsc[deg_sig=="deg_hto"]$gene
res_lga_ctrl_hto_hsc[gene=="HES1"]

ggsave(fp(out,"2D-volcano_LGA_vs_Ctrl_HSC_activation_response_signature_highlight.pdf"))


#2Dsupp- for all lineage
res_lga_ctrl_hto_lin<-fread("../singlecell/outputs/04-HSC_stimulation_response_lga_vs_ctrl/pseudobulk_deseq2/res_all_lineages_all_genes.csv")

res_lga_ctrl_hto_lin[padj<0.05&abs(log2FoldChange)>0.6,deg_sig:="deg"]
res_lga_ctrl_hto_lin[padj<0.05&abs(log2FoldChange)>0.6&gene %in% hto_signature,deg_sig:="deg_hto"]
res_lga_ctrl_hto_lin[is.na(deg_sig),deg_sig:="other_gene"]

ggplot(res_lga_ctrl_hto_lin[!lineage%in%c("18","B cell","DC","Mk/Er","T cell")],aes(x=log2FoldChange,y=-log10(padj),col=deg_sig))+
  geom_point(size=0.3)+
  facet_wrap("lineage",scales = "free_y")+
  scale_color_manual(values = c("red","gold","grey")) +
  theme_minimal() +
  theme(legend.position = "bottom")

ggsave(fp(out,"2Dsupp-volcano_LGA_vs_Ctrl_All_Lineages_activation_response_signature_highlight.pdf"))

#2E : pathways HSC up and dn
library(clusterProfiler)
library(enrichplot)
kegg_dn<-readRDS("../singlecell/outputs/04-HSC_stimulation_response_lga_vs_ctrl/pseudobulk_deseq2/res_kegg_dn_padj0.2.rds")
kegg_dn_dt<-fread("../singlecell/outputs/04-HSC_stimulation_response_lga_vs_ctrl/pseudobulk_deseq2/res_kegg_dn_padj0.2.csv")
kegg_dn_dt[p.adjust<0.2]
emapplot(pairwise_termsim(kegg_dn,showCategory = 44),showCategory = 44)
ggsave(fp(out,"2E-emmaplot_kegg_dn_padj0.2.pdf"))


go_mf_up<-readRDS("../singlecell/outputs/04-HSC_stimulation_response_lga_vs_ctrl/pseudobulk_deseq2/res_go_mf_up.rds")
go_mf_up_dt<-fread("../singlecell/outputs/04-HSC_stimulation_response_lga_vs_ctrl/pseudobulk_deseq2/res_go_mf_up.csv.gz")
go_mf_up_dt[p.adjust<0.15]
emapplot(pairwise_termsim(go_mf_up,showCategory = 32),showCategory = 32)
ggsave(fp(out,"2E-emmaplot_go_mf_up.pdf"))



go_mf_dn<-readRDS("../singlecell/outputs/04-HSC_stimulation_response_lga_vs_ctrl/pseudobulk_deseq2/res_go_mf_dn.rds")
go_mf_dn_dt<-fread("../singlecell/outputs/04-HSC_stimulation_response_lga_vs_ctrl/pseudobulk_deseq2/res_go_mf_dn.csv.gz")
go_mf_dn_dt[p.adjust<0.15]
emapplot(pairwise_termsim(go_mf_dn,showCategory = 20),showCategory = 20)
ggsave(fp(out,"2E-emmaplot_go_mf_dn.pdf"))


go_bp_dn<-readRDS("../singlecell/outputs/04-HSC_stimulation_response_lga_vs_ctrl/pseudobulk_deseq2/res_go_bp_dn.rds")
go_bp_dn_dt<-fread("../singlecell/outputs/04-HSC_stimulation_response_lga_vs_ctrl/pseudobulk_deseq2/res_go_bp_dn.csv.gz")
go_bp_dn_dt[p.adjust<0.05]
emapplot(pairwise_termsim(go_bp_dn,showCategory = 36),showCategory = 36)
ggsave(fp(out,"2E-emmaplot_go_bp_dn.pdf"))

go_bp_up<-readRDS("../singlecell/outputs/04-HSC_stimulation_response_lga_vs_ctrl/pseudobulk_deseq2/res_go_bp_up.rds")
go_bp_up_dt<-fread("../singlecell/outputs/04-HSC_stimulation_response_lga_vs_ctrl/pseudobulk_deseq2/res_go_bp_up.csv.gz")
go_bp_up_dt[p.adjust<0.2]
emapplot(pairwise_termsim(go_bp_up,showCategory = 36),showCategory = 36)
ggsave(fp(out,"2E-emmaplot_go_bp_dn.pdf"))

#=>same differentiation signalling affected than methylation, 
#+ methyltranferaqe activity upreg in HSC LGA
#programmation of activation response ?


#figure3 :
out0<-"outputs/figures_epi_response/"
out<-fp(out0,"figure3")
dir.create(out)
#3A : Heatmap TF by celltype => STAT3, JUN, EGR1 HSC-1 specific
cbps<-readRDS("../singlecell/outputs/cbps0_8.rds")
cbps<-subset(cbps,lineage_hmap!="18"&group%in%c("ctrl","lga")&ambigous==F&orig.ident!="cd34_hto1_0C1I1L")



regul_ct <- sapply(split(colnames(cbps), cbps@meta.data$cell_type_hmap),
                                     function(cells) {
                                       if(length(cells)>1){
                                         return(rowMeans(as.matrix(cbps@assays$SCENIC@data)[,cells]))
                                       }else if (length(cells)==1){
                                         return(cbps@assays$SCENIC@data[,cells])
                                           }else return(NA)

                                       })
regul_ct_scaled<-t(scale(t(regul_ct), center = T, scale=T))

regul_ct_dt<-melt(regul_ct_scaled)
colnames(regul_ct_dt)<-c("regulon","cell_type","relative_activity")
regul_ct_dt<-data.table(regul_ct_dt)

regul_ct_dt[,top15:=relative_activity>=sort(relative_activity,decreasing = T)[15],by='cell_type']


pdf(fp(out,"3A-heatmap_top15_tf_by_cell_type.pdf"),width = 8,height = 20)
ComplexHeatmap::Heatmap(regul_ct_scaled[unique(regul_ct_dt[top15==T]$regulon),], name="Regulon activity")
dev.off()

#by lin
regul_lin <- sapply(split(colnames(cbps), cbps@meta.data$lineage_hmap),
                                     function(cells) {
                                       if(length(cells)>1){
                                         return(rowMeans(as.matrix(cbps@assays$SCENIC@data)[,cells]))
                                       }else if (length(cells)==1){
                                         return(cbps@assays$SCENIC@data[,cells])
                                           }else return(NA)

                                       })
regul_lin_scaled<-t(scale(t(regul_lin), center = T, scale=T))

regul_lin_dt<-melt(regul_lin_scaled)
colnames(regul_lin_dt)<-c("regulon","lineage","relative_activity")
regul_lin_dt<-data.table(regul_lin_dt)

regul_lin_dt[,top15:=relative_activity>=sort(relative_activity,decreasing = T)[15],by='lineage']

pdf(fp(out,"3A-heatmap_top15_tf_by_lineage.pdf"),width = 6,height = 18)
ComplexHeatmap::Heatmap(regul_lin_scaled[unique(regul_lin_dt[top15==T]$regulon),], name="Regulon activity")
dev.off()


DefaultAssay(cbps)<-"SCENIC"
FeaturePlot(cbps,c("STAT3","EGR1","JUN","FOSB"),min.cutoff = 0.2,max.cutoff = "q95")
ggsave(fp(out,"3A2-umap_key_tf_activity_in_hsc.pdf"))


#03Asupp- heatmap ct_sample cluster bien ct entre eux ? [to finish]
cbps$lineage_hmap_sample<-paste(cbps$lineage_hmap,cbps$sample,sep="_")
regul_lin_sample <- sapply(split(colnames(cbps), cbps@meta.data$lineage_hmap_sample),
                                     function(cells) {
                                       if(length(cells)>1){
                                         return(rowMeans(as.matrix(cbps@assays$SCENIC@data)[,cells]))
                                       }else if (length(cells)==1){
                                         return(cbps@assays$SCENIC@data[,cells])
                                           }else return(NA)

                                       })

regul_lin_sample_scaled <- t(scale(t(regul_lin_sample), center = T, scale=T))
ComplexHeatmap::Heatmap(regul_lin_scaled, name="Regulon activity",
                        top_annotation ="" #to finish
                        show_row_names = F,show_column_names = F)


#3B: TF activity bias of STAT3 and JUN ?

res_tf_diff<-fread("../singlecell/outputs/05-SCENIC/cbps0-8_clean/regulon_activity_lga_vs_ctrl_HTO_by_cell_type.csv.gz")
  
res_tf_diff[,auc_change:=avg_log2FC]
ggplot(res_tf_diff[p_val_adj<0.001&lineage=="HSC"&abs(avg_log2FC)>0.05])+geom_col(aes(x=regulon,y=auc_change,fill=-log10(p_val_adj)))
ggsave(fp(out,"3B1-barplot_top7_tf_change_hsc_lga_vs_ctrl.pdf"))


tf_int_dt_mtd<-fread("../singlecell/outputs/05-SCENIC/cbps0-8_clean/activated_tf_of_interest_by_cells.csv.gz")
tf_int_dt_mtsl<-unique(tf_int_dt_mtd[lineage_hmap=="HSC"],by=c("sample_hto","lineage_hmap","regulon"))

ggplot(tf_int_dt_mtsl[lineage_hmap=="HSC"])+
  geom_boxplot(aes(x=hto,y=pct.activ,fill=group))+facet_wrap("regulon")
ggsave(fp(out,"3B2-boxplot_pct_activ_tf_of_interest_by_sample_hsc.pdf"))


tf_int_dt_mtsl[,pvalue:=wilcox.test(pct.activ[group=="lga"],pct.activ[group=="ctrl"])$p.value,by=.(regulon,hto,lineage_hmap)]
tf_int_dt_mtsl[pvalue<0.1]
unique(tf_int_dt_mtsl[pvalue<0.2],by=c("regulon","lineage_hmap","hto"))[,.(regulon,hto,lineage_hmap,pvalue)]
#    regulon   hto lineage_hmap     pvalue
# 1:  ARID5A FALSE          HSC 0.18065268
# 2:    EGR1  TRUE          HSC 0.10789211
# 3:     JUN  TRUE          HSC 0.18115218
# 4:    FOSB  TRUE          HSC 0.14185814
# 5:  ARID5A  TRUE          HSC 0.05927406


#3C: regulons epigenetically affected 
#=> JUN /STAT3 /EGR1 pathway activation epigenetically altered,
#which genes exactly ?

#3D : genes of this pathways hypermet and downregulated in HSC 

#INFLUENCE ON DIFF/PROLIF
#figure4 
#4A: cell pop switch

#4B : Pseudotime

#4supp : velocity



#SEXUAL DIMORPHISM
#figure 5 
#5A : volcano
#5B : DEGs expression
#5C : alteration ++ regulons in LGA F


