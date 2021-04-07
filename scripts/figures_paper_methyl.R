#figures paper methyl
source("../methyl/scripts/utils/new_utils.R")
library(limma)
library(pheatmap)
out<-"../methyl/outputs/figures_paper_meth"
dir.create(out)
#validation cohorts
meth<-fread("../methyl/datasets/cd34/2020-05-25_methyl_data_before_limma.csv")
mtd<-fread("../methyl/datasets/cd34/cleaned_batch_CD34_library_date_220620.csv")

var_to_model<-c("Group_Sex","Mat.Age","latino","Group_Complexity_Fac")

mtdf<-mtd[rowSums(is.na(mtd[,.SD,.SDcols=var_to_model]))==0][batch==2]
formule<- ~0 + Group_Sex  + latino + Mat.Age + Group_Complexity_Fac 

design<-model.matrix(formule,data = data.frame(mtdf,row.names = "sample"))
fit <- lmFit(data.frame(meth,row.names = "locisID")[,mtdf$sample], design)

cont.matrix <- makeContrasts(C.L = "(Group_SexC_F+Group_SexC_M)-(Group_SexL_F+Group_SexL_M)",
                             F.M="(Group_SexC_F+Group_SexL_F)-(Group_SexC_M+Group_SexL_M)",
                             CF.CM="Group_SexC_F-Group_SexC_M",
                             CF.LF="Group_SexC_F-Group_SexL_F",
                             CF.LM="Group_SexC_F-Group_SexL_M",
                             CM.LM="Group_SexC_M-Group_SexL_M",
                             CM.LF="Group_SexC_M-Group_SexL_F",
                             LM.LF="Group_SexL_M-Group_SexL_F",
                             levels=design)


fit2  <- contrasts.fit(fit, cont.matrix)
fit2  <- eBayes(fit2)

res<-Reduce(rbind,lapply(colnames(cont.matrix), function(comp)data.table(topTable(fit2,coef = comp,n = Inf),keep.rownames = "cpg_id")[,compa:=comp]))
fwrite(res,fp(out,"res_limma_cohort2.tsv.gz"),sep="\t")
res<-fread(fp(out,"res_limma_cohort2.tsv.gz"),sep="\t")

table(res[adj.P.Val<0.2&abs(logFC)>30][,hyper_meth:=logFC>0][,.(hyper_meth,compa)])
p<-ggplot(res[compa%in%c("C.L","CF.LF",'CM.LM',"LM.LF")])+
  geom_point(aes(x=logFC,y=-log10(P.Value),col=adj.P.Val<0.2&abs(logFC)>30))+
  facet_wrap("compa")+
  scale_color_manual(values = c("grey","red"))
ggsave(fp(out,"volcano_plot_cohort2.png"),plot=p,height=5,width=7)

#comp to cohort 1 
mtdf<-mtd[rowSums(is.na(mtd[,.SD,.SDcols=var_to_model]))==0][batch==1]
formule<- ~0 + Group_Sex  + latino + Mat.Age + Group_Complexity_Fac 

design<-model.matrix(formule,data = data.frame(mtdf,row.names = "sample"))
fit <- lmFit(data.frame(meth,row.names = "locisID")[,mtdf$sample], design)

cont.matrix <- makeContrasts(C.L = "(Group_SexC_F+Group_SexC_M)-(Group_SexL_F+Group_SexL_M)",
                             F.M="(Group_SexC_F+Group_SexL_F)-(Group_SexC_M+Group_SexL_M)",
                             CF.CM="Group_SexC_F-Group_SexC_M",
                             CF.LF="Group_SexC_F-Group_SexL_F",
                             CF.LM="Group_SexC_F-Group_SexL_M",
                             CM.LM="Group_SexC_M-Group_SexL_M",
                             CM.LF="Group_SexC_M-Group_SexL_F",
                             LM.LF="Group_SexL_M-Group_SexL_F",
                             levels=design)


fit2  <- contrasts.fit(fit, cont.matrix)
fit2  <- eBayes(fit2)

#merge res cohort 1 and 2

res<-rbind(res[,cohort:=2],
           Reduce(rbind,lapply(colnames(cont.matrix),
                               function(comp)data.table(topTable(fit2,coef = comp,n = Inf),
                                                        keep.rownames = "cpg_id")[,compa:=comp]))[,cohort:=1])
fwrite(res,fp(out,"res_limma_cohorts.tsv.gz"),sep="\t")

#cl
ggplot(res[compa%in%c("C.L")])+
  geom_point(aes(x=logFC,y=-log10(P.Value),col=P.Value<10^-4&abs(logFC)>30),size=)+
  facet_wrap("cohort")+
  scale_color_manual(values = c("grey","red"))

#cflf
ggplot(res[compa%in%c("CF.LF")])+
  geom_point(aes(x=logFC,y=-log10(P.Value),col=P.Value<10^-4&abs(logFC)>30))+
  facet_wrap("cohort")+
  scale_color_manual(values = c("grey","red"))

#cmlm
ggplot(res[compa%in%c("CM.LM")])+
  geom_point(aes(x=logFC,y=-log10(P.Value),col=P.Value<10^-4&abs(logFC)>30))+
  facet_wrap("cohort")+
  scale_color_manual(values = c("grey","red"))


table(res[P.Value<10^-4&abs(logFC)>30][,hyper_meth:=logFC>0][,.(hyper_meth,compa,cohort)])

#genescore heatmap
meth_df<-melt(meth,id.vars = 1,
              variable.name = "sample",
              value.name = "methylation")
cpgs_weight<-fread("../methyl/ref/2020-06-29_All_CpG-Gene_links.csv")
cpgs_weight[,cpg_weight:=RegWeight+LinksWeight]
meth_df<-merge(meth_df,cpgs_weight[,.(locisID,gene,cpg_weight)],by="locisID")
meth_df[,gene_meth:=sum(methylation*cpg_weight)/sum(cpg_weight),by=c("sample","gene")]
fwrite(meth_df,fp(out,"methylation_sum_by_gene_datatable.tsv.gz"),sep="\t")
meth_genes<-unique(meth_df,by=c("sample","gene"))

meth_genes[,gene_meth_scaled:=scale(gene_meth),by=c("sample")]

meth_genes_mat<-dcast(meth_genes[,-c("locisID","hi.conf.genes","methylation","gene_meth")],sample~gene)
meth_genes_mat[1:10,1:10]

meth_genes_mat<-as.matrix(data.frame(meth_genes_mat,row.names = "sample"))
meth_genes_mat<-t(meth_genes_mat)
meth_genes_mat[1:10,1:10]
dim(meth_genes_mat)


pheatmap(meth_genes_mat,
         annotation_col = data.frame(mtd,row.names = "sample")[colnames(meth_genes_mat),c("Group_Sex","Group_Complexity_Fac")],
         show_rownames = F,show_colnames = F
        ) #doesn work
#see other test on figures_paper_methyl_annexe

#see key genes

ggplot(meth_genes[gene%in%c("SOCS3","HES1","JUN")&Group_name%in%c("C","L")])+
geom_boxplot(aes(x=Group_Sex,y=gene_meth))+facet_wrap("gene")

#res
res<-fread("../methyl/outputs/model14_without_iugr/2020-10-08_all_res_with_perm.csv")
resg<-unique(res[order(gene,pval)],by="gene")
ggplot(resg)+geom_density(aes(x=GeneScore))+facet_wrap("compa")

resg[,gs_scaled:=scale(GeneScore,center=F),by="compa"]
ggplot(resg)+geom_density(aes(x=gs_scaled))+facet_wrap("compa")

p1<-ggplot(resg[compa%in%c("C.L","CF.LF","CM.LM")],aes(x=compa,y=gs_scaled))+
  geom_jitter(aes(col=abs(gs_scaled)>1))+
  scale_color_manual(values = c("grey","red"))+
  geom_boxplot(aes(fill=compa),outlier.shape = NA)

p2<-ggplot(resg[compa%in%c("C.L","CF.LF","CM.LM")][abs(gs_scaled)>1][,hyper_meth:=gs_scaled>0])+geom_bar(aes(x=compa,fill=hyper_meth))
p1+p2


#Integration basal
source("../methyl/scripts/utils/new_utils.R")
library(Seurat)

out<-"../methyl/outputs/figures_paper_meth"

  #Lineage bias
cbps<-readRDS("../singlecell/outputs/02-hematopo_datasets_integration/cbps0-8_clean/cbps0-8_clean.rds")
mtd<-data.table(cbps@meta.data,keep.rownames ="bc")
mtd[,n.sample:=.N,by=c("sample","orig.ident")]
mtd[,pct.ct:=.N/n.sample,by=c("sample","orig.ident","cell_type")]
mtd[,pct.lin:=.N/n.sample,by=c("sample","orig.ident","lineage")]

cbps_cl<-subset(cbps,hto==F&group%in%c("ctrl","lga")&ambigous==F&orig.ident!="cd34_hto1_0C1I1L")
mtd_cl<-mtd[hto==F&group%in%c("ctrl","lga")&ambigous==F&orig.ident!="cd34_hto1_0C1I1L"]



ggplot(unique(mtd_cl[!str_detect(lineage,"unknown")][,.SD,.SDcols=unique(colnames(mtd))],by=c("sample","orig.ident","lineage")))+
  geom_boxplot(aes(x=group,y=pct.lin,fill=group))+facet_wrap("lineage",scales = "free")



  #DEGs
degs_cl<-fread("outputs/08-DEGs_LGA_no_stress/pseudobulk_deseq2_all_cbps/res_de_analysis_all_genes.csv")

ggplot(degs_cl,aes(x=log2FoldChange,y=-log10(pvalue),col=padj<0.1&abs(log2FoldChange)>0.6))+
  geom_point()+ 
  geom_label_repel(aes(label = ifelse(padj<0.1&
                                        abs(log2FoldChange)>0.6,gene,"")),
                   max.overlaps=3000,
                   box.padding   = 0.35,
                   point.padding = 0.5,
                   segment.color = 'grey50')+
  scale_color_manual(values = c("grey","red")) +
  theme_minimal() +
  theme(legend.position = "bottom")

#by lineage
degs_lin<-fread("../singlecell/outputs/08-DEGs_LGA_no_stress/pseudobulk_deseq2_by_lin/res_de_analysis_all_genes.csv")
ggplot(degs_lin,aes(x=log2FoldChange,y=-log10(pvalue),col=padj<0.1&abs(log2FoldChange)>0.6))+
  geom_point()+ 
  geom_label_repel(aes(label = ifelse(padj<0.1&
                                        abs(log2FoldChange)>0.6,gene,"")),
                   max.overlaps=3000,
                   box.padding   = 0.35,
                   point.padding = 0.5,
                   segment.color = 'grey50')+
  scale_color_manual(values = c("grey","red")) +
  facet_wrap("lineage")+
  theme_minimal() +
  theme(legend.position = "bottom")

#single cell degs

degs_cl_sc<-fread("outputs/08-DEGs_LGA_no_stress/sc_edger_deseq2_all_cbps/res_de_analysis_all_genes.csv")

degs_cl_sc[p_val_adj<0.01&abs(avg_logFC)>0.6]

ggplot(degs_cl_sc,aes(x=avg_logFC,y=-log10(p_val),col=p_val_adj<0.1&abs(avg_logFC)>0.6))+
  geom_point()+
  scale_color_manual(values = c("grey","red")) +
  theme_minimal() +
  theme(legend.position = "bottom")



#correl scdegs with methyl

res_cl<-fread("../methyl/outputs/model14_without_iugr/2020-09-16_res_C.L_with_GeneScore_and_permut.csv")

res_cl[gene=="HES1"]

res_cl_merge<-merge(degs_cl_sc,res_cl)
res_cl_merge<-unique(res_cl_merge,by="gene")
ggplot(res_cl_merge) + geom_point(aes(x=-log10(p_val_adj)*avg_logFC,y=GeneScore))
ggplot(res_cl_merge) + geom_point(aes(x=avg_logFC,y=GeneScore))

#Over repre test GO BP and KEGG pathway

res_cl_merge[,gs_scaled:=scale(GeneScore,center = F)]
ggplot(res_cl_merge)+geom_density(aes(x=gs_scaled))
ggplot(res_cl_merge)+geom_density(aes(x=avg_logFC))

res_cl_merge[,correl:=ifelse(gs_scaled>1&avg_logFC>1,"Hypermet-Upreg",
                             ifelse(gs_scaled>1&avg_logFC<(-1),"Hypermet-Downreg",
                              ifelse(gs_scaled<(-1)&avg_logFC<(-1),"Hypomet-Downreg",
                               ifelse(gs_scaled<(-1)&avg_logFC>1,"Hypomet-Upreg","no change"))))]
ggplot(res_cl_merge) + geom_point(aes(x=avg_logFC,y=gs_scaled,col=correl))

#kegg
pathways<-read.gmt("../methyl/ref/msigdb/c2.cp.kegg.v7.1.symbols.gmt")
pathways_list<-split(pathways$gene,f = pathways$term)

res_kegg<-over_repr_test_multi(genes_of_interest =res_cl_merge[correl=="Hypermet-Upreg"]$gene ,
                     terms_list = pathways_list ,
                     size_universe = nrow(res_cl_merge)
                    )

res_kegg[padj<0.5]#0

#with hypermet-downreg
res_kegg2<-over_repr_test_multi(genes_of_interest =res_cl_merge[correl=="Hypermet-Downreg"]$gene ,
                     terms_list = pathways_list ,
                     size_universe = nrow(res_cl_merge)
                    )

res_kegg2[padj<0.5]#0


#GO
go_term<-read.gmt("../methyl/ref/msigdb/c5.go.bp.v7.3.symbols.gmt")
go_list<-split(go_term$gene,f = go_term$term)

res_go<-over_repr_test_multi(genes_of_interest =res_cl_merge[correl=="Hypermet-Upreg"]$gene ,
                     terms_list = go_list ,
                     size_universe = nrow(res_cl_merge)
                    )

res_go[padj<0.5]#0
#with hypermet-downreg
res_go<-over_repr_test_multi(genes_of_interest =res_cl_merge[correl=="Hypermet-Downreg"]$gene ,
                     terms_list = go_list ,
                     size_universe = nrow(res_cl_merge)
                    )

res_go[padj<0.5]#0


#GSEA

summary(lm(res_cl_merge$avg_logFC~res_cl_merge$GeneScore))

