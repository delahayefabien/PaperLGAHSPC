
#Gene score calculation and validation
source("scripts/utils/new_utils.R")
source("scripts/utils/methyl_utils.R")
library(limma)
out<-"outputs/02-gene_score_calculation_and_validation"
dir.create(out)


res<-fread("outputs/01-lga_vs_ctrl_limma_DMCs_analysis/res_limma.tsv.gz",sep="\t",
           select = c("cpg_id","P.Value","adj.P.Val","AveExpr","logFC"),
           col.names = c("cpg_id","pval","padj","avg.meth","meth.change"))
mtd<-fread("datasets/cd34/metadata_pcs_cl_190421.csv",sep=";")


cpgs_ref<-fread("ref/2020-06-29_All_CpG-Gene_links.csv")
cpgs_ref<-cpgs_ref[,cpg_id:=locisID][,-"locisID"]

res_anno<-merge(res,cpgs_ref,by="cpg_id")

#add useful info :
res_anno[,n.cpg.gene:=.N,by=.(gene)]
res_anno[,n.cpg.gene.region:=.N,by=.(gene,region_type)]
res_anno[,avg.meth.change.region:=mean(meth.change),by=c('region_type',"gene")]
res_anno[,avg.m.log10.pval.region:=mean(-log10(pval)),by=c('region_type',"gene")]
res_anno[,avg.dmc_score.region:=mean(-log10(pval)*meth.change),by=c('region_type',"gene")]

res_anno[,avg.meth.change:=mean(meth.change),by=c("gene")]
res_anno[,avg.m.log10.pval:=mean(-log10(pval)),by=c("gene")]
res_anno[,avg.dmc_score:=mean(-log10(pval)*meth.change),by=c("gene")]


n_min_group<-min(table(mtd$group))
res_anno[,cpg_score:=(-log10(pval)/sqrt(n_min_group)*meth.change)*RegWeight*LinksWeight]


res_anno[,n.cpg_weight:=(1/sum(1/(abs(cpg_score)+1)))^(1/4),by="gene"]
ggplot(unique(res_anno,by="gene"))+geom_point(aes(x=n.cpg.gene,y=n.cpg_weight))
res_anno[,gene_score:=sum(cpg_score)*n.cpg_weight,by="gene"]
ggplot(unique(res_anno,by="gene"))+geom_boxplot(aes(x = as.factor(n.cpg.gene),y =gene_score )) #ok

res_anno[,region_type:=ifelse(abs(tss_dist)<=2000,"promoter","other"),by="gene"]
res_anno[,n_cpg_weight_region:=(1/sum(1/(abs(cpg_score)+1)))^(1/4),by=c('region_type',"gene")]
res_anno[,gene_score_region:=sum(cpg_score)*n_cpg_weight_region,by=c("gene","region_type")]

res_anno[,gene_score_add:=sum(unique(gene_score_region),na.rm = T),by="gene"]
res_anno[is.na(gene_score_add)] 
res_anno<-res_anno[!is.na(gene_score_add)] #rm 54 genes with no tss dist so can not calculate genescore
ggplot(unique(res_anno,by="gene"))+geom_boxplot(aes(x = as.factor(n.cpg.gene),y =gene_score_add )) #ok

#not correlated to ncpg.gene, bu tocrrelated to n cpg sig :
res_anno[,n.cpg.sig.gene:=sum(pval<0.01),by=.(gene)]
ggplot(unique(res_anno,by="gene"))+
      geom_boxplot(aes(x = as.factor(n.cpg.sig.gene),y =gene_score )) #ok
ggplot(unique(res_anno,by="gene"))+
      geom_boxplot(aes(x = as.factor(n.cpg.sig.gene),y =gene_score_add )) #ok


fwrite(res_anno,fp(out,"res_gene_score.tsv.gz"),sep="\t")

#VALIDATION  Gene Score
resg<-unique(res_anno[order(gene,region_type,pval)],by=c("gene"))

res_de_cl<-fread("../singlecell/outputs/04-DEG_in_LGA/2020-09-01_pseudo_bulk_DEseq2_LgaVsCtrl_CBP1andcbp558_559_samples_excluded_regr_on_batch_and_sex_all_genes.csv")

res_de_cl[is.na(padj),padj:=1]
resg_de<-merge(resg,res_de_cl,by=c("gene"))
resg_de[padj.y<0.1]

p2<-ggplot(resg_de)+geom_boxplot(aes(x=padj.y<0.1,y=gene_score_add))

p3<-ggplot(resg_de)+geom_boxplot(aes(x=padj.y<0.1,y=gene_score))

p2+p3

wilcox.test(resg_de[padj.y<=0.1]$gene_score_add,resg_de[padj.y>0.1]$gene_score_add)
#p=0.0012

wilcox.test(resg_de[padj.y<=0.1]$gene_score,resg_de[padj.y>0.1]$gene_score)
#p=0.0015

 
p4<-ggplot(unique(resg_de,by="gene"))+geom_boxplot(aes(x=padj.y<0.1,y=avg.meth.change))

p5<-ggplot(unique(resg_de,by="gene"))+geom_boxplot(aes(x=padj.y<0.1,y=avg.m.log10.pval))

p6<-ggplot(unique(resg_de,by="gene"))+geom_boxplot(aes(x=padj.y<0.1,y=avg.dmc_score))

p2+p4+p5+p6

p7<-ggplot(unique(resg_de,by="gene"))+geom_boxplot(aes(x=padj.y<0.1,y=abs(meth.change.x)))

p8<-ggplot(unique(resg_de,by="gene"))+geom_boxplot(aes(x=padj.y<0.1,y=-log10(pval.x)))

p9<-ggplot(unique(resg_de,by="gene"))+geom_boxplot(aes(x=padj.y<0.1,y=-log10(pval.x)*abs(meth.change.x)))

p2+p7+p8+p9

wilcox.test(resg_de[padj.y<=0.1]$meth.change.x,resg_de[padj.y>0.1]$meth.change.x)
#p=0.04

res_cl<-fread("outputs/model14_without_iugr/2020-09-16_res_C.L_with_GeneScore_and_permut.csv")
res_cl_merge<-merge(res_cl,res_de_cl[,.(gene,log2FoldChange,pvalue,padj)])

p11<-ggplot(unique(res_cl_merge,by="gene"))+geom_boxplot(aes(x=padj<0.1,y=GeneScore))

p2<-p2+coord_cartesian(ylim = c(0,200))
p11<-p11+coord_cartesian(ylim = c(0,200))
p2+p11

wilcox.test(unique(res_cl_merge,by="gene")[padj<=0.1]$GeneScore,unique(res_cl_merge,by="gene")[padj>0.1]$GeneScore)
#p=0.0012
