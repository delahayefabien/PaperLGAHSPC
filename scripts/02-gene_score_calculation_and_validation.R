
#Gene score calculation and validation
source("scripts/utils/new_utils.R")
source("scripts/utils/methyl_utils.R")
library(limma)
out<-"outputs/02-gene_score_calculation_and_validation"
dir.create(out)

res<-fread("outputs/01-lga_vs_ctrl_limma_DMCs_analysis/res_limma.tsv.gz",sep="\t",
           select = c("cpg_id","P.Value","adj.P.Val","AveExpr","logFC"),
           col.names = c("cpg_id","pval","padj","avg.meth","meth.change"))

cpgs_ref<-fread("ref/2020-06-29_All_CpG-Gene_links.csv")
cpgs_ref<-cpgs_ref[,cpg_id:=locisID][,-"locisID"]

res_anno<-merge(res,cpgs_ref,by="cpg_id")
res_anno[,cpg_score:=(-log10(pval)/4*meth.change)*RegWeight*LinksWeight]


res_anno[,n.cpg_weight:=(1/sum(1/(abs(cpg_score)+1)))^(1/4),by="gene"]
res_anno[,gene_score:=sum(cpg_score)*n.cpg_weight,by="gene"]

res_anno[,region_type:=ifelse(abs(tss_dist)<=2000,"promoter","other"),by="gene"]
res_anno[,n_cpg_weight_region:=(1/sum(1/(abs(cpg_score)+1)))^(1/4),by=c('region_type',"gene")]
res_anno[,gene_score_prom:=sum(cpg_score[region_type=="promoter"])*n_cpg_weight_region,by=c("gene")]
res_anno[,gene_score_enh:=sum(abs(cpg_score[region_type=="other"]))*n_cpg_weight_region,by=c("gene")]

res_anno[is.na(gene_score_prom),gene_score_add:=gene_score_enh]
res_anno[is.na(gene_score_enh),gene_score_add:=gene_score_prom]
res_anno[!(is.na(gene_score_prom)|is.na(gene_score_enh)),gene_score_add:=gene_score_prom+gene_score_enh]

res_anno<-res_anno[!is.na(gene_score_add)] #rm 54 genes with no tss dist so can not calculate genescore

#add useful info :
res_anno[,n.cpg.gene:=.N,by=.(gene)]
res_anno[,n.cpg.gene.region:=.N,by=.(gene,region_type)]
res_anno[,avg.meth.change.region:=mean(meth.change),by=c('region_type',"gene")]
res_anno[,avg.m.log10.pval.region:=mean(-log10(pval)),by=c('region_type',"gene")]
res_anno[,avg.dmc_score.region:=mean(-log10(pval)*meth.change),by=c('region_type',"gene")]

res_anno[,avg.meth.change:=mean(meth.change),by=c("gene")]
res_anno[,avg.m.log10.pval:=mean(-log10(pval)),by=c("gene")]
res_anno[,avg.dmc_score:=mean(-log10(pval)*meth.change),by=c("gene")]


fwrite(res_anno,fp(out,"res_gene_score.tsv.gz"),sep="\t")

#VALIDATION  Gene Score
#check gene_score not correlated with ncpg.gene..
resg<-unique(res_anno[order(gene,region_type,pval)],by=c("gene"))

ggplot(resg)+
      geom_boxplot(aes(x = as.factor(n.cpg.gene),y =gene_score )) #ok

ggplot(resg)+
      geom_boxplot(aes(x = as.factor(n.cpg.gene),y =gene_score_add )) #ok

ggplot(resg)+
      geom_boxplot(aes(x = as.factor(n.cpg.gene.region),y =gene_score_prom )) #ok

ggplot(resg)+
      geom_boxplot(aes(x = as.factor(n.cpg.gene.region),y =gene_score_enh )) #ok

resg[gene_score_add>200]$gene
resg[gene=="SOCS3"]
resg[gene=="HES1"]

#..but correlated with ncpg sig
res_anno[,n.cpg.sig.gene:=sum(pval<0.01),by=.(gene)]
ggplot(unique(res_anno,by="gene"))+
      geom_boxplot(aes(x = as.factor(n.cpg.sig.gene),y =gene_score )) #ok
res_anno[n.cpg.sig.gene==4084] #non gene linked cpgs

res_anno[,n.cpg.sig.gene.region:=sum(pval<0.01),by=.(gene,region_type)]
ggplot(unique(res_anno[region_type=="promoter"],by="gene"))+
      geom_boxplot(aes(x = as.factor(n.cpg.sig.gene.region),y =gene_score_region )) #ok

ggplot(unique(res_anno[region_type=="other"],by="gene"))+
      geom_boxplot(aes(x = as.factor(n.cpg.sig.gene.region),y =gene_score_region )) #ok


#Compared with classical metrics
res_de_cl<-fread("../singlecell/outputs/04-DEG_in_LGA/2020-09-01_pseudo_bulk_DEseq2_LgaVsCtrl_CBP1andcbp558_559_samples_excluded_regr_on_batch_and_sex_all_genes.csv")

res_de_cl[is.na(padj),padj:=1]
resg_de<-merge(resg,res_de_cl,by=c("gene"))
resg_de[padj.y<0.1]
p1<-ggplot(resg_de)+geom_boxplot(aes(x=padj.y<0.1,y=gene_score_prom),outlier.shape = NA)+scale_y_continuous(limits = c(0,300))
p1.5<-ggplot(resg_de)+geom_boxplot(aes(x=padj.y<0.1,y=gene_score_enh))+scale_y_continuous(limits = c(0,100))
p2<-ggplot(unique(resg_de,by="gene"))+geom_boxplot(aes(x=padj.y<0.1,y=gene_score_add))+scale_y_continuous(limits = c(0,300))

p3<-ggplot(unique(resg_de,by="gene"))+geom_boxplot(aes(x=padj.y<0.1,y=gene_score))+scale_y_continuous(limits = c(0,300))

p1+p1.5+p2+p3

wilcox.test(resg_de[padj.y<=0.1]$gene_score_add,resg_de[padj.y>0.1]$gene_score_add)
#p=0.0004782

wilcox.test(resg_de[padj.y<=0.1]$gene_score,resg_de[padj.y>0.1]$gene_score)
#p=0.001755

 
p4<-ggplot(unique(resg_de,by="gene"))+geom_boxplot(aes(x=padj.y<0.1,y=avg.meth.change))

p5<-ggplot(unique(resg_de,by="gene"))+geom_boxplot(aes(x=padj.y<0.1,y=avg.m.log10.pval))

p6<-ggplot(unique(resg_de,by="gene"))+geom_boxplot(aes(x=padj.y<0.1,y=avg.dmc_score))

p2+p4+p5+p6

p7<-ggplot(unique(resg_de,by="gene"))+geom_boxplot(aes(x=padj.y<0.1,y=abs(meth.change.x)))

p5<-ggplot(unique(resg_de,by="gene"))+geom_boxplot(aes(x=padj.y<0.1,y=-log10(pval.x)))

p6<-ggplot(unique(resg_de,by="gene"))+geom_boxplot(aes(x=padj.y<0.1,y=-log10(pval.x)*abs(meth.change.x)))

p2+p4+p5+p6

res_de_cl2<-fread("../singlecell/outputs/08-DEGs_LGA_no_stress/pseudobulk_deseq2_all_cbps/res_de_analysis_all_genes.csv")

resg_de2<-merge(resg,res_de_cl2,by=c("gene"))
resg_de2[padj.y<0.5]
p1<-ggplot(resg_de2[region_type!=""])+geom_boxplot(aes(x=padj.y<0.05,y=gene_score_region))+facet_wrap("region_type")

p2<-ggplot(unique(resg_de2,by="gene"))+geom_boxplot(aes(x=padj.y<0.05,y=gene_score_region_add))

p3<-ggplot(unique(resg_de2,by="gene"))+geom_boxplot(aes(x=padj.y<0.05,y=GeneScore))

p1+p2+p3

wilcox.test(unique(resg_de2,by="gene")[padj.y<=0.05]$gene_score_region_add,unique(resg_de2,by="gene")[padj.y>0.05]$gene_score_region_add)
#p=0.008036

wilcox.test(unique(resg_de2,by="gene")[padj.y<=0.05]$gene_score,unique(resg_de2,by="gene")[padj.y>0.05]$gene_score)
#p=0.02
