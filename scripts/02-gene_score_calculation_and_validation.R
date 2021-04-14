
#Gene score calculation and validation
source("scripts/utils/new_utils.R")
source("scripts/utils/methyl_utils.R")
library(limma)
out<-"outputs/02-gene_score_calculation_and_validation"
dir.create(out)

res<-fread(fp(out,"res_limma.tsv.gz"),sep="\t",
           select = c("cpg_id","P.Value","adj.P.Val","AveExpr","logFC"),
           col.names = c("cpg_id","pval","padj","avg.meth","meth.change"))

cpgs_ref<-fread("ref/2020-06-29_All_CpG-Gene_links.csv")
cpgs_ref<-cpgs_ref[,cpg_id:=locisID][,-"locisID"]

res<-merge(res,cpgs_ref,all.x=T,by="cpg_id")
res[,cpg_score:=(-log10(pval)/4*meth.change)*RegWeight*LinksWeight]


res[,n.cpg_weight:=(1/sum(1/(abs(cpg_score)+1)))^(1/4),by="gene"]
res[,gene_score:=sum(cpg_score)*n.cpg_weight,by="gene"]

res[,region_type:=ifelse(abs(tss_dist)<=2000,"promoter","other"),by="gene"]
res[,n_cpg_weight_region:=(1/sum(1/(abs(cpg_score)+1)))^(1/4),by=c('region_type',"gene")]
res[region_type=="promoter",gene_score_region:=sum(cpg_score)*n_cpg_weight_region,by=c("gene")]
res[region_type=="other",gene_score_region:=sum(abs(cpg_score))*n_cpg_weight_region,by=c("gene")]


#check gene_score not correlated with ncpg.gene..
res[,n.cpg.gene:=.N,by=.(gene)]

ggplot(unique(res,by="gene"))+
      geom_boxplot(aes(x = as.factor(n.cpg.gene),y =gene_score )) #ok

res[,n.cpg.gene.region:=.N,by=.(gene,region_type)]
ggplot(unique(res[region_type=="promoter"],by="gene"))+
      geom_boxplot(aes(x = as.factor(n.cpg.gene.region),y =gene_score_region )) #ok

ggplot(unique(res[region_type=="other"],by="gene"))+
      geom_boxplot(aes(x = as.factor(n.cpg.gene.region),y =gene_score_region )) #ok

unique(res[region_type=="promoter"][gene_score_region>150]$gene)
unique(res,by=c('region_type',"gene"))[gene=="SOCS3"]
unique(res,by=c('region_type',"gene"))[gene=="HES1"]

#..but correlated with ncpg sig
res[,n.cpg.sig.gene:=sum(pval<0.01),by=.(gene)]
ggplot(unique(res,by="gene"))+
      geom_boxplot(aes(x = as.factor(n.cpg.sig.gene),y =gene_score )) #ok
res[n.cpg.sig.gene==4084] #non gene linked cpgs

res[,n.cpg.sig.gene.region:=sum(pval<0.01),by=.(gene,region_type)]
ggplot(unique(res[region_type=="promoter"],by="gene"))+
      geom_boxplot(aes(x = as.factor(n.cpg.sig.gene.region),y =gene_score_region )) #ok

ggplot(unique(res[region_type=="other"],by="gene"))+
      geom_boxplot(aes(x = as.factor(n.cpg.sig.gene.region),y =gene_score_region )) #ok


length(unique(res[gene_score>70]$gene))

fwrite(res,fp(out,"res_gene_score.tsv.gz"),sep="\t")


#VALIDATION  Gene Score
#Compared with classical metrics
res[,avg.meth.change:=mean(meth.change),by=c('region_type',"gene")]
res[,avg.m.log10.pval:=mean(-log10(pval)),by=c('region_type',"gene")]
res[,avg.dmc_score:=mean(-log10(pval)*meth.change),by=c('region_type',"gene")]


#DMR 


#

