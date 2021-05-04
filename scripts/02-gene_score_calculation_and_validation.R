
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

#link cpg to gene and weitgh link confidence
#see 02A 
cpgs_ref<-fread("ref/2020-06-29_All_CpG-Gene_links.csv")
max(abs(cpgs_ref[in_eQTR==F]$tss_dist))
cpgs_ref<-cpgs_ref[,cpg_id:=locisID][,-"locisID"]

res<-merge(res,cpgs_ref,by="cpg_id")

#add useful info :
res[,region_type:=ifelse(abs(tss_dist)<=2000,"promoter","other"),by="gene"]
res[,n.cpg.gene:=.N,by=.(gene)]
res[,n.cpg.gene.region:=.N,by=.(gene,region_type)]
res[,avg.meth.change.region:=mean(meth.change),by=c('region_type',"gene")]
res[,avg.m.log10.pval.region:=mean(-log10(pval)),by=c('region_type',"gene")]
res[,avg.dmc_score.region:=mean(-log10(pval)*meth.change),by=c('region_type',"gene")]

res[,avg.meth.change:=mean(meth.change),by=c("gene")]
res[,avg.m.log10.pval:=mean(-log10(pval)),by=c("gene")]
res[,avg.dmc_score:=mean(-log10(pval)*meth.change),by=c("gene")]

#GeneScore calculation
n_min_group<-min(table(mtd$group))
  #first, calculate CpGScore :
res[,cpg_score:=(-log10(pval)/sqrt(n_min_group)*meth.change)*RegWeight*LinksWeight] #divided by n_sample ro normalized gene score ~ n_sample


  #then, the GeneScore :
res[,n.cpg_weight:=(1/sum(1/(abs(cpg_score)+1)))^(1/4),by="gene"] #n.cpg_weight to reduce the influence of the n_cpg by gene to the GeneScore
ggplot(unique(res,by="gene"))+geom_point(aes(x=n.cpg.gene,y=n.cpg_weight))
res[,gene_score:=sum(cpg_score)*n.cpg_weight,by="gene"] 

ggplot(unique(res,by="gene"))+geom_boxplot(aes(x = as.factor(n.cpg.gene),y =gene_score )) #ok

# the GeneScore by region :
res[,n_cpg_weight_region:=(1/sum(1/(abs(cpg_score)+1)))^(1/4),by=c('region_type',"gene")]
res[region_type=="promoter",gene_score_region:=sum(cpg_score)*n_cpg_weight_region,by=c("gene")]
res[region_type=="other",gene_score_region:=sum(abs(cpg_score))*n_cpg_weight_region,by=c("gene")]

res[,gene_score_add:=sum(unique(abs(gene_score_region)),na.rm = T),by="gene"]
res[is.na(gene_score_add)] 
unique(res[gene=="A1BG",.(gene,gene_score,region_type,gene_score_region,gene_score_add)] )
res<-res[!is.na(gene_score_add)] #rm 54 genes with no tss dist so can not calculate genescore
ggplot(unique(res,by="gene"))+geom_boxplot(aes(x = as.factor(n.cpg.gene),y =gene_score_add )) #ok

#not correlated to ncpg.gene, but crrelated to n cpg sig :
res[,n.cpg.sig.gene:=sum(pval<0.01),by=.(gene)]
plot(density(unique(res,by="gene")$gene_score_add))
abline(v=70)
unique(res,by="gene")[gene_score_add>70]

ggplot(unique(res,by="gene"))+
      geom_boxplot(aes(x = as.factor(n.cpg.sig.gene),y =gene_score )) #ok
ggplot(unique(res,by="gene"))+
      geom_boxplot(aes(x = as.factor(n.cpg.sig.gene),y =gene_score_add )) #ok


#VALIDATION  Gene Score
#valid wieight
source("scripts/utils/new_utils.R")

res<-fread("outputs/02-gene_score_calculation_and_validation/res_gene_score.tsv.gz")
res[,ncpg.sig.gene:=sum(pval<0.01),by="gene"]
resg<-unique(res[order(gene,region_type,pval)],by="gene")
summary(resg$gene_score_add)
#see correl covariates with genescore
resg
summary(lm(gene_score~n.cpg.gene+ncpg.sig.gene+pval+meth.change+type+EnsRegScore+in_eQTR+abs(tss_dist),data = resg))

summary(lm(gene_score_add~n.cpg.gene+ncpg.sig.gene+pval+meth.change+type+EnsRegScore+in_eQTR+abs(tss_dist),data = resg)) #best gene_score_add than gene_score

summary(lm(gene_score_add~n.cpg.gene,data = resg))$r.squared


#deter if ^1/4 for n_cpg_weight were the opti one to reduce influence of n_cpg_gene (werease kept influences of n_cpg_sig)
res_test<-copy(res)
rs_ncpg<-rep(0,10)
rs_ncpg_sig<-rep(0,10)
ps<-list()
for(i in 1:10){
  print(i)
  res_test[,n_cpg_weight_region:=(1/sum(1/(abs(cpg_score)+1)))^(1/i),by=c('region_type',"gene")]
  res_test[region_type=="promoter",gene_score_region:=sum(cpg_score)*n_cpg_weight_region,by=c("gene")]
  res_test[region_type=="other",gene_score_region:=sum(abs(cpg_score))*n_cpg_weight_region,by=c("gene")]
  
  res_test[,gene_score_add:=sum(unique(abs(gene_score_region)),na.rm = T),by="gene"]
  resg_test<-unique(res_test[order(gene,region_type,pval)],by="gene")
  print(summary(lm(gene_score_add~n.cpg.gene+ncpg.sig.gene+pval+meth.change+type+EnsRegScore+in_eQTR+abs(tss_dist),data = resg_test))) 

  rs_ncpg[i]<-summary(lm(gene_score_add~n.cpg.gene,data = resg_test))$r.squared
  rs_ncpg_sig[i]<-summary(lm(gene_score_add~ncpg.sig.gene,data = resg_test))$r.squared
  ps[[i]]<-ggplot(resg_test)+geom_boxplot(aes(x = as.factor(n.cpg.gene),y =gene_score_add )) +ggtitle(ps("x = ",i))


}
wrap_plots(ps)

plot(1:10,rs_ncpg)
plot(1:10,rs_ncpg_sig)

rs_ncpg2<-rep(0,10)
rs_ncpg_sig2<-rep(0,10)
xs<-3+(1:10/10)
for(i in 1:length(xs)){
  x<-xs[i]
  print(x)
  res_test[,n_cpg_weight_region:=(1/sum(1/(abs(cpg_score)+1)))^(1/x),by=c('region_type',"gene")]
  res_test[region_type=="promoter",gene_score_region:=sum(cpg_score)*n_cpg_weight_region,by=c("gene")]
  res_test[region_type=="other",gene_score_region:=sum(abs(cpg_score))*n_cpg_weight_region,by=c("gene")]
  
  res_test[,gene_score_add:=sum(unique(abs(gene_score_region)),na.rm = T),by="gene"]
  resg_test<-unique(res_test[order(gene,region_type,pval)],by="gene")
  print(summary(lm(gene_score_add~n.cpg.gene+ncpg.sig.gene+pval+meth.change+type+EnsRegScore+in_eQTR+abs(tss_dist),data = resg_test))) 

  rs_ncpg2[i]<-summary(lm(gene_score_add~n.cpg.gene,data = resg_test))$r.squared
  rs_ncpg_sig2[i]<-summary(lm(gene_score_add~ncpg.sig.gene,data = resg_test))$r.squared


}
plot(xs,rs_ncpg2)
plot(xs,rs_ncpg_sig2)
#3.4 is best : 
# [1] 3.4
# 
# Call:
# lm(formula = gene_score_add ~ n.cpg.gene + ncpg.sig.gene + pval + 
#     meth.change + type + EnsRegScore + in_eQTR + abs(tss_dist), 
#     data = resg_test)
# 
# Residuals:
#     Min      1Q  Median      3Q     Max 
# -421.86  -14.97   -6.11   11.93  202.40 
# 
# Coefficients:
#                 Estimate Std. Error t value Pr(>|t|)    
# (Intercept)    1.108e+01  6.097e-01  18.176  < 2e-16 ***
# n.cpg.gene    -9.621e-04  8.930e-03  -0.108   0.9142    
# ncpg.sig.gene  2.334e+01  1.250e-01 186.726  < 2e-16 ***
# pval          -8.026e+00  1.084e+00  -7.405 1.37e-13 ***
# meth.change    6.203e-02  1.142e-02   5.434 5.58e-08 ***
# type           7.337e-01  9.591e-02   7.650 2.09e-14 ***
# EnsRegScore    1.452e+01  1.074e+00  13.520  < 2e-16 ***
# in_eQTRTRUE   -1.304e+00  5.737e-01  -2.273   0.0231 *  
# abs(tss_dist) -4.914e-06  1.160e-06  -4.237 2.28e-05 ***
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 27.83 on 21189 degrees of freedom
#   (54 observations deleted due to missingness)
# Multiple R-squared:  0.728,	Adjusted R-squared:  0.7279 
# F-statistic:  7090 on 8 and 21189 DF,  p-value: < 2.2e-16

res[,n_cpg_weight_region:=(1/sum(1/(abs(cpg_score)+1)))^(1/3.4),by=c('region_type',"gene")]
res[region_type=="promoter",gene_score_region:=sum(cpg_score)*n_cpg_weight_region,by=c("gene")]
res[region_type=="other",gene_score_region:=sum(abs(cpg_score))*n_cpg_weight_region,by=c("gene")]

res[,gene_score_add:=sum(unique(abs(gene_score_region)),na.rm = T),by="gene"]

fwrite(res,fp(out,"res_gene_score.tsv.gz"),sep="\t")


resg<-unique(res[order(gene,region_type,pval)],by="gene")


res_de_cl<-fread("../singlecell/outputs/04-DEG_in_LGA/2020-09-01_pseudo_bulk_DEseq2_LgaVsCtrl_CBP1andcbp558_559_samples_excluded_regr_on_batch_and_sex_all_genes.csv")

res_de_cl[is.na(padj),padj:=1]
resg_de<-merge(resg,res_de_cl,by=c("gene"))
resg_de[padj.y<0.05]

p2<-ggplot(resg_de)+geom_boxplot(aes(x=padj.y<0.1,y=gene_score_add))

p3<-ggplot(resg_de)+geom_boxplot(aes(x=padj.y<0.1,y=gene_score))

p2+p3

wilcox.test(resg_de[padj.y<=0.1]$gene_score_add,resg_de[padj.y>0.1]$gene_score_add)
#p=0.0009198

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
#p=0.001231



