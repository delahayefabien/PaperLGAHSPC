#genescore by prom/ enh, best correl with Foldchange of DEGs ?
source("scripts/utils/new_utils.R")

res_cl<-fread("outputs/model14_without_iugr/2020-09-16_res_C.L_with_GeneScore_and_permut.csv")


#gene_score calculation
res_cl[,region_type:=ifelse(abs(tss_dist)<=2000,"promoter","other"),by="gene"]

res_cl[,n_cpg_weight:=(1/sum(1/(abs(CpGScore)+1)))^(1/4),by=c('region_type',"gene")]
  

res_cl[,gene_score:=sum(CpGScore)*n_cpg_weight,by=c('region_type',"gene")]

#r2 avec DEGs (top1000)

res_de_cl<-fread("../singlecell/outputs/08-DEGs_LGA_no_stress/pseudobulk_deseq2_all_cbps/res_de_analysis_all_genes.csv")
res_de_cl[,top1000:=pvalue<=sort(pvalue)[1000]]

res_cl_1000<-merge(unique(res_cl[order(pval)],by=c("gene","region_type")),res_de_cl[top1000==T],by=c("gene"))
res_cl_1000[,r2:=cor(abs(gene_score),abs(log2FoldChange))^2,by=c("region_type")]

res_cl_1000
