#TF signature enriched in DMC  ?
source("../methyl/scripts/utils/new_utils.R")
out<-here("../methyl/outputs/05-regulons_enrichment_genescore")
dir.create(out)
regulons_list<-readRDS("../singlecell/outputs/05-SCENIC/cbps0-8_clean/regulons_list.rds")

#in GeneScore
res_methg<-fread("../methyl/outputs/02-gene_score_calculation_and_validation/res_genes.csv.gz")

#GSEA
library(fgsea)

res_methg[,gs_rank:=rank(gene_score_add)]

genes_rank<-res_methg$gs_rank
names(genes_rank)<-res_methg$gene
res_gsea<-fgsea(pathways=regulons_list,
      stats=genes_rank,eps=0)

fwrite(res_gsea,fp(out,"res_gsea_genescore_regulons.csv.gz"))
res_gsea[padj<0.001]
res_gsea[order(padj)]$pathway[1:10]
res_gsea[pathway=="STAT3e"]
res_gsea[pathway=="STAT3"]


#over representation test
res_methg[,gene_altered:=abs(gene_score_add)>150]

res_or<-data.table(regulon=names(regulons_list),regulon.size=sapply(regulons_list,length))
genes_altered<-res_methg[gene_altered==T]$gene
res_or[,n.genes.altered:=length(genes_altered)]
res_or[,n.enriched:=sum(genes_altered%in%regulons_list[[regulon]]),by="regulon"]
res_or[,genes.enriched:=paste(genes_altered[genes_altered%in%regulons_list[[regulon]]],collapse="|"),by="regulon"]
res_or[,pct.enriched:=n.enriched/regulon.size]

size_universe<-length(res_methg$gene)

res_or[,pval:=phyper(q=n.enriched-1, 
     m=n.genes.altered, 
     n=size_universe-n.genes.altered, 
     k=regulon.size, 
     lower.tail=FALSE),
     by="regulon"]

res_or[,padj:=p.adjust(pval,method = 'BH')]
  
fwrite(res_or,fp(out,"res_or_genescore150_regulons.csv.gz"))



#pca on the cpgs determining TF #[to update]

genes_determining_tf<-Reduce(union,regulons_list[unique(res_or_all[padj<0.001&!str_detect(regulon,"e$")]$regulon)])

res_meth[,tf_determining_genes:=gene%in%genes_determining_tf]

res_meth[,tf_determining_cpgs:=tf_determining_genes&pval<0.001&abs(meth.change)>25]

cpgs_determining_tf<-unique(res_meth[tf_determining_cpgs==T]$locisID)
length(cpgs_determining_tf)
meth<-fread("../methyl/datasets/cd34/2020-05-25_methyl_data_before_limma.csv")
meth_deter<-meth[locisID%in%cpgs_determining_tf] #1413 cpgs

pca<-prcomp(t(as.matrix(meth_deter[,.SD,.SDcols=mtd[Group_name%in%c("C","L")]$sample])))

pca$x
vars_pcs<-GetVarPCs(pca)
round(vars_pcs*100)

res_pca_deter<-as.data.table(pca$x,keep.rownames = "sample")

mtd<-fread("../methyl/datasets/cd34/cleaned_batch_CD34_250121.csv")

res_pca_deter<-merge(res_pca_deter,mtd,by="sample")

ggplot(res_pca_deter)+geom_point(aes(x=PC1,y=PC2,col=Group_name))
ggplot(res_pca_deter)+geom_point(aes(x=PC1,y=PC2,col=Group_Sex))
ggplot(res_pca_deter)+geom_point(aes(x=PC1,y=PC3,col=Group_Sex))
ggplot(res_pca_deter)+geom_point(aes(x=PC1,y=PC3,col=Group_Sex))

ggplot(res_pca_deter)+geom_boxplot(aes(x=Group_Sex,y=PC1,fill=Group_Sex))
