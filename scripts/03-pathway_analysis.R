#Gene score calculation and validation
source("scripts/utils/new_utils.R")
source("scripts/utils/methyl_utils.R")
library(clusterProfiler)
library(enrichplot)
library(org.Hs.eg.db)

out<-"outputs/03-pathway_analysis"
dir.create(out)


resg<-fread("outputs/02-gene_score_calculation_and_validation/res_genes.csv.gz")

#OR KEGG
plot(density(resg$gene_score_add))
abline(v=300)
resg[gene_score_add>300]
res_or_kegg<-enrichKEGG(bitr(resg[gene_score_add>300]$gene,fromType = "SYMBOL",toType = "ENTREZID",OrgDb = org.Hs.eg.db)$ENTREZID,
                 organism = "hsa",pvalueCutoff = 0.05)
as.data.frame(res_or_kegg)


emapplot(pairwise_termsim(res_or_kegg))

#GSE KEGG
gene_scores<-resg$gene_score_add
names(gene_scores)<-resg$gene
gene_scores<-sort(gene_scores,decreasing = T)
head(gene_scores)

genes.df<-bitr(names(gene_scores),
               fromType = 'SYMBOL',
               toType = 'ENTREZID',
               OrgDb = org.Hs.eg.db)
head(genes.df)
head(gene_scores)
gene_scores<-gene_scores[genes.df$SYMBOL]
names(gene_scores)<-genes.df$ENTREZID
res_gsea_kegg<- gseKEGG(geneList     = rank(gene_scores),
                        organism     = 'hsa', 
                        minGSSize    = 50,
                        pvalueCutoff = 0.001,
                        verbose = FALSE)
nrow(as.data.frame(res_gsea_kegg))#63

dotplot(res_gsea_kegg,showCategory=63)
gsea_kegg<-data.table(as.data.frame(res_gsea_kegg))
gsea_kegg[,gene_score.avg:=mean(resg$gene_score_add[resg$gene %in% tr(core_enrichment,tradEntrezInSymbol = T)],na.rm=T),.(ID)]

dotplot(res_gsea_kegg,x=gsea_kegg$gene_score.avg,showCategory=63)

emapplot(pairwise_termsim(res_gsea_kegg,showCategory = 63))
saveRDS(res_gsea_kegg,fp(out,"res_gsea_kegg.rds"))
fwrite(gsea_kegg,fp(out,"res_gsea_kegg.csv"))

#GSE GO

res_gsea_go<- gseGO(geneList     = rank(gene_scores), 
                        minGSSize    = 50,
                        pvalueCutoff = 0.001,
                        eps = 0,
                        OrgDb = org.Hs.eg.db)
nrow(as.data.frame(res_gsea_go))#699

dotplot(res_gsea_go,showCategory=20)
gsea_go<-data.table(as.data.frame(res_gsea_go))
gsea_go[,gene_score.avg:=mean(resg$gene_score_add[resg$gene %in% tr(core_enrichment,tradEntrezInSymbol = T)],na.rm=T),.(ID)]

saveRDS(res_gsea_go,fp(out,"res_gsea_go.rds"))
fwrite(gsea_go[order(p.adjust)],fp(out,"res_gsea_go.csv"))

dotplot(res_gsea_go,x=gsea_go[order(p.adjust)]$gene_score.avg[1:40],showCategory=40)

emapplot(pairwise_termsim(res_gsea_go,showCategory = 40),showCategory = 40)

#GSE GWAS

source("scripts/utils/new_utils.R")
library(clusterProfiler)
library(enrichplot)
out<-"outputs/03-pathway_analysis"

resg<-fread("outputs/02-gene_score_calculation_and_validation/res_genes.csv.gz")
gwas_genes<-fread("../methyl/ref/gwas/list_ref_GWAS_042721.csv",
                  select = c(2,3,4),
                  col.names = c("reported_gene","mapped_genes","disease_trait"),
                  skip = 1)

#☺1) with reported genes
gwas_genes[reported_gene=="NR",reported_gene:=NA]
gwas_genes_rep<-unique(gwas_genes[!is.na(reported_gene)],by=c("reported_gene","disease_trait"))
gwas_genes_rep #92k gene-trait associations
length(unique(gwas_genes_rep$disease_trait)) #4k5 disease traits
gwas_genes_rep[,n.gene.trait:=.N,by="disease_trait"]
summary(unique(gwas_genes_rep[,.(disease_trait,n.gene.trait)])$n.gene.trait)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 1.00    2.00    5.00   20.38   11.00 1532.00 

# remove traits if < 10 genes 
gwas_genes_rep10<-gwas_genes_rep[n.gene.trait>=10]
gwas_genes_rep10 #81k gene-trait associations
fwrite(gwas_genes_rep10,"../methyl/ref/gwas/reported_gene_traits_GWAS_10genes.csv")
gwas_genes_rep10<-fread("../methyl/ref/gwas/reported_gene_traits_GWAS_10genes.csv")

length(unique(gwas_genes_rep10$disease_trait)) #1299 disease traits
gwas_genes_rep10[,n.gene.trait:=.N,by="disease_trait"]
summary(unique(gwas_genes_rep10[,.(disease_trait,n.gene.trait)])$n.gene.trait)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 10.00   14.00   23.00   62.75   50.00 1532.00 

gene_scores<-resg$gene_score_add
names(gene_scores)<-resg$gene
gene_scores<-sort(gene_scores,decreasing = T)

res_gsea_gwas<-GSEA(geneList = rank(gene_scores),
                    TERM2GENE = gwas_genes_rep10[,.(disease_trait,reported_gene)],
                    maxGSSize    = 500,
                    eps = 0,
                    pvalueCutoff = 0.01
                    )
nrow(as.data.frame(res_gsea_gwas))#86/1300
dotplot(res_gsea_gwas,showCategory=20)
emapplot(pairwise_termsim(res_gsea_gwas,showCategory = 86),showCategory = 86)

gsea_gwas<-data.table(as.data.frame(res_gsea_gwas))
gsea_gwas[,gene_score.avg:=mean(resg$gene_score_add[resg$gene %in% tr(core_enrichment,tradEntrezInSymbol = F)],na.rm=T),.(ID)]

saveRDS(res_gsea_gwas,fp(out,"res_gsea_gwas.rds"))
fwrite(gsea_gwas[order(p.adjust)],fp(out,"res_gsea_gwas.csv"))


#2) if not satisfying, with mapped gene #satisfying

#PERMUT Patway====
library(limma)
library(parallel)
source("scripts/utils/new_utils.R")
set.seed(1234)
methf<-fread("datasets/cd34/meth_data_filtered.csv.gz")

mtd_f<-fread("datasets/cd34/metadata_cl_pcs_040521.csv")
cpgs_score<-fread("outputs/02-gene_score_calculation_and_validation/cpgs_genes_annot_and_weight.csv.gz")

res_perm<-Reduce(rbind,mclapply( 1:1000,function(i){

  mtd_f$group_sex<-sample(mtd_f$group_sex)
  
  #limma
  design<-model.matrix(~0 + group_sex   + batch+ group_complexity_fac +mat.age  + latino + PC2,
                       data = data.frame(mtd_f,row.names = "sample"))
  
  fit <- lmFit(data.frame(methf,row.names = "cpg_id")[,mtd_f$sample], design)
  
  cont.matrix <- makeContrasts(C.L = "(group_sexCTRL_F+group_sexCTRL_M)-(group_sexLGA_F+group_sexLGA_M)",
                               levels=design)
  fit2  <- contrasts.fit(fit, cont.matrix)
  fit2  <- eBayes(fit2)
  
  res<-data.table(topTable(fit2,coef = "C.L",n = Inf),keep.rownames = "cpg_id")
  res[,cpg_id:=as.numeric(cpg_id)]
  cols1<-c("cpg_id","P.Value","adj.P.Val","AveExpr","logFC")
  cols2<-c("cpg_id","pval","padj","avg.meth","meth.change")
  res<-res[,(cols2):=.SD,.SDcols=cols1][,.SD,.SDcols=cols2]
  
  res_anno<-merge(res,cpgs_score,by="cpg_id")

  #genescore
  res_anno[,cpg_score:=-log10(pval)*meth.change*links_weight*regul_weight] #divided by n_sample ro normalized gene score ~ n_sampleres_anno[,n.cpg_weight:=(1/sum(1/(abs(cpg_score)+1)))^(1/5),by="gene"] #n.cpg_weight to reduce the influence of the n_cpg by gene to the GeneScore
  
  res_anno[,region_type:=ifelse(abs(tss_dist)<=2000,"promoter","other")]
  res_anno[is.na(tss_dist),region_type:="other"]
  res_anno[,n_cpg_weight_region:=(1/sum(1/(abs(cpg_score)+1)))^(1/3.8),by=c('region_type',"gene")]
  
  res_anno[region_type=="promoter",gene_score_region:=sum(cpg_score)*n_cpg_weight_region,by=c("gene")]
  res_anno[region_type=="other",gene_score_region:=sum(abs(cpg_score))*n_cpg_weight_region,by=c("gene")]
  res_anno[,gene_score_add:=sum(unique(abs(gene_score_region)),na.rm = T),by="gene"]
  
  return(unique(res_anno,by="gene")[,perm:=i][order(gene)][,.(gene,gene_score_add,perm)])
  
  
},mc.cores = 20))

fwrite(res_perm,"outputs/02-gene_score_calculation_and_validation/res_1000perm_genescore_add.csv.gz")

 #pathway
library(clusterProfiler)
library(org.Hs.eg.db)

#kegg
res_k_perm<-Reduce(rbind,mclapply(1:1000,function(i){
  resg<-res_perm[perm==i]
  gene_scores<-resg$gene_score_add
  names(gene_scores)<-resg$gene
  gene_scores<-sort(gene_scores,decreasing = T)
  
  genes.df<-bitr(names(gene_scores),
                 fromType = 'SYMBOL',
                 toType = 'ENTREZID',
                 OrgDb = org.Hs.eg.db)
  gene_scores<-gene_scores[genes.df$SYMBOL]
  names(gene_scores)<-genes.df$ENTREZID
  res_gsea_kegg<- data.table(as.data.frame(gseKEGG(geneList     = rank(gene_scores),
                        organism     = 'hsa', 
                        minGSSize    = 50,
                        pvalueCutoff = 1,
                        verbose = FALSE)))
  
  return(res_gsea_kegg[,perm:=i][,.(ID,p.adjust,perm)])}))

fwrite(res_k_perm,"outputs/03-pathway_analysis/res_1000perm_kegg.csv.gz")

#go bp
res_g_perm<-Reduce(rbind,mclapply(1:1000,function(i){
  resg<-res_perm[perm==i]
  gene_scores<-resg$gene_score_add
  names(gene_scores)<-resg$gene
  gene_scores<-sort(gene_scores,decreasing = T)
  
  genes.df<-bitr(names(gene_scores),
                 fromType = 'SYMBOL',
                 toType = 'ENTREZID',
                 OrgDb = org.Hs.eg.db)
  gene_scores<-gene_scores[genes.df$SYMBOL]
  names(gene_scores)<-genes.df$ENTREZID

  res_gsea_go<- data.table(as.data.frame(gseGO(geneList     = rank(gene_scores), 
                                               ont = "BP",
                                                minGSSize    = 50,
                                                pvalueCutoff = 1,
                                                eps = 0,
                                                OrgDb = org.Hs.eg.db)))
  
  return(res_gsea_go[,perm:=i][,.(ID,p.adjust,perm)])}))

fwrite(res_g_perm,"outputs/03-pathway_analysis/res_1000perm_go.csv.gz")

#gwas
gwas_genes_rep10<-fread("ref/gwas/reported_gene_traits_GWAS_10genes.csv")

res_g_perm<-Reduce(rbind,mclapply(1:1000,function(i){
  resg<-res_perm[perm==i]
  gene_scores<-resg$gene_score_add
  names(gene_scores)<-resg$gene
  gene_scores<-sort(gene_scores,decreasing = T)
  
  genes.df<-bitr(names(gene_scores),
                 fromType = 'SYMBOL',
                 toType = 'ENTREZID',
                 OrgDb = org.Hs.eg.db)
  gene_scores<-gene_scores[genes.df$SYMBOL]
  names(gene_scores)<-genes.df$ENTREZID

  res_gsea_gwas<- data.table(as.data.frame(GSEA(geneList = rank(gene_scores),
                                              TERM2GENE = gwas_genes_rep10[,.(disease_trait,reported_gene)],
                                              maxGSSize    = 500,
                                              eps = 0,
                                              pvalueCutoff = 1
                                              )))
  
  return(res_gsea_gwas[,perm:=i][,.(ID,p.adjust,perm)])}))

fwrite(res_g_perm,"outputs/03-pathway_analysis/res_1000perm_gwas.csv.gz")
