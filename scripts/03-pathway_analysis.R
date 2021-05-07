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
nrow(as.data.frame(res_gsea_go))#685

dotplot(res_gsea_go,showCategory=20)
gsea_go<-data.table(as.data.frame(res_gsea_go))
gsea_go[,gene_score.avg:=mean(resg$gene_score_add[resg$gene %in% tr(core_enrichment,tradEntrezInSymbol = T)],na.rm=T),.(ID)]

saveRDS(res_gsea_go,fp(out,"res_gsea_go.rds"))
fwrite(gsea_go[order(p.adjust)],fp(out,"res_gsea_go.csv"))

dotplot(res_gsea_go,x=gsea_go[order(p.adjust)]$gene_score.avg[1:40],showCategory=40)

emapplot(pairwise_termsim(res_gsea_go,showCategory = 40),showCategory = 40)
