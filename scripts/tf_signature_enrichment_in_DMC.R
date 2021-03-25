
#TF signature enriched in DMC  ?
source("../methyl/scripts/utils/new_utils.R")
out<-here("outputs/tf_signature_enrichment_in_DMC")
dir.create(out)
regulons_list<-readRDS("../singlecell/outputs/06-SCENIC/cbps0-8_clean/regulons_list.rds")

#in GeneScore
res_meth<-fread("../methyl/outputs/model14_without_iugr/2020-10-08_all_res_with_perm.csv")

res_meth[compa=="CF.LF"&gene=="HES1"]

#GSEA
library(fgsea)
?fgsea
res_methg<-unique(res_meth[order(compa,gene,pval)],by=c("compa","gene"))
res_methg[GeneScore>0,gs_rank:=rank(GeneScore),by="compa"]
res_methg[GeneScore<=0,gs_rank:=-rank(abs(GeneScore)),by="compa"]


genes_rank<-res_methg[compa=="CF.LF"]$gs_rank
names(genes_rank)<-res_methg[compa=="CF.LF"]$gene
res_gsea<-fgsea(pathways=regulons_list,
      stats=genes_rank,eps=0)

fwrite(res_gsea,f)
res_gsea[padj<0.001]
res_gsea[order(padj)]$pathway[1:10]
res_gsea[pathway=="STAT3e"]
res_gsea[pathway=="STAT3"]

fwrite(res_gsea,fp(out,"res_cflf_gsea_genescore_ranked.csv.gz"),sep=";")

#for each compa
#run tf1-

res_gsea_all<-fread("outputs/tf_signature_enrichment_in_DMC/res_gsea_genescore_ranked.csv.gz")

res_gsea_all[,top10:=padj<=sort(padj)[10],by="compa"]

res_gsea_all[top10==T,.(pathway,compa,padj,NES)]
res_gsea_all[pathway=="STAT3",.(pathway,compa,padj,NES)]

#=> TF epigen alteration by lineage[a finir]

regulons_lineage<-


#over representation test [à finir]
res_meth[,gene_score_scaled:=scale(GeneScore,center = F),by="compa"]
unique(res_meth[gene=="SOCS3",.(gene,compa,GeneScore,gene_score_scaled)])
unique(res_meth[gene=="HES1",.(gene,compa,GeneScore,gene_score_scaled)])

