
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

res_gsea_all<-fread("../methyl/outputs/tf_signature_enrichment_in_DMC/res_gsea_genescore_ranked.csv.gz")

res_gsea_all[,top10:=padj<=sort(padj)[10],by="compa"]

res_gsea_all[top10==T,.(pathway,compa,padj,NES)]
res_gsea_all[pathway=="STAT3",.(pathway,compa,padj,NES)]

#C vs L
#top 20 regulons by gsea in C vs L
res_gsea<-res_gsea_all[compa=="C.L"][,regulon:=pathway][,-"pathway"]
regulons_ord<-res_gsea[padj<=sort(padj)[15]][order(padj)]$regulon

ggplot(res_gsea[padj<=sort(padj)[15]])+geom_col(aes(x=regulon,y=-log10(padj),fill=NES))+scale_x_discrete(limits=regulons_ord)

regulons_ord<-res_gsea[!str_detect(regulon,'e$')][padj<=sort(padj)[15]][order(padj)]$regulon
ggplot(res_gsea[!str_detect(regulon,'e$')][padj<=sort(padj)[15]])+
  geom_col(aes(x=regulon,y=-log10(padj),fill=NES))+
  scale_x_discrete(limits=regulons_ord)

#=> TF epigen alteration by lineage
res_gsea<-res_gsea_all[compa=="C.L"][,regulon:=pathway][,-"pathway"]

regulons_lineage<-fread("../singlecell/outputs/06-SCENIC/cbps0-8_clean/tf_markers_lineage.csv.gz")

regulons_lineage_f<-regulons_lineage[p_val_adj<10^-30&avg_log2FC>0]
table(regulons_lineage_f$lineage)
res_gsea_lineage<-merge(res_gsea,regulons_lineage_f[,.(regulon,lineage,p_val_adj)])

table(res_gsea_lineage[padj<0.001]$lineage)
res_gsea_lineage[padj<0.001&!str_detect(regulon,'e$')]

ggplot(res_gsea_lineage[padj<0.001&!str_detect(regulon,'e$')])+geom_bar(aes(x=lineage))
ggplot(res_gsea_lineage[padj<0.001&!str_detect(regulon,'e$')])+geom_boxplot(aes(x=lineage,y=-log10(padj)))



#CF vs LF
#top 20 regulons by gsea in C vs L
res_gsea<-res_gsea_all[compa=="CF.LF"][,regulon:=pathway][,-"pathway"]
regulons_ord<-res_gsea[padj<=sort(padj)[15]][order(padj)]$regulon

ggplot(res_gsea[padj<=sort(padj)[15]])+geom_col(aes(x=regulon,y=-log10(padj),fill=NES))+scale_x_discrete(limits=regulons_ord)

regulons_ord<-res_gsea[!str_detect(regulon,'e$')][padj<=sort(padj)[15]][order(padj)]$regulon
ggplot(res_gsea[!str_detect(regulon,'e$')][padj<=sort(padj)[15]])+
  geom_col(aes(x=regulon,y=-log10(padj),fill=NES))+
  scale_x_discrete(limits=regulons_ord)

#=> TF epigen alteration by lineage
res_gsea<-res_gsea_all[compa=="CF.LF"][,regulon:=pathway][,-"pathway"]

res_gsea_lineage<-merge(res_gsea,regulons_lineage_f[,.(regulon,lineage,p_val_adj)])
res_gsea_lineage<-res_gsea_lineage[padj<0.001&p_val_adj<0.001]
ggplot(res_gsea_lineage[!str_detect(regulon,'e$')&lineage%in%c("HSC","MPP","Erythroid","Lymphoid","Myeloid")])+
  geom_col(aes(x=regulon,y=-log10(padj),fill=-log10(p_val_adj)))+facet_grid( .~ lineage, space = "free_x",scales="free_x")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,size = 9),
        axis.title=element_text(size=12,face="bold"))

table(res_gsea_lineage[padj<0.001]$lineage)
res_gsea_lineage[padj<0.001&!str_detect(regulon,'e$')]

ggplot(res_gsea_lineage[padj<0.001&!str_detect(regulon,'e$')])+geom_bar(aes(x=lineage))
ggplot(res_gsea_lineage[padj<0.001&!str_detect(regulon,'e$')])+geom_boxplot(aes(x=lineage,y=-log10(padj)))




#over representation test [à finir]
res_meth[,gene_score_scaled:=scale(GeneScore,center = F),by="compa"]
unique(res_meth[gene=="SOCS3",.(gene,compa,GeneScore,gene_score_scaled)])
unique(res_meth[gene=="HES1",.(gene,compa,GeneScore,gene_score_scaled)])


