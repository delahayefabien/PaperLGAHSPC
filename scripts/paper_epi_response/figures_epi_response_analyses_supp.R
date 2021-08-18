
#2B-HTO_signature
#scDEGs
res_dup<-fread("../singlecell/outputs/03-HTOs_stimulation/duplicates/scdegs_hto_vs_not_wilcoxon.csv")
#volcano
genes_of_interest<-c("SOCS3","HES1","JUN","FOS","JUNB","ZFP36","EGR1",
                      "DUSP2","DUSP1","FOSB","SOCS1","KLF2","KLF4",
                       "PLK2","PLK3","ID1","MYC","","ID2","IDS","RGCC")

ggplot(res_dup,aes(x=avg_log2FC,y=-log10(p_val_adj),col=p_val_adj<0.01&abs(avg_log2FC)>0.4))+
  geom_point()+ 
  geom_label_repel(aes(label = ifelse(p_val_adj<0.01&
                                        abs(avg_log2FC)>0.4&
                                        gene%in%genes_of_interest,gene,"")),
                   max.overlaps = 3000,
                   box.padding   = 0.35,
                   point.padding = 0.5,
                   segment.color = 'grey50')+
  scale_color_manual(values = c("grey","red")) +
  theme_minimal() +
  theme(legend.position = "bottom")+facet_wrap("sample")
ggsave(fp(out,"2B-volcano_scDEGs_antibody_activation_signature_in_duplicates.pdf"))

table(res_dup[p_val_adj<0.01&abs(avg_log2FC)>0.4]$sample)
# ctrlM518 ctrlM537 ctrlM555 
#      572      675      331 

res_dup[,degs_in_3:=sum(p_val_adj<0.01&abs(avg_log2FC)>0.4)==3,by="gene"]
length(unique(res_dup[degs_in_3==T]$gene)) #210
cat(paste(unique(res_dup[degs_in_3==T]$gene),collapse = "\n"))

#supp : signature enrichment by lineage 
source("scripts/utils/new_utils.R")
res_sign_kegg<-fread("outputs/08-HTO_signature/by_lineage/res_hto_signature_kegg_by_lineage.csv")

res_sign_kegg[,n.enriched:=Count]
res_sign_kegg[,n_gene_set:=as.numeric(str_extract(BgRatio,"^[0-9]+"))]
res_sign_kegg[,pct.enriched:=n.enriched/n_gene_set]
ggplot(res_sign_kegg[lineage%in%c("HSC","MPP/LMPP","Erythro-Mas","Lymphoid","Myeloid")&p.adjust<0.05])+
  geom_point(aes(x=Description,y=pct.enriched,size=n.enriched,col=-log10(p.adjust)))+
  facet_grid(~lineage,scales = "free_x",space = "free")+
  theme(axis.text.x = element_text(angle = 60,vjust =1,hjust = 1))
ggsave("outputs/figures_epi_response/figure2/2Bsupp-kegg_dotplot_antibody_activation_signature_by_lineage.pdf")



res_sign_go_list<-readRDS("outputs/08-HTO_signature/by_lineage/res_hto_signature_go_bp_by_lineage.rds")
res_sign_go<-Reduce(function(x,y)rbind(x,y,fill=TRUE),lapply(names(res_sign_go_list), function(lin)data.table(as.data.frame(res_sign_go_list[[lin]]))[,lineage:=lin]))

res_sign_go[,n.enriched:=Count]
res_sign_go[,n_gene_set:=as.numeric(str_extract(BgRatio,"^[0-9]+"))]
res_sign_go[,pct.enriched:=n.enriched/n_gene_set]
res_sign_go[,top10:=p.adjust<=max(sort(p.adjust)[1:10],na.rm = T),by="lineage"]
ggplot(res_sign_go[top10==T&lineage%in%c("HSC","MPP/LMPP","Erythro-Mas","Lymphoid","Myeloid")&p.adjust<0.05])+
  geom_point(aes(x=Description,y=pct.enriched,size=n.enriched,col=-log10(p.adjust)))+
  facet_grid(~lineage,scales = "free_x",space = "free")+
  theme(axis.text.x = element_text(angle = 66,vjust =1,hjust = 1))
ggsave("outputs/figures_epi_response/figure2/2Bsupp-go_dotplot_antibody_activation_signature_by_lineage.pdf")


