
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
