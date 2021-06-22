source("scripts/utils/new_utils.R")


res_b<-fread("outputs/07-LGA_vs_Ctrl_Basal/res_pseudobulkDESeq2_all_cbps.csv.gz")[,stimulation:=F]

res_h<-fread("outputs/09-LGA_vs_Ctrl_Activated/res_pseudobulkDESeq2_all_cbps.csv.gz")[,stimulation:=T]

res<-rbind(res_b_hsc,res_h_hsc)

ggplot(res)+
  geom_point(aes(x=log2FoldChange,y=-log10(padj),col=padj<0.05&abs(log2FoldChange)>0.6))+
  facet_wrap("stimulation")+scale_color_manual(values = c("grey","red"))+
  theme_minimal()
                           

ggplot(res_hsc)+
  geom_point(aes(x=log2FoldChange,y=-log10(padj),col=padj<0.05&abs(log2FoldChange)>0.6))+
  facet_wrap("stimulation")+scale_color_manual(values = c("grey","red"))+
  theme_minimal()+theme(legend.position = 'bottom')
                           


res_b_hsc<-fread("outputs/07-LGA_vs_Ctrl_Basal/res_pseudobulkDESeq2_by_lineage.csv.gz")[lineage=="HSC"][,stimulation:=F]

res_h_hsc<-fread("outputs/09-LGA_vs_Ctrl_Activated/res_pseudobulkDESeq2_by_lineage.csv.gz")[lineage=="HSC"][,stimulation:=T]

res_hsc<-rbind(res_b_hsc,res_h_hsc)

ggplot(res_hsc)+
  geom_point(aes(x=log2FoldChange,y=-log10(padj),col=padj<0.05&abs(log2FoldChange)>0.6))+
  facet_wrap("stimulation")+scale_color_manual(values = c("grey","red"))+
  theme_minimal()
                           

ggplot(res_hsc)+
  geom_point(aes(x=log2FoldChange,y=-log10(padj),col=padj<0.05&abs(log2FoldChange)>0.6))+
  facet_wrap("stimulation")+scale_color_manual(values = c("grey","red"))+
  theme_minimal()+theme(legend.position = 'bottom')
                           
