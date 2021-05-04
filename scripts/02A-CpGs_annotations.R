#ANNOT CPGs

source("scripts/utils/new_utils.R")
cpgs<-fread("datasets/cd34/CD34_angle_119_noEmptyLocis_withConfScore_withoutChrXY.txt",
            select = c("id","chr","start"),
            col.names = c("cpg_id","chr","pos"))
out<-"outputs/02A-CpGs_annotations"
dir.create(out)

#1)cpgs_genes_links
#by distance to tss
#first gene within +/-200kb

genes_in_200kb<-bed_inter(cpgs[,start:=pos-200e3][start<0,start:=0][,end:=pos+200e3][,.(chr,start,end,cpg_id)][order(chr,start)],
          "ref/hg19/tss_genes.bed",
          select = c(5,6,8,9,4),col.names = c("chr","tss_pos","gene","strand","cpg_id"))

genes_in_200kb
cpgs_genes_tss<-merge(cpgs,genes_in_200kb,all.x=T)[,.(chr,pos,cpg_id,gene,tss_pos,strand)]
cpgs_genes_tss[strand=="+",tss_dist:=pos-tss_pos][strand=="-",tss_dist:=tss_pos-pos]
cpgs_genes_tss[,n.gene.cpg:=length(unique(gene)),by="cpg_id"]
summary(cpgs_genes_tss$n.gene.cpg)
  # Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
  #  1.00    6.00   10.00   12.36   17.00  101.00 
summary(abs(cpgs_genes_tss$tss_dist))
   # Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA
   #    0   45620   96440   96946  147918  200000  131093 

#why a lot of NA ?
cpgs_genes_tss[is.na(tss_dist)] #131093 CpGs with no gene at +/-200kb

cpgs_genes_tss[!is.na(tss_dist),closest.rank:=rank(abs(tss_dist)),by="cpg_id"]

summary(abs(cpgs_genes_tss[closest.rank==1]$tss_dist))
  # Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
  #     0    4455   18193   35849   50208  200000  131093 

summary(abs(cpgs_genes_tss[closest.rank==2]$tss_dist))
   # Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
   #    0   10139   31030   48263   72131  199996 

fwrite(cpgs_genes_tss,fp(out,"cpg_tss_200kb.csv.gz"))


cpgs_genes_1<-cpgs_genes_tss[!is.na(tss_dist)&closest.rank==1]

#for >1 cpg_gene link (because multiple TSS), keep most closest cpg_gene link
cpgs_genes_1<-unique(cpgs_genes_1[order(cpg_id,gene,abs(tss_dist))],by=c("cpg_id","gene")) 

cpgs_genes_1
rm(cpgs_genes_tss)

fwrite(cpgs_genes_1,fp(out,"cpgs_closest_gene_tss_linked_within_200kb_around.tsv"),sep="\t")

#by presence in eQTL region
#need first create eQTR
#see 02A1-create_eQTR
cpgs_eQTR<-bed_inter(a=cpgs[,start:=pos][,end:=pos+1][,.(chr,start,end,cpg_id)][order(chr,start)],
          b=eqtr[,.(chr,start,end,eqtr_id)],
          select = c(4,8),col.names = c("cpgID","gene"))

cpgs_eQTR #887k  cpgs - snp+/-500pb match
unique(cpgs_eQTR) #554k cpg-gene match
unique(cpgs_eQTR,by="cpgID")#246k cpgs linked to a gene
unique(cpgs_eQTR,by="gene")#13k genes linked to a cpg
cpgs_eQTR[,in_eQTR:=T]

#merge the 2 links
cpgs_genes_4[,in_eQTR:=F]
cpgs_genes<-merge(cpgs_genes_4[,.(cpgID,gene,in_eQTR,tss_dist)],cpgs_eQTR,all=T)
cpgs_genes[in_eQTR==T]

#merge with coord
cpgs_genes<-merge(cpgs[,.(cpgID,chr,pos)],cpgs_genes,all.x = T,by="cpgID")

fwrite(cpgs_genes,fp(out,"cpgs_4closest_gene_tss_linked_within_200kb_around_and_eQTR_linked_genes.tsv"),sep="\t")

#2)cpgs_regulatory region
#ensembl regulatory reg matching
cpgs_ensembl<-bed_inter(a=unique(cpgs,by="cpgID")[,chr:=str_remove(chr,'chr')][,start:=pos][,end:=pos+1][,.(chr,start,end,cpgID)][order(chr,start)],
          b="ref/ensembl_regulatory/ensembl_regulatory_hg19.bed",
          select = c(4,8),col.names = c("cpgID","ensembl_regulatory_domain"))



cpgs_ensembl[,ensembl_regulatory_domain:=paste(ensembl_regulatory_domain,collapse = "/"),by="cpgID"] 
cpgs_ensembl<-unique(cpgs_ensembl)
cpgs_ensembl

cpgs_chromatin<-res[,.(cpgID,chr,pos,chromatin_state,chromatin_feature)]

cpgs_reg<-merge(cpgs_chromatin,cpgs_ensembl,all.x=T,by="cpgID")

fwrite(cpgs_reg,fp(out,"cpgs_annot_ensembl_regulatory_domain_and_chromHMM_chromatin_features.tsv"),sep="\t")

cpgs_anno<-merge(cpgs_reg,cpgs_genes,all=T,by="cpgID")

fwrite(unique(cpgs_anno[order(cpgID,gene)][,.(cpgID,chr,pos,in_eQTR,tss_dist,gene,chromatin_state,chromatin_feature,ensembl_regulatory_domain)]),
       fp(out,"cpgs_annot.tsv"),sep="\t")


# Calculate GeneScore

#1) calculate cpg weight
cpgs_genes<-cpgs_anno[!is.na(gene)&gene!=""]
unique(cpgs_genes,by="cpgID") #741408/~800k cpgs link to a gene

unique(cpgs_genes[in_eQTR==T],by="cpgID") #dont 246k linked thx to eQTR
cpgs_genes[,double.linked:=any(in_eQTR==T)&any(in_eQTR==F),by=c("cpgID","gene")]
unique(cpgs_genes[double.linked==T],by="cpgID")#dont 63504 also gene linked based on tss
#=> + ~180k cpgs gene link thx to eQTR

# a) linksWeight
cpgs_genes[in_eQTR==F,links_score:=sapply(abs(tss_dist),function(x){
  if(x<1000)return(1)
  else if(x<20000)return(0.5+0.5*sqrt(1000/x))
  else return(0.5*sqrt(20000/x))
    })]

cpgs_genes[in_eQTR==T,links_score:=1]

cpgs_genes[,links_weight:=max(links_score),by=c("cpgID","gene")]


#b) Regulatory Weights
cpgs_reg<-unique(cpgs_genes,by=c("cpgID"))
unique(cpgs_reg$chromatin_feature)
cpgs_reg[,chromatin_score:=sapply(chromatin_feature,function(x){
      if(x%in%c("Promoter","Active Enhancer"))return(1)
      else if(x=="Inactive Enhancer")return(0.75)
      else if(x%in%"Gene Body")return(0.5)
      else return(0)
    })]


unique(cpgs_reg$chromatin_feature)
cpgs_reg[,ensembl_reg_score:=sapply(ensembl_regulatory_domain,function(x){
    vecX<-strsplit(as.character(x),"/")[[1]]
    if(any(c("CTCF_binding_site","promoter","enhancer")%in%vecX)){
      score<-0.5
    }else if(any(c("open_chromatin_region","promoter_flanking_region")%in%vecX)){
      score<-0.25
    }else{
      score<-0
    }
    if("TF_binding_site"%in%vecX){
      score<-score+0.5
    }
    return(score)
  })]

cpgs_reg[,regul_weight:=(0.5+1.5*((chromatin_score+ensembl_reg_score)/2))]

cpgs_score<-merge(cpgs_genes,cpgs_reg[,.(cpgID,chromatin_score,ensembl_reg_score,regul_weight)],by="cpgID")
fwrite(cpgs_score,fp(out,"cpgs_genes_annot_and_weight_reference.tsv"),sep="\t")

