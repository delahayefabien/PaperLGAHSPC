#DMR 
library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library("org.Hs.eg.db")
source("scripts/utils/new_utils.R")

out<-"outputs/02-gene_score_calculation_and_validation/dmr/"
dir.create(out)
fwrite(res[,start:=pos][,end:=pos+1][order(chr,start)][!is.na(start)][,.(chr,start,end,pval)],fp(out,"res_cpgs.bed"),sep="\t")

system("tools/comb-p/combined-pvalues/cpv/pipeline.py -p outputs/02-gene_score_calculation_and_validation/ctrl_lga --region-filter-p 0.1 --anno hg19 --seed 0.01 --dist 1000 outputs/02-gene_score_calculation_and_validation/res_cpgs.bed > outputs/02-gene_score_calculation_and_validation/ctrl_vs_lga_pvals.dmr.bed")
#
res_dmr<-fread("outputs/02-gene_score_calculation_and_validation/ctrl_lga.regions-p.bed.gz")
res_dmr<-res_dmr[,chr:=`#chrom`][,-1][,.(chr,start,end,min_p,n_probes,z_p,z_sidak_p)]
#merge with gene

dmr<-makeGRangesFromDataFrame(res_dmr[,1:3])
anno_dmr <- annotatePeak(dmr, tssRegion=c(-2000, 2000),
                       TxDb=TxDb.Hsapiens.UCSC.hg19.knownGene,
                       annoDb="org.Hs.eg.db")
anno_dmr<-as.data.frame(anno_dmr)
head(anno_dmr,10)
fwrite(anno_dmr,fp(out,"chipseaker_anno_dmr.csv.gz"),sep=";")
anno_dmr<-fread(fp(out,"chipseaker_anno_dmr.csv.gz"),sep=";")

dmr_tss<-data.table(anno_dmr)[,chr:=seqnames][,tss_dist:=distanceToTSS][,gene:=SYMBOL][,.(chr,start,end,tss_dist,gene)]

res_dmr<-merge(res_dmr,dmr_tss,all.x=T)
res_dmr[z_sidak_p<0.05&abs(tss_dist)<2000]$gene
res_dmr[gene=="SOCS3"]
fwrite(res_dmr,fp(out,"res_dmr_anno.csv.gz"),sep=";")


