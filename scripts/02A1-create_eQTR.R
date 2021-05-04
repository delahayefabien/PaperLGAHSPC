#02A1-create_eQTR
#2020-05-30
#config et library
options(stringsAsFactors=F)
set.seed(12345)
#library(limma)
library(data.table)
out<-"outputs/02A1-create_eQTR"
dir.create(out)
#I) with whole_blood eQTL
eqtls_wb<-fread("ref/eQTL/GTEx_Analysis_v8_eQTL/Whole_Blood.v8.signif_variant_gene_pairs.txt.gz")
eqtls_wb #2 414 653 signif_variant_gene_pairs
unique(eqtls_wb$gene_id) #1000 genes
eqtls_wb[,n_eqtl.gene:=.N,by="gene_id"]

summary(unique(eqtls_wb,by="gene_id")$n_eqtl.gene)
   # Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
   #  1.0    22.0    87.0   195.4   230.0  8571.0 

#extract chr and pos hg38 of snps
eqtls_wb[,chr.hg38:=str_extract(variant_id,"^chr[0-9XYM]{1,2}") ]
eqtls_wb<-eqtls_wb[!is.na(chr.hg38),]
eqtls_wb[,pos.hg38:=as.numeric(str_sub(str_extract(variant_id,"_[0-9]+"),2)) ]
eqtls_wb

#trans in hg19 :
translator<-fread("ref/eQTL/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze.lookup_table.txt.gz",
                  select = c(1,8))
translator


eqtls_wb<-merge(eqtls_wb,translator,all.x=T,by="variant_id")


eqtls_wb[,chr.hg19:=paste0("chr",str_extract(variant_id_b37,"^[0-9XY]{1,2}"))]

eqtls_wb[,pos.hg19:=as.numeric(str_sub(str_extract(variant_id_b37,"_[0-9]+"),2))]
eqtls_wb

eqtls_wb[,chr:=chr.hg19]
eqtls_wb[,pos:=pos.hg19]

#calculate dist snps asso to a same gene
GetDistToNext<-function(positions){
  return(c(sapply(1:(length(positions)-1), function(i)abs(positions[i]-positions[i+1])),NA))
}

eqtls_wb[,dist_to_next:=GetDistToNext(pos),by="gene_id"]
summary(eqtls_wb$dist_to_next)
# Min.   1st Qu.    Median      Mean   3rd Qu.      Max.      NA's 
#    0       114       397     12418      1158 180679628     39827 

#definition of eQTR : 1 to 10kb region around most signif eQTL locally 

# 1) find min.local, i.e eqtl with pval in top25% of eQTL of the gene, and first eQTL at +/-5kb

eqtls_wb[,minLocal.cand:=pval_nominal<quantile(pval_nominal,0.25),by=c("gene_id")]

is.minLocal<-function(dists,pvals){
  isMins<-sapply(1:length(dists), function(i){
    locisA5kb<-which(dists>(dists[i]-5000)&dists<(dists[i]+5000))
    
    if(all(pvals[locisA5kb]>=pvals[i])){
      return(T)
    }else{
      return(F)
    }
  })
  
  return(isMins)
}

eqtls_wb[minLocal.cand==T,minLocal:=is.minLocal(tss_distance,pval_nominal),by=c("gene_id")]
eqtls_wb[minLocal==T] #158904 min local

eqtls_wb[,n.local.gene:=sum(minLocal,na.rm = T),by="gene_id"]
summary(unique(eqtls_wb,by="gene_id")$n.local.gene)
  # Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
  #  0.00    2.00    6.00   12.86   15.00  320.00 


# 2) find eQTR, i.e regions at +/-500pb around min local iteractively expanded if eQTLs at +/-2000pb and log10(pval) difference < 4
findeQTR<-function(minLocal,pvals,dists,seed=2000,log10.pval.dt.thr=4,length.max=10000,length.min=1000){
  reg<-0
  length.seed.reg<-1
  eQTRs<-rep(NA,length(pvals))
  iMinsLocals<-which(minLocal==T)
  for(i in iMinsLocals){
    #pour chaque minlocal, si pas deja affectÃ© a une reg : 
    if(!(i%in%which(!is.na(eQTRs)))){
      # => une nouvelle reg, avec un nouveau pval.thr de seed : 
      reg<-reg+1
      
      pvals.thr<-10^(log10(pvals[i])+log10.pval.dt.thr)
      # si eQTL a une dist <seed =>  inclus dans reg
      idansSeedsReg<-which(dists>(dists[i]-seed)&dists<(dists[i]+seed))
      idansReg<-idansSeedsReg
      eQTRs[idansReg]<-reg
      #pour chaque nouvelle reg si length.seed.reg < length.max, regarde si eQTL a une dist <seed, si oui, :
      length.seed.reg<- max(dists[idansSeedsReg])-min(dists[idansSeedsReg])
      
      if(length.seed.reg>0){
        while(length.seed.reg<length.max){
          
          start.seed.reg<-min(dists[idansSeedsReg])
          end.seed.reg<-max(dists[idansSeedsReg])
          candIdansSeedsReg<-setdiff(which(dists>(start.seed.reg-seed)&dists<(end.seed.reg+seed)),idansSeedsReg)
          if(length(candIdansSeedsReg)>0){
            
            #si length.seed.reg <length.min => inclus dans seed.reg
            if(length.seed.reg<length.min){
              idansSeedsReg<-c(idansSeedsReg,candIdansSeedsReg)
              eQTRs[idansSeedsReg]<-reg
              start.seed.reg<-min(dists[idansSeedsReg])
              end.seed.reg<-max(dists[idansSeedsReg])
              length.seed.reg<-end.seed.reg-start.seed.reg
            }
            #si length.seed.reg > length.min => inclus mais stopseed si diff pval min <log10.pval.dt.thr
            else{
              #ajout eQTL dans reg mais pas oblogatoirement dans seeed : 
              idansReg<-c(idansSeedsReg,candIdansSeedsReg)
              
              #inclus dans seed si assez signif :
              sigIdansSeedsReg<-candIdansSeedsReg[which(pvals[candIdansSeedsReg]<pvals.thr)]
              #sil y a des signifs on continue le seed :
              if(length(sigIdansSeedsReg)>0){
                idansSeedsReg<-c(idansSeedsReg,sigIdansSeedsReg)
                eQTRs[idansSeedsReg]<-reg
                start.seed.reg<-min(dists[idansSeedsReg])
                end.seed.reg<-max(dists[idansSeedsReg])
                length.seed.reg<-end.seed.reg-start.seed.reg
              }else{
                #sinon, on inclus les eQTL dans reg et on stop seed
                eQTRs[idansReg]<-reg
                length.seed.reg<-length.max
                
              }
              
            }
            
            
          }else{
            length.seed.reg<-length.max
          }
        }
        #sinon, on inclus les eQTL dans reg et on stop seed
        eQTRs[idansReg]<-reg
        
        
      }
      
    }
    
    
  }
  
  
  
  
  return(eQTRs)
}

#cluster snps by region
eqtls_wb[,eQTR:=findeQTR(minLocal,pval_nominal,tss_distance,seed = 2000,log10.pval.dt.thr = 4,length.max=10000,length.min = 1000),by="gene_id"]
eqtls_wb[!is.na(eQTR)] #1.2M / 2.4M eQTL on eQTR

eqtls_wb[!is.na(eQTR),start.eQTR:=min(pos),by=.(gene_id,eQTR)]
eqtls_wb[!is.na(eQTR),end.eQTR:=max(pos),by=.(gene_id,eQTR)]

#n eQTL in eQTR
eqtls_wb[!is.na(eQTR),n.eQTL:=.N,by=.(gene_id,eQTR)]
summary(unique(eqtls_wb[!is.na(eQTR)],by=c("gene_id","eQTR"))$n.eQTL)
 # Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
 #   1.00    2.00    6.00   10.73   13.00  689.00 

#for eQTR with only 1 eQTL, increase regions at +/-500pb from the eQTL
eqtls_wb[n.eQTL==1,start.eQTR:=start.eQTR-500]
eqtls_wb[n.eQTL==1,end.eQTR:=end.eQTR+500]



#translate gene_id in gene symbol
trans<-TransEnsemblVerstoSymbol(unique(eqtls_wb$gene_id))
trans<-trans[,gene_id:=ensembl_gene_id_version][,gene:=hgnc_symbol][,-c("ensembl_gene_id_version","hgnc_symbol")]
trans[gene=="",gene:=NA]
eqtls_wb<-merge(eqtls_wb,trans,by="gene_id")
unique(eqtls_wb,by="gene_id")[is.na(gene)]#1261/3179

trans2<-TransEnsembltoSymbol(str_extract(unique(eqtls_wb$gene_id),"ENSG[0-9]+"))
trans2[hgnc_symbol==""] #1261 too

#summary :
#eqtr by gene
eqtls_wb[,n.eqtr.gene:=length(unique(na.omit(eQTR))),by="gene_id"]

summary(unique(eqtls_wb,by='gene_id')$n.eqtr.gene)
  # Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
  # 0.000   2.000   5.000   8.894  11.000 178.000 

eqtls_wb[n.eqtr.gene==178] # ENSG00000281831.1 , https://www.genecards.org/cgi-bin/carddisp.pl?gene=HCP5B
#lncrna gene

eqtls_wb[n.eqtr.gene==0] #936 eQTL so negligeable

fwrite(eqtls_wb,fp(out,"eqtls_wb_with_eQTR.csv.gz"))

#eQTR_id for eQTR link to gene symbol and save in a bed files :
eqtls_wb[!is.na(gene),eqtr_id:=paste(gene,eQTR,sep="-")]
eqtrs_symbol<-unique(eqtls_wb[!is.na(eQTR)&!is.na(gene),.(chr,start.eQTR,end.eQTR,eqtr_id,gene)])
eqtrs_symbol#20k eqtr
unique(eqtrs_symbol,by="gene") #1800 genes
fwrite(eqtrs_symbol,fp(out,"whole_blood_eQTR_symbol.bed"),sep = "\t",col.names = F)

#II) With meta eqtls
eqtls_wb2<-fread("ref/eQTL/whole_blood_and_tissue_wide_eQTR_hgnc_gene_hg19.bed",sep="\t")
length(unique(eqtls_wb2$gene))

#III) merge eQTRs
