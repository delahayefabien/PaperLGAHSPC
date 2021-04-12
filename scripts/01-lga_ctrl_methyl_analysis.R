source("scripts/utils/new_utils.R")
source("scripts/utils/methyl_utils.R")

out<-"outputs/gene_score_by_prom_enh"
dir.create(out)


meth_file<-here("datasets/cd34/CD34_angle_119_noEmptyLocis_withConfScore_withoutChrXY.txt")
mtd_file<-here("datasets/cd34/cleaned_batch_CD34_library_date_220620.csv")


#Loading

mtd<-fread(mtd_file)[Group_name%in%c("C","L")]
meth<-fread(meth_file,
            select = c("id","chr","start","msp1c",mtd$sample),
            col.names = c("cpg_id","chr","pos","msp1c",mtd$sample))

all(mtd$sample%in%colnames(meth))
nrow(mtd)
table(mtd$Group_Sex)
# C_F C_M L_F L_M 
#  17  20  21  21 

#data exploration
plot(density(na.omit(as.matrix(meth[,.SD,.SDcols=mtd$sample]))))
abline(v=10)
#CpG filtering

# based on msp1 count, nb of na, pct zero and nb of high mehtylation score different of zero
# build some quality scores :  
  #nb of NA
meth[,n.na:=rowSums(is.na(.SD)),.SDcols=mtd$sample]
  #nb of zero
meth[,n.zero:=rowSums(.SD==0,na.rm = T),.SDcols=mtd$sample]

  #nb of high mehtylation score different of zero (cpg methylated with confidence)
meth[,n.conf.meth:=rowSums(.SD>0&.SD<15,na.rm = T),.SDcols=mtd$sample]


# find good threshold
  #msp1c
thr_to_test<-quantile(meth$msp1c,1:9/40,na.rm=T)

quals_meth<-sapply(thr_to_test,function(t){
  meth_f<-meth[msp1c<t]
  return(sum(meth_f[,rowSums(.SD>0&.SD<5,na.rm = T),.SDcols=mtd$sample])/nrow(meth_f))
  })
plot(thr_to_test,quals)

#Significant covariates findings


#Clinical na imputation


#Limma


#Gene score calculation


#Others approach calculation


#Validation Gene Score


