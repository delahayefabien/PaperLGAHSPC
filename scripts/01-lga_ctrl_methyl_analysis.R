source("scripts/utils/new_utils.R")
source("scripts/utils/methyl_utils.R")

out<-"outputs/gene_score_by_prom_enh"
dir.create(out)


meth_file<-here("datasets/cd34/CD34_angle_119_noEmptyLocis_withConfScore_withoutChrXY.txt")
mtd_file<-here("datasets/cd34/cleaned_batch_CD34_library_date_220620.csv")


#Loading

mtd<-fread(mtd_file)
meth<-fread(meth_file,
            select = c("id","chr","start",mtd$sample,"msp1c","confidenceScore"),
            col.names = c("cpg_id","chr","pos",mtd$sample,"msp1c","confidence_score"))

all(mtd$sample%in%colnames(meth))
nrow(mtd)
table(mtd$Group_Sex)
# C_F C_M I_F I_M L_F L_M 
#  17  20  20  20  21  21 

#data exploration
plot(density(na.omit(as.matrix(meth[,.SD,.SDcols=mtd$sample]))))
abline(v=10)

#CpG filtering
# add useful quality metrics 

meth[,n.na:=rowSums(is.na(.SD)),.SDcols=mtd$sample]
meth[,n.not.methyl:=rowSums(.SD>10,na.rm = T),.SDcols=mtd$sample]
meth[,pct.zero:=rowSums(.SD==0,na.rm = T)/length(mtd$sample),.SDcols=mtd$sample]
meth[,n.methyl.not.zero:=rowSums(.SD>0&.SD<10,na.rm = T),.SDcols=mtd$sample]
# + msp1count (msp1c) and confidencescore, already present in the data
#   msp1c : count for reference library digest with msp1
#   confidencescore : sum of all count for hpa2 libraries and msp1 library, normalized by library size

#   see 01A-CpG_filtering to see how we determine the threshold,

meth<-meth[msp1c>quantile(msp1c,0.125,na.rm=T)&n.na==0]

meth<-meth[n.not.methyl>4]
meth<-meth[confidence_score>quantile(confidence_score,0.2,na.rm=T)]
meth<-meth[!(pct.zero>0.7&n.methyl.not.zero==0)]

meth
#Focus on Ctrl and LGA

mtd<-mtd[Group_name%in%c("C","L")]
meth<-meth[,.SD,.SDcols=c("cpg_id",mtd$sample)]

all(mtd$sample%in%colnames(meth))
nrow(mtd)
table(mtd$Group_Sex)
# C_F C_M L_F L_M 
#  17  20  21  21 

#data exploration
plot(density(na.omit(as.matrix(meth[,.SD,.SDcols=mtd$sample]))))
abline(v=10)




#Significant covariates findings


#Clinical na imputation


#Limma


#Gene score calculation


#Others approach calculation


#Validation Gene Score


