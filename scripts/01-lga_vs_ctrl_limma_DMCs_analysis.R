source("scripts/utils/new_utils.R")
source("scripts/utils/methyl_utils.R")

out<-"outputs/01-lga_vs_ctrl_limma_DMCs_analysis"
dir.create(out)


meth_file<-here("datasets/cd34/CD34_angle_119_noEmptyLocis_withConfScore_withoutChrXY.txt")
mtd_file<-here("datasets/cd34/cleaned_batch_CD34_library_date_220620.csv")


#Loading

mtd<-fread(mtd_file,
            select = c("sample","Group_name","Gender","Group_Sex",
                       "Weight..g.","Weight.at.term..lbs.","GA..wk.","Length..cm.","Mat.Age","HC..cm.","PI","SeqDepth",
                       "latino","Preterm","GDM","Drugs","Etoh", "Smoking",
                       "Race","Labor","batch","date","DNA.extraction","sequencing",
                       "Library_Complexity","Group_Complexity","Group_Complexity_Fac","GroupBatch_Complexity","GroupBatch_Complexity_Fac"),
            col.names = c("sample","group","sex","group_sex",
                       "weight.g","weight_at_term.lbs","gest.age.wk","length.cm","mat.age","head_circumf.cm.","ponderal_index","seq.depth",
                       "latino","preterm","GDM","drugs","etoh", "smoking",
                       "ethnicity","lab","batch","date","DNA.extraction","sequencing",
                       "library_complexity","group_complexity","group_complexity_fac","groupbatch_complexity","groupbatch_complexity_fac"))

mtd[,group:=ifelse(group=="L","LGA",ifelse(group=="I","IUGR","CTRL"))]
mtd[,group_sex:=paste(group,sex,sep = "_")]
mtd[,year:=str_extract(date,"20[0-9]+$")]

bool_vars<-c("latino","preterm","GDM","drugs","etoh", "smoking")
mtd[,(bool_vars):=lapply(.SD,as.logical),.SDcols=bool_vars]

categorical_vars<-c("group","sex","group_sex","ethnicity","lab","batch","date","DNA.extraction","sequencing","year")
mtd[,(c(bool_vars,categorical_vars)):=lapply(.SD,as.factor),.SDcols=c(bool_vars,categorical_vars)]


fwrite(mtd,"datasets/cd34/metadata_140421.csv",sep=";")
meth<-fread(meth_file,
            select = c("id","chr","start",mtd$sample,"msp1c","confidenceScore"),
            col.names = c("cpg_id","chr","pos",mtd$sample,"msp1c","confidence_score"))

all(mtd$sample%in%colnames(meth))
nrow(mtd) #119
table(mtd[,.(group,sex)])
#       sex
# group   F  M
#   CTRL 17 20
#   IUGR 20 20
#   LGA  21 21

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

meth #754931 CpGs after filtering

#Focus on Ctrl and LGA
mtd<-mtd[group%in%c("CTRL","LGA")]
meth<-meth[,.SD,.SDcols=c("cpg_id",mtd$sample)]

all(mtd$sample%in%colnames(meth))
nrow(mtd) #79
table(mtd[,.(group,sex)])
#       sex
# group   F  M
#   CTRL 17 20
#   LGA  21 21


#Significant covariates findings
meth_mat<-as.matrix(data.frame(meth,row.names = "cpg_id")[,mtd$sample])
pca<-prcomp(t(meth_mat))

pval_mat<-CorrelCovarPCs(pca =pca ,mtd,rngPCs =1:10,res = "pval") #batch, seq depth, sequencing,group_complexity,group , ethnicity

meth.influencing.vars<-rownames(pval_mat)[rowSums(pval_mat<0.01)>0]

fwrite(data.table(pval_mat,keep.rownames = "covar"),fp(out,"res_correl_covar_pcs.csv"),sep=";")

#

#Clinical na imputation


#Limma


#Gene score calculation


#Others approach calculation


#Validation Gene Score


