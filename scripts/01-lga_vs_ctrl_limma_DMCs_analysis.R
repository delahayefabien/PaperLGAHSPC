source("scripts/utils/new_utils.R")
source("scripts/utils/methyl_utils.R")
library(limma)
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
mtd[ethnicity%in%c("Declined",""," "),ethnicity:=NA]
bool_vars<-c("latino","preterm","GDM","drugs","etoh", "smoking")
mtd[,(bool_vars):=lapply(.SD,as.logical),.SDcols=bool_vars]

categorical_vars<-c(bool_vars,"group","sex","group_sex","ethnicity","lab","batch","date","DNA.extraction","sequencing","year")
mtd[,(categorical_vars):=lapply(.SD,as.factor),.SDcols=categorical_vars]


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



#FOCUS ON CTRL LGA
mtd<-mtd[group%in%c("CTRL","LGA")]
meth<-meth[,.SD,.SDcols=c("cpg_id",mtd$sample)]

all(mtd$sample%in%colnames(meth))
nrow(mtd) #79
table(mtd[,.(group,sex)])
#       sex
# group   F  M
#   CTRL 17 20
#   LGA  21 21


#COVARIATES ANALYSIS
# check no vars too much unique individual by levels (singleton)

lapply(mtd[,.SD,.SDcols=categorical_vars],function(x)sum(table(x)==1)) #exclude sequencing and date

mtd<-mtd[,-c("sequencing","date")]

meth_mat<-as.matrix(data.frame(meth,row.names = "cpg_id")[,mtd$sample])
pca<-prcomp(t(meth_mat))

pval_mat<-CorrelCovarPCs(pca =pca ,mtd,rngPCs =1:10,res = "pval",seuilP = 0.05) #batch, seq depth,group_complexity,group , ethnicity

meth.influencing.vars<-rownames(pval_mat)[rowSums(pval_mat<0.01)>0]

fwrite(data.table(pval_mat,keep.rownames = "covar"),fp(out,"res_correl_covar_pcs.csv"),sep=";")

# correl between covars and group
infl.vars.fac<-meth.influencing.vars[meth.influencing.vars%in%categorical_vars]
infl.vars.num<-setdiff(meth.influencing.vars,infl.vars.fac)

infl.vars.fac<-union(infl.vars.fac,"group")
cor_nums<-sapply(infl.vars.num,function(var1){
  r2s<-sapply(infl.vars.num,function(var2){
    f<-as.formula(paste(var1,"~",var2))
    return(summary(lm(f,mtd))$adj.r.squared)
    })
  return(r2s)
  })
cor_nums<-data.matrix(cor_nums)
pheatmap(cor_nums,display_numbers = T,cluster_rows = F,cluster_cols = F) #group_complexity fac and seq depth

pvals_num_fac<-sapply(infl.vars.num,function(var1){
  pvals<-sapply(infl.vars.fac,function(var2){
    f<-as.formula(paste(var1,"~",var2))
    return(anova(lm(f,mtd))$Pr[1])
    })
  return(pvals)
  })
pvals_num_fac<-data.matrix(pvals_num_fac)
plotPvalsHeatMap(pvals_num_fac) #group_complexity fac and seq depth

pvals_fac<-sapply(infl.vars.fac,function(var1){
  pvals<-sapply(infl.vars.fac,function(var2)correl(mtd[,.SD,.SDcols=c(var1,var2)],verbose = T))
  return(pvals)
  })
pvals_fac<-data.matrix(pvals_fac)
pvals_fac[is.na(pvals_fac)]<-1
pheatmap(-log10(pvals_fac),display_numbers = T,cluster_rows = F,cluster_cols = F) #no correl
table(mtd[,.(group,GDM)])

# so include batch, ethnicity, GDM, group_complexity_fac,sequencing, and seq.depth in the model 
vars_to_include<-c("batch","ethnicity","GDM","group_complexity_fac","seq.depth")


#DATA MODELING and DMC analysis with limma
#Clinical na imputation
library(FactoMineR)
library(missMDA)
sum(rowSums(sapply(mtd[,.SD,.SDcols=vars_to_include], function(x)is.na(x)))>0) #8/79 with na values
mtd_cat_df<-data.frame(mtd[,.SD,.SDcols=c("sample",intersect(colnames(mtd),categorical_vars))],row.names = "sample")
res.mca<-MCA(mtd_cat_df, quali.sup = which(!colnames(mtd_cat_df)%in%intersect(vars_to_include,categorical_vars)))
res.impute <- imputeMCA(mtd_cat_df[,which(colnames(mtd_cat_df)%in%intersect(vars_to_include,categorical_vars))], ncp=2)
mtd_imp<-data.table(res.impute$comp[,c("ethnicity","GDM")],keep.rownames = "sample")
colnames(mtd_imp)[2:3]<-paste0(colnames(mtd_imp)[2:3],"_imp")

mtd<-merge(mtd,mtd_imp,by="sample")
mtd[,.(ethnicity,ethnicity_imp,GDM,GDM_imp)]

#Limma
formule<- ~0 + group  + ethnicity_imp + GDM_imp + batch+ group_complexity_fac + seq.depth
mtd[,group:=factor(group,levels = unique(mtd$group))]
design<-model.matrix(formule,data = data.frame(mtd,row.names = "sample"))
fit <- lmFit(data.frame(meth,row.names = "cpg_id")[,mtd$sample], design)

cont.matrix <- makeContrasts(C.L = "groupCTRL-groupLGA",
                             levels=design)


fit2  <- contrasts.fit(fit, cont.matrix)
fit2  <- eBayes(fit2)

res<-data.table(topTable(fit2,coef = "C.L",n = Inf),keep.rownames = "cpg_id")

fwrite(res,fp(out,"res_limma.tsv.gz"),sep="\t")

p<-ggplot(res)+
  geom_point(aes(x=logFC,y=-log10(P.Value),col=P.Value<0.001&abs(logFC)>20))+
  scale_color_manual(values = c("grey","red"))
ggsave(fp(out,"volcano_plot_cohort2.png"),plot=p,height=5,width=7)




