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
mtd[,ethnicity:=factor(ethnicity,levels=unique(ethnicity))]
bool_vars<-c("latino","preterm","GDM","drugs","etoh", "smoking")
mtd[,(bool_vars):=lapply(.SD,as.logical),.SDcols=bool_vars]

categorical_vars<-c(bool_vars,"group","sex","group_sex","ethnicity","lab","batch","date","DNA.extraction","sequencing","year")
mtd[,(categorical_vars):=lapply(.SD,as.factor),.SDcols=categorical_vars]

numerical_vars<-c("weight.g","weight_at_term.lbs","gest.age.wk","length.cm","mat.age","head_circumf.cm.","ponderal_index","seq.depth","library_complexity","group_complexity","group_complexity_fac","groupbatch_complexity","groupbatch_complexity_fac")


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
fwrite(mtd,"datasets/cd34/metadata_cl_190421.csv",sep=";")

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
library(ggrepel)
meth_mat<-as.matrix(data.frame(meth,row.names = "cpg_id")[,mtd$sample])
pca<-prcomp(t(meth_mat))
pc_mtd<-merge(mtd,data.table(pca$x,keep.rownames = "sample"))

fwrite(pc_mtd,"datasets/cd34/metadata_pcs_cl_190421.csv",sep=";")

ggplot(pc_mtd)+geom_point(aes(x=PC1,y=PC2,col=group))
ggplot(pc_mtd)+geom_point(aes(x=PC1,y=PC2,col=ponderal_index))
ggplot(pc_mtd)+geom_point(aes(x=PC1,y=PC3,col=GDM))
#PC2, unexplained variation so include in the model

ggplot(pc_mtd)+geom_point(aes(x=PC1,y=PC6,col=mat.age))

ggplot(pc_mtd)+geom_point(aes(x=PC2,y=PC9,col=latino))

mtd[is.na(mat.age)]
mtd[is.na(latino)]

ggplot(pc_mtd,aes(x=PC1,y=PC2,col=GDM))+geom_point()+ggrepel::geom_label_repel(aes(label=sample))
mtd[sample%in%c("CBP186","CBP253","CBP252","CBP205")]

ggplot(pc_mtd)+geom_point(aes(x=PC1,y=PC2,col=library_complexity))


pc_mtd[is.na(ponderal_index)]
pval_mat<-CorrelCovarPCs(pca =pca ,mtd,rngPCs =1:30,res = "pval",seuilP = 0.1) #batch, seq depth,group_complexity,group , ethnicity

r2_mat<-CorrelCovarPCs(pca =pca ,mtd,rngPCs =1:30,res = "r2",seuilP = 0.1) #batch, seq depth,group_complexity,group , ethnicity

meth.influencing.vars<-rownames(pval_mat)[rowSums(pval_mat<0.01)>0]
meth.influencing.vars<-c(meth.influencing.vars,"mat.age","latino","group")
fwrite(data.table(pval_mat,keep.rownames = "covar"),fp(out,"res_correl_covar_pcs.csv"),sep=";")

# correl between covars and group
infl.vars.fac<-meth.influencing.vars[meth.influencing.vars%in%categorical_vars]
infl.vars.num<-meth.influencing.vars[meth.influencing.vars%in%numerical_vars]

cor_nums<-sapply(infl.vars.num,function(var1){
  r2s<-sapply(infl.vars.num,function(var2){
    f<-as.formula(paste(var1,"~",var2))
    return(summary(lm(f,mtd))$adj.r.squared)
    })
  return(r2s)
  })
cor_nums<-data.matrix(cor_nums)
pheatmap(cor_nums,display_numbers = T,cluster_rows = F,cluster_cols = F) 

pvals_num_fac<-sapply(infl.vars.num,function(var1){
  pvals<-sapply(infl.vars.fac,function(var2){
    f<-as.formula(paste(var1,"~",var2))
    return(anova(lm(f,mtd))$Pr[1])
    })
  return(pvals)
  })
pvals_num_fac<-data.matrix(pvals_num_fac)
plotPvalsHeatMap(pvals_num_fac) 
summary(lm(mtd$seq.depth~mtd$group))

pvals_fac<-sapply(infl.vars.fac,function(var1){
  pvals<-sapply(infl.vars.fac,function(var2)correl(mtd[,.SD,.SDcols=c(var1,var2)],verbose = T))
  return(pvals)
  })
pvals_fac<-data.matrix(pvals_fac)
pvals_fac[is.na(pvals_fac)]<-1
pheatmap(-log10(pvals_fac),display_numbers = T,cluster_rows = F,cluster_cols = F) #no correl
table(mtd[,.(group,GDM)])

#all this covars necessary to explain PC1?
summary(lm(PC1~group_complexity_fac+group+mat.age+latino,data = pc_mtd)) #seq.depth and mat age not sig in PC1
summary(lm(PC1~group_complexity_fac+group+mat.age+latino+seq.depth,data = pc_mtd)) #seq.depth to rm
summary(lm(PC1~group_complexity_fac+group+latino,data = pc_mtd)) #mat.age does not help to explain PC1 but explain well PC6 and PC26
ggplot(pc_mtd)+geom_point(aes(x=PC6,y=PC26,col=mat.age))
#and corralted with pct0
pc_mtd[,pct_zero:=colSums(meth_mat[,pc_mtd$sample]==0)/nrow(meth_mat)]
ggplot(pc_mtd)+geom_point(aes(x=mat.age,y=pct_zero,col=mat.age))


#include group_sex in the model insteag of group ?
summary(lm(PC1~group_complexity_fac+group_sex+latino,data = pc_mtd)) #can see that effect ++ in lga F than Lga M, and CtrlM not sig



# so include batch,mat.age, group_complexity_fac,group_sex and PC2 in the model 
vars_to_include<-c("batch","mat.age","group_complexity_fac","group","latino","PC2")


#DATA MODELING and DMC analysis with limma
# rm samples without all necessary clinical infos
mtd_f<-pc_mtd[,to_keep:=rowSums(is.na(.SD))==0,.SDcols=vars_to_include][to_keep==T]

formule<- ~0 + group_sex   + batch+ group_complexity_fac +mat.age  + latino + PC2

mtd_f[,group_sex:=factor(group_sex,levels = unique(mtd_f$group_sex))]
design<-model.matrix(formule,data = data.frame(mtd_f,row.names = "sample"))
fit <- lmFit(data.frame(meth,row.names = "cpg_id")[,mtd_f$sample], design)

cont.matrix <- makeContrasts(C.L = "(group_sexCTRL_F+group_sexCTRL_M)-(group_sexLGA_F+group_sexLGA_M)",
                             levels=design)


fit2  <- contrasts.fit(fit, cont.matrix)
fit2  <- eBayes(fit2)

res<-data.table(topTable(fit2,coef = "C.L",n = Inf),keep.rownames = "cpg_id")
res[adj.P.Val<=0.05] #43 cpgs
res[P.Value<0.001&abs(logFC)>30] #4774 cpgs
fwrite(res,fp(out,"res_limma.tsv.gz"),sep="\t")

p<-ggplot(res)+
  geom_point(aes(x=logFC,y=-log10(P.Value),col=P.Value<0.001&abs(logFC)>30))+
  scale_color_manual(values = c("grey","red"))
ggsave(fp(out,"volcano_plot.png"),plot=p,height=5,width=7)





