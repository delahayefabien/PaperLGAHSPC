#Concept : Find comethylated genes across samples 

#postulats : 
  #des TFs peuvent réguler l'activité de la DNMT/TET autour des sites/motifs sur lesquelles ils se fixent (https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6923324/)

#question biologique : y a til des groupes (modules) de genes cométhylés a travers les échantillons, qui seraient la signature de l'action d'un TF ?
#question subsidiaire : ces modules sont ils associés avec un trait phénotypique particulier ?
# ==> question technique : Comment identifier des modules de genes co-méthylés ?

#PLAN DANALYSE
  #main steps (Verifier/Valider chaque étape avant de passer à la suivante)
    #1) aggréger la methylation des CpGs à l'échelle du gène (genescore)
      #a) lier les CpGs au genes
      #b) assigner un poids pour chaque lien CpG-gene
      #c) pour chaque gène, faire la moyenne pondéré des valeurs de méthylation des CpGs liés à lui

    #2) cross correlation de ce genescore à travers les samples
      #a) caractérisation du score pour choisir le bon test de correlation
      #b) analyse de correlation
      #c) clusteriser la matrice de correlation pour identifier des modules

    #3) caractérisations/validation des modules 
      #a) pathway/Biological process/TF signature enrichment analysis in these modules
      #b) TF binding motifs analysis 
      #c) correlation avec des traits phénotypiques (Birth weight, maternal age..)
      #d) Integration avec les modules de co-expression identifies en scRNA-seq




#ANALYSIS
source("scripts/utils/new_utils.R")
out<-"../methyl/outputs/04-comethylation_analysis"
dir.create(out)
    #1) aggréger la methylation des CpGs à l'échelle du gène (genescore)
      #a) lier les CpGs au genes et b) assigner un poids pour chaque lien CpG-gene
#see 02A 
cpgs_weight<-fread("../methyl/ref/2020-06-29_All_CpG-Gene_links.csv")
cpgs_weight[,cpg_weight:=RegWeight+LinksWeight]

      #c) pour chaque gène, faire la moyenne pondéré des valeurs de méthylation des CpGs liés à lui
meth_df<-merge(meth_df,cpgs_weight[,.(locisID,gene,cpg_weight)],by="locisID")
meth_df[,gene_meth:=sum(methylation*cpg_weight)/sum(cpg_weight),by=c("sample","gene")]
fwrite(meth_df,fp(out,"methylation_sum_by_gene_datatable.tsv.gz"),sep="\t")
meth_df<-fread(fp(out,"methylation_sum_by_gene_datatable.tsv.gz"),sep="\t")
#caracterisation score : 
methg<-unique(meth_df,by=c("sample","gene"))

plot(density(methg$gene_meth))

#tecnic bias by sample ?
ggplot(methg)+geom_density(aes(gene_meth))+facet_wrap("sample") #yes

#QC : rm genes with 0 values
methg0<-methg[!(gene%in%gene[gene_meth==0])]
length(unique(methg0$gene)) #7401

plot(density(methg0$gene_meth))
#tecnic bias by sample ?
ggplot(methg0)+geom_density(aes(gene_meth))+facet_wrap("sample") #yes

methg0[,gene_meth_scaled:=scale(gene_meth),by=c("sample")]
ggplot(methg0)+geom_density(aes(gene_meth_scaled))+facet_wrap("sample") #ok

    #2) cross correlation de ce genescore à travers les samples
      #a) caractérisation du score pour choisir le bon test de correlation
#between 0-100 (-4 and 8 for the scaled version), 
#moyenne pondéré of the percentage of HpA2 (methylation sensitive) cleaving counts compared to Msp1 (reference, unsensitive) cleaving counts
#=> based on https://statsandr.com/blog/correlation-coefficient-and-correlation-test-in-r/, use spearman correl because :
#"Spearman correlation (which is actually similar to Pearson but based on the ranked values for each variable rather than on the raw data) is often used to evaluate relationships involving qualitative ordinal variables or quantitative variables if the link is partially linear"

      #b) analyse de correlation with spearman
#use Hmisc because calculate correlation value and significance in one 
renv::install("Hmisc")
library(Hmisc)
#need a matrix obs/var
methg0_mat<-as.matrix(dcast(methg0,sample~gene,value.var = "gene_meth_scaled"),rownames = "sample")
dim(methg0_mat) #108 7401
7401^2 #54M de tests de correl
head(methg0_mat[,1:10])
res<-rcorr(methg0_mat,type="spearman")
res$P.adj<-p.adjust(res$P,method = "BH")
sum(res$P.adj<0.05,na.rm = T) #3,7M of signif correl

saveRDS(res,fp(out,"gene_correl_spearman.rds"))

108^2 #12k tests de correl
res_sample<-rcorr(t(methg0_mat),type="spearman")
res_sample$P.adj<-p.adjust(res_sample$P,method = "BH")
sum(res_sample$P.adj<0.05,na.rm = T) #12k of signif correl
saveRDS(res_sample,fp(out,"sample_correl_spearman.rds"))

      #c) clusteriser la matrice de correlation pour identifier des modules
res<-readRDS(fp(out,"gene_correl_spearman.rds"))
res_sample<-readRDS(fp(out,"sample_correl_spearman.rds"))

#1) test with pheatmap : is there module of genes ?
renv::install("pheatmap")
mtd<-fread("datasets/cd34/cleaned_batch_CD34_library_date_220620.csv")
library(pheatmap)
pheatmap(
  t(methg0_mat), 
  annotation_col = data.frame(mtd,row.names = "sample")[colnames(methg0_mat),c("Group_Sex","Group_Complexity_Fac")],
  clustering_distance_cols = as.dist(1 - res_sample$r),
  clustering_distance_rows = as.dist(1 - res$r),
  show_rownames = F,show_colnames = F
  ) #doesnt work


#2) test with louvain alg/meth of neurons jc paper
#i) meth based on neurons jc paper
#use MEGENA ()
#-test for planarity with PMFG 
renv::install("songw01/MEGENA")

library(MEGENA)
# input parameters
n.cores <- 5; # number of cores/threads to call for PCP
doPar <-TRUE; # do we want to parallelize?
method = "spearman" # method for correlation. either pearson or spearman. 
FDR.cutoff = 0.05 # FDR threshold to define significant correlations upon shuffling samples. 
module.pval = 0.05 # module significance p-value. Recommended is 0.05. 
hub.pval = 0.05 # connectivity significance p-value based random tetrahedral networks
cor.perm = 10; # number of permutations for calculating FDRs for all correlation pairs. 
hub.perm = 100; # number of permutations for calculating connectivity significance p-value. 
# annotation to be done on the downstream
annot.table=NULL
id.col = 1
symbol.col= 2
###########

ijw <- calculate.correlation(t(methg0_mat),
                             doPerm = cor.perm,
                             output.corTable = FALSE,
                             output.permFDR = FALSE,
                             method=method)

# calculate PFN
#In this step, Planar Filtered Network (PFN) is calculated by taking significant correlation pairs, ijw. In the case of utilizing a different similarity measure, one can independently format the results into 3-column data frame with column names c("row","col","weight"), and make sure the weight column ranges within 0 to 1. Using this as an input to calculate.PFN() will work just as fine. 
#### register multiple cores if needed: note that set.parallel.backend() is deprecated. 
run.par = doPar & (getDoParWorkers() == 1) 
if (run.par)
{
  cl <- parallel::makeCluster(n.cores)
  registerDoParallel(cl)
  # check how many workers are there
  cat(paste("number of cores to use:",getDoParWorkers(),"\n",sep = ""))
}
##### calculate PFN
el <- calculate.PFN(ijw[,1:3],doPar = doPar,num.cores = n.cores,keep.track = FALSE)
g <- graph.data.frame(el,directed = FALSE)
```

# perform clustering
#MCA clustering is performed to identify multiscale clustering analysis. "MEGENA.output" is the core output to be used in the down-stream analyses for summarization and plotting.

##### perform MCA clustering.
MEGENA.output <- do.MEGENA(g,
 mod.pval = module.pval,hub.pval = hub.pval,remove.unsig = TRUE,
 min.size = 10,max.size = vcount(g)/2,
 doPar = doPar,num.cores = n.cores,n.perm = hub.perm,
 save.output = FALSE)
###### unregister cores as these are not needed anymore.
if (getDoParWorkers() > 1)
{
  env <- foreach:::.foreachGlobals
  rm(list=ls(name=env), pos=env)
}
saveRDS(MEGENA.output,fp(out,"megena_spearman_pfn_mca_outputs.rds"))


    #3) caractérisations/validation des modules 
      #a) pathway/Biological process/TF signature enrichment analysis in these modules
      #b) TF binding motifs analysis 
      #c) correlation avec des traits phénotypiques (Birth weight, maternal age..)
      #d) Integration avec les modules de co-expression identifies en scRNA-seq



