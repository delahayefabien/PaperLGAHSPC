
source("scripts/utils/new_utils.R")
set.seed(1234)
old_path <- Sys.getenv("PATH")
Sys.setenv(PATH = paste(old_path, "/home/apelletier/HSPC_EpiStress/Data/tools/homer/bin/", sep = ":"))

out0<-"outputs/03B-motif_analysis_DMC_pval0.001_FC30"
dir.create(out0)
n_perm<-100

res<-fread("outputs/01-lga_vs_ctrl_limma_DMCs_analysis/res_limma.tsv.gz")

n_cpgs<-nrow(res[P.Value<0.001&abs(logFC)>25])

cpgs<-fread("datasets/cd34/CD34_angle_119_noEmptyLocis_withConfScore_withoutChrXY.txt",
            select = c("id","chr","start"),
            col.names = c("cpg_id","chr","pos"))
res<-merge(res,cpgs)
res[,start:=pos-20][start<1,start:=1][,end:=pos+20]
res<-res[!is.na(chr)&!is.na(pos)]


for(i in 1:n_perm){
  
  message("perm ",i,"/",n_perm)
  out<-fp(out0,paste0("perm",i))
  dir.create(out)
  
  bed_path<-fp(out,paste0("random",i,".bed"))
  
  random_cpgs<-res[cpg_id%in%sample(cpg_id,n_cpgs),.(chr,start,end)][order(chr,start)]

  fwrite(random_cpgs[,.(chr,start,end)][order(chr,start)],
       bed_path,sep="\t",
       col.names=F)
  
  cmd<-paste("findMotifsGenome.pl",bed_path,"hg19",out,"-size 41 -p 4 -cpg")
  system(cmd)
  
  res_k<-fread(fp(out,"knownResults.txt"),
           select = c(1,2,3,5,6,7,8,9),
           col.names = c("motif","consensus","pval","padj","n_dmc_with_motif","pct_dmc_with_motif","n_background_with_motif","pct_background_with_motif"))
  res_k[,pct_dmc_with_motif:=as.numeric(str_remove(pct_dmc_with_motif,"%"))]
  res_k[,pct_background_with_motif:=as.numeric(str_remove(pct_background_with_motif,"%"))]
  res_k[,permut:=i]
  
  if(i==1){
    res_all<-copy(res_k)
  }else{
    res_all<-rbind(res_all,res_k)
  }
  
  if(i>5){
    message("removing res")
    unlink(out,recursive = T)
  }
}

fwrite(res_all,fp(out0,"res_known_motif_all_perm.csv"))
message("Success !")


