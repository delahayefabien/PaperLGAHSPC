# methyl

Analyst: Alexandre Pelletier
##Config
``` bash
git config --global user.email "apelletier@egid.local"
git config --global user.name "Alexandre Pelletier"
```
## install dependencies
``` r
renv::hydrate()
renv::install(c("YuLab-SMU/clusterProfiler"))
renv::install("missMDA")
renv::snapshot() #git commit
```

## add symbolic link to datasets, ref and old project directory

``` bash
cd 
ln -s /disks/DATATMP/PhD_AlexandrePelletier/methyl/methyl/datasets datasets
ln -s /disks/DATATMP/PhD_AlexandrePelletier/methyl/methyl/ref ref
ln -s /disks/DATATMP/HSPC_EpiStress_AlexandrePelletier/methyl/ old 
```

<!--
## Design

``` bash
nohup Rscript scripts/01-design.R > logs/01.log &
```

## Quality Control

``` bash
nohup Rscript -e 'rmarkdown::render(input = here::here("scripts", "02-qc.Rmd"), output_file = "QC.html", output_dir = here::here("reports"), encoding = "UTF-8")' > logs/02.log &
```

## Statistical Analyses

``` bash
nohup Rscript scripts/03-analysis.R > logs/03.log &
```

## Meeting Slides

### 2021-02-22

``` bash
nohup Rscript -e 'rmarkdown::render(input = here::here("scripts", "20210222-meeting.Rmd"), output_dir = here::here("reports"))' > logs/20210222-meeting.log &
```
-->
