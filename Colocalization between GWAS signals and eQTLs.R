# Colocalization
if(!require("remotes"))
  install.packages("remotes")
  install.packages("dplyr")
library(remotes)
install_github("chr1swallace/coloc",build_vignettes=TRUE)
library("coloc")
library(dplyr)

setwd("<SET_YOUR_DIR>")
gwas <- read.table(file="GWAS.txt", header=F,stringsAsFactors=FALSE)
colnames(gwas)<-c("chr","rs_id", "pos", "A1", "A2", "N","AF1", "beta", "SE", "pval_nominal")
gwas$varbeta <- gwas$SE^2
gwas_sort <- gwas[order(gwas$pval_nominal), ]
eqtl <- read.table(file="eQTL.txt", header=T)
colnames(eqtl)<-c("gene_id","rs_id", "dst", "ma_samples", "ma_count", "maf","pval_nominal", "slope", "slope_se")
input <- merge(eqtl, gwas, by="rs_id", all=FALSE, suffixes=c("_eqtl","_gwas"))
head(input)

library(remotes)
library(coloc)
library(dplyr)
result <- coloc.abf(dataset1=list(pvalues=input$pval_nominal_gwas, type="quant", s=0.33, N=nrow(gwas)), dataset2=list(pvalues=input$pval_nominal_eqtl, type="quant", N=nrow(eqtl)), MAF=input$maf)
result$summary #若PP.H4.abf >0.75,则认为该区域存在共定位证据
result$results
need_result=result$results %>% filter(SNP.PP.H4 > 0.70)



