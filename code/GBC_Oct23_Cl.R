rm(list=ls()) #clear the environment

commanarg <- commandArgs(trailingOnly = TRUE)
i1 <- as.numeric(commanarg[1])
i2 <- as.numeric(commanarg[2])

library(ggplot2)
library(CGEN)
library(GenomicRanges)
library(rtracklayer)
setwd("../data")
small.head <- function(x){
  return(x[1:20,1:20])
}


cgen.logit <- function(pheno,geno,i){
  
  snp <- as.numeric(geno[,i])
  data <- cbind(pheno,as.numeric(snp))
  colnames(data) <- c(colnames(pheno),"snp")
  model <- snp.logistic(data,"affected","snp",main.vars=~age+gender+EV3+EV4+EV5+EV8+EV9)
  model.sum <- getSummary(model)
  return(c(i,model.sum$UML[9,4],model.sum$CML[9,4],model.sum$EB[9,4]))
}


pheno <- read.table("pheno_ev",header=T,sep="\t")  ###read the phenotype file
fam <- read.table("GBC_final.fam")    ###read the fam file
### change the column affected_status from "case" and "control" into 1 and 0
pheno$X <- NULL
pheno$affected_status.1 <- NULL
pheno$affected <- 0
pheno$affected[pheno$affected_status=="Case"] <- 1
###create my own pheno.txt
my.pheno <- fam[,1:2]
my.pheno <- cbind(my.pheno,pheno$affected,rbinom(nrow(my.pheno),1,0.5),rbinom(nrow(my.pheno),1,0.5))
write.table(my.pheno,"pheno.txt",quote=F,col.names=F,row.names = F)


EV.pheno <- pheno[,17:26]
pheno.EV.cor <- apply(EV.pheno,2,function(x) {cor.test(pheno$affected,x)})

ptm <- proc.time()
all.geno <- read.plink("GBC_final.bed","GBC_final.bim","GBC_final.fam")
proc.time()-ptm
geno <- all.geno$genotypes
fam <- all.geno$fam
map <- all.geno$map
loop.seq <- c(i1:i2)

cgen.logit.result<-sapply(loop.seq,function(x){cgen.logit(pheno,geno,x) })
setwd("../result")
save(cgen.logit.result,file=paste0("i1_i2.rda"))


