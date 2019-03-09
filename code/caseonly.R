setwd("../data")
##  build snp list for function GxE.scan
snp.list <- list()
snp.list$format <- "tped"
file.name <- "caseonlysnp.tped"

snp.list$file <- file.name
#snp.list$start.vec <- 1
#snp.list$stop.vec <- 10
snp.list$subject.list <- "caseonlysnp.tfam"
##  build phenotype list for function GxE.scan
f = "pheno.txt"
pheno.list <- list(file=f, header=1, delimiter="\t", id.var=c("familyID","ID"), response.var="affected") 
pheno.list$file.type <- 3
pheno.list$main.vars <- ~age+gender+EV3+EV4+EV5+EV8+EV9
pheno.list$int.vars <- ~stoneyn
##  build op list for funciton GxE.scan
op<- list()


result <- GxE.scan(snp.list,pheno.list,op)
caseonly_result <- read.table("GxE.scan.output.txt",header=T)