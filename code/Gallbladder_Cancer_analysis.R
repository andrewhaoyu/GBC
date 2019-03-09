##  Research goal is to analyze the association of GBC with Snp.
##  We used package CGEN to do the analysis
##  First step is to transform the bed data into tped data with plink
##so that the package CGEN could read it
##  ./plink bfile --GBC_final --recode --transpose --out GBC_final
##  Second step is to split the tped file into 50 small tped file to increase
##computation speed, we used linux command to do this
##  split -l 14298 GBC_final.tped GBC_final
##  use R to perform CGEN analysze
rm(list=ls()) #clear the environment
export(cache)
export(cache.project)
export(create.project)
export(get.project)
export(load.project)
export(reload.project)
export(require.package)
export(run.project)
export(show.project)
export(stub.tests)
export(test.project)
export(translate.dcf)
db.reader <- function(data.file, filename, variable.name)
{
  require.package('RSQLite')
  
  sqlite.driver <- dbDriver("SQLite")
  connection <- dbConnect(sqlite.driver,
                          dbname = filename)
  
  tables <- dbListTables(connection)
  for (table in tables)
  {
    message(paste('  Loading table:', table))
    
    data.parcel <- dbReadTable(connection,
                               table,
                               row.names = NULL)
    
    assign(clean.variable.name(table),
           data.parcel,
           envir = .TargetEnv)
  }
  
  disconnect.success <- dbDisconnect(connection)
  if (! disconnect.success)
  {
    warning(paste('Unable to disconnect from database:', filename))
  }
}


##  summary data analysis
##  read the data in R and change affected status from case and control 
##to 1 and 0
library(ggplot)
setwd("../data")
pheno <- read.table("pheno_ev",header=T,sep="\t")
db.reader("TW_Adipose_Subcutaneous_0.5.db","/Users/haoyuzhang/Documents/study/project/Nilanjan/Gallballadar_Cancer/data",coef)
extracted.snp <- read.csv("extract_snps.csv",header = F)
extracted.snp <- t(extracted.snp)
snpvalue <- extracted.snp[5:nrow(extracted.snp),]
snpvalue <- matrix(as.numeric(snpvalue),nrow(snpvalue),3)
snpinfo <- extracted.snp[1:4,]
allele <- snpinfo[3:4,]
position <- as.numeric(snpinfo[2,])
snpid <- gsub(":(.*)","",snpinfo[1,])
snpvaluenew <- matrix(0,nrow(snpvalue)/3,3)
for(i in 1:(nrow(snpvalue)/3)){
  snpvaluenew[i,] <- 0*snpvalue[(3*i-2),]+1*snpvalue[(3*i-1),]+2*snpvalue[(3*i),]
}

extracted.snp.result <- data.frame(pheno$ID,snpvaluenew)
colnames(extracted.snp.result) <- c("individual_id",snpid)
snpinfonew <- data.frame(snpid,position,t(allele))
colnames(snpinfonew) <- c("snpid","position","allele1","allele2")
library(xlsx)
write.xlsx(extracted.snp.result,"extractedsnp.xlsx",sheetName = "snpvalue",row.names = F)
write.xlsx(snpinfonew,"extractedsnp.xlsx",sheetName = "snpinformation",row.names = F,append = T)



pheno <- read.table("pheno_ev",header=T,sep="\t")
fam <- read.table("GBC_final.fam")
bim <- read.table("GBC_final.bim")
save(bim,file="locusmap.rda")
pheno$X <- NULL
pheno$affected_status.1 <- NULL
pheno$affected <- 0
pheno$affected[pheno$affected_status=="Case"] <- 1
pheno$familyID <- fam$V1
idx <- pheno.new$id %in% fam$V2
pheno.new.c <- pheno.new[idx,]
pheno.new.c <- pheno.new.c[match(fam$V2,pheno.new.c$id),]
pheno$highrisk <- pheno.new.c$birth_place_in_high_risk_area
pheno$stoneyn <- pheno.new.c$stone_yn
write.table(pheno,"pheno.txt",quote=F,row.names=F,sep="\t")
##  test correlation between afftected status and principle components
##  found EV3+EV4+EV5+EV8+EV9 were significant correlated with affected status
EV.pheno <- pheno[,17:26]
pheno.EV.cor <- apply(EV.pheno,2,function(x) {cor.test(pheno$affected,x)})



## the code to run on cluster

rm(list=ls()) 
commanarg <- commandArgs(trailingOnly = TRUE)
i1 <- as.numeric(commanarg[1])


library(ggplot2)
library(CGEN)
library(qqman)
setwd("../data")
##  build snp list for function GxE.scan
snp.list <- list()
snp.list$format <- "tped"
file.name <- paste0("GBC_final_a",letters)
file.name <- c(file.name,paste0("GBC_final_b",letters[1:24]))
snp.list$file <- file.name[i1]
#snp.list$start.vec <- 1
#snp.list$stop.vec <- 10
snp.list$subject.list <- "GBC_final.tfam"
##  build phenotype list for function GxE.scan
f = "pheno.txt"
pheno.list <- list(file=f, header=1, delimiter="\t", id.var=c("familyID","ID"), response.var="affected") 
pheno.list$file.type <- 3
pheno.list$main.vars <- ~age+gender+EV3+EV4+EV5+EV8+EV9
pheno.list$strata.var <- ~gender
##  build op list for funciton GxE.scan
op<- list()
op$out.file <- paste0("GxEout_strat_gender",i1)

result <- GxE.scan(snp.list,pheno.list,op)


##  result analysze
##  use GxE.scan.combine to combine all the 50 files together
out.file <- "../result/nov_2/all_output.txt"
dir <- "../result/nov_2"
GxE.scan.combine(out.file,dir) 
result <- read.table("../result/nov_2/all_output.txt",header=T)
##  rematch the snp order as bim file
##  clean the 0 in chromosome
result2 <- result[match(bim$V2,result$SNP),]
result2$chromosome <- bim$V1
result2$chromosome[result2$chromosome==0] <- NA
result2$location <- bim$V4
result.UML <- result2[!is.na(result2$chromosome)&!is.na(result2$UML.Omnibus.Pvalue),]

##build the UML, CML, EB three data frame for manhattan and qqplot

UML.result <- data.frame(CHR = result.UML$chromosome,BP = result.UML$location, 
                         P = result.UML$UML.Omnibus.Pvalue,SNP =result.UML$SNP )

manhattan(UML.result,main="UML.result",chr="CHR",bp="BP",p="P",snp="SNP")
qq(UML.result$P,main="UML.result")
result.CML <- result2[!is.na(result2$chromosome)&!is.na(result2$CML.Omnibus.Pvalue),]
CML.result <- data.frame(CHR = result.CML$chromosome,BP = result.CML$location, 
                         P = result.CML$CML.Omnibus.Pvalue,SNP =result.CML$SNP )

manhattan(CML.result,main="CML.result",chr="CHR",bp="BP",p="P",snp="SNP")
qq(CML.result$P,main="CML.result")

result.EB <- result2[!is.na(result2$chromosome)&!is.na(result2$EB.Omnibus.Pvalue),]
EB.result <- data.frame(CHR = result.EB$chromosome,BP = result.EB$location, 
                         P = result.EB$CML.Omnibus.Pvalue,SNP =result.EB$SNP )

manhattan(EB.result,main="EB.result",chr="CHR",bp="BP",p="P",snp="SNP")
qq(CML.result$P,main="EB.result")



##  my question
##  1.some snp freq is 0
##  2.the chromosome question
##  3.which variable should be called strata.var andwhich should be 
##called int.var
## Sharayu's question double check


pheno.new <- read.csv("pheno.csv")




covar <- data.frame(FamID=fam$V1,SampleID=fam$V2,age=pheno.new.c$age,
                    gender=as.numeric(pheno.new.c$gender),
                    highriskarea=pheno.new.c$birth_place_in_high_risk_area,
                    stoneyn=pheno.new.c$stone_yn,
                    EV3=pheno$EV3,EV4=pheno$EV4,EV5=pheno$EV5,
                    EV8=pheno$EV8,EV9=pheno$EV9)

write.table(covar,"covar.txt",sep=" ",quote=F,row.names=F,col.names=F)
freq <- read.table("GBC_final.frq",header=T)






out.file <- "../result/nov_4/all_output_strat_highrisk.txt"
dir <- "../result/nov_4"
pattern <- "GxEout_strat_highrisk"
GxE.scan.combine(out.file,dir) 
result <- read.table("../result/nov_4/all_output_strat_highrisk.txt",header=T,fill=T)
##  rematch the snp order as bim file
##  clean the 0 in chromosome
result2 <- result[match(bim$V2,result$SNP),]
result2$chromosome <- bim$V1
result2$chromosome[result2$chromosome==0] <- NA
result2$location <- bim$V4
result.UML <- result2[!is.na(result2$chromosome)&!is.na(result2$UML.Omnibus.Pvalue),]

##build the UML, CML, EB three data frame for manhattan and qqplot

UML.result <- data.frame(CHR = result.UML$chromosome,BP = result.UML$location, 
                         P = result.UML$UML.Omnibus.Pvalue,SNP =result.UML$SNP )
# pdf("UMLGWAS.pdf",width=8,height=6, compress = TRUE)
png("UMLGWAS.png",width=8,height=6, 
    units = "in", res = 600)
manhattan(UML.result,main="UML.result",chr="CHR",bp="BP",p="P",snp="SNP")
dev.off()
png("UMLQQ.png",width=8,height=6,units="in",res=600)
qq(UML.result$P,main="UML.result")
dev.off()
result.CML <- result2[!is.na(result2$chromosome)&!is.na(result2$CML.Omnibus.Pvalue),]
CML.result <- data.frame(CHR = result.CML$chromosome,BP = result.CML$location, 
                         P = result.CML$CML.Omnibus.Pvalue,SNP =result.CML$SNP )
png("CMLGWAS.png",width=8,height=6,units="in",res=600)
manhattan(CML.result,main="CML.result",chr="CHR",bp="BP",p="P",snp="SNP")
dev.off()
png("CMLQQ.png",width=8,height=6,units="in",res=600)
qq(CML.result$P,main="CML.result")
dev.off()
result.EB <- result2[!is.na(result2$chromosome)&!is.na(result2$EB.Omnibus.Pvalue),]
EB.result <- data.frame(CHR = result.EB$chromosome,BP = result.EB$location, 
                        P = result.EB$CML.Omnibus.Pvalue,SNP =result.EB$SNP )
png("EBGWAS.png",width=8,height=6,units="in",res=600)
manhattan(EB.result,main="EB.result",chr="CHR",bp="BP",p="P",snp="SNP")
dev.off()
png("EBQQ.png",width=8,height=6,units="in",res=600)
qq(EB.result$P,main="EB.result")
dev.off()
