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


##  summary data analysis
##  read the data in R and change affected status from case and control 
##to 1 and 0
samplesnp <- read.table("sampleSNPs.txt",header=F)
samplesnp <- as.vector(samplesnp$V1)
total <- as.vector(bim$V2)
result <- samplesnp%in% total
sum(result)
library(ggplot2)
setwd("../data")
pheno <- read.table("pheno_ev",header=T,sep="\t")
fam <- read.table("GBC_final.fam")
bim <- read.table("GBC_final.bim")
# idx <- sample(c(1:nrow(bim)),100)
# bim_sample <- bim[idx,]
# write.table(bim_sample,"bim_sample.txt")
pheno.new <- read.csv("pheno.csv")
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
pheno_GCTA <- pheno[,c(28,1,27)]
pheno_GCTA_test <- pheno_GCTA
pheno_GCTA_test$affected <- rnorm(nrow(pheno_GCTA))
write.table(pheno_GCTA_test,"pheno_GCTA_test.txt",sep=" ",quote=F,row.names = F,col.names = F)
write.table(pheno_GCTA,"pheno_GCTA.txt",sep=" ",quote=F,row.names = F,col.names = F)
covar_GCTA <- data.frame(FamID=fam$V1,SampleID=fam$V2,
                    gender=as.numeric(pheno.new.c$gender))
qcovar_GCTA <- data.frame(FamID=fam$V1,SampleID=fam$V2,
                         EV1=pheno$EV1,EV2=pheno$EV2,
                         EV3=pheno$EV3,EV4=pheno$EV4,
                         EV5=pheno$EV5,EV6=pheno$EV6,
                         EV7=pheno$EV7,EV8=pheno$EV8,
                         EV9=pheno$EV9,EV10=pheno$EV10)
write.table(covar_GCTA,"covar_GCTA.txt",sep=" ",quote=F,row.names=F,col.names=F)
write.table(qcovar_GCTA,"qcovar_GCTA.txt",sep=" ",quote=F,row.names = F,col.names = F)
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
png("../result/nov_4/UMLGWAS.png",width=8,height=6, 
    units = "in", res = 600)
manhattan(UML.result,main="UML.result",chr="CHR",bp="BP",p="P",snp="SNP")
dev.off()
png("../result/nov_4/UMLQQ.png",width=8,height=6,units="in",res=600)
qq(UML.result$P,main="UML.result")
dev.off()
result.CML <- result2[!is.na(result2$chromosome)&!is.na(result2$CML.Omnibus.Pvalue),]
CML.result <- data.frame(CHR = result.CML$chromosome,BP = result.CML$location, 
                         P = result.CML$CML.Omnibus.Pvalue,SNP =result.CML$SNP )
png("../result/nov_4/CMLGWAS.png",width=8,height=6,units="in",res=600)
manhattan(CML.result,main="CML.result",chr="CHR",bp="BP",p="P",snp="SNP")
dev.off()
png("../result/nov_4/CMLQQ.png",width=8,height=6,units="in",res=600)
qq(CML.result$P,main="CML.result")
dev.off()
result.EB <- result2[!is.na(result2$chromosome)&!is.na(result2$EB.Omnibus.Pvalue),]
EB.result <- data.frame(CHR = result.EB$chromosome,BP = result.EB$location, 
                        P = result.EB$CML.Omnibus.Pvalue,SNP =result.EB$SNP )
png("../result/nov_4/EBGWAS.png",width=8,height=6,units="in",res=600)
manhattan(EB.result,main="EB.result",chr="CHR",bp="BP",p="P",snp="SNP")
dev.off()
png("../result/nov_4/EBQQ.png",width=8,height=6,units="in",res=600)
qq(EB.result$P,main="EB.result")
dev.off()






##
out.file <- "../result/nov_9/all_output_strat_gender.txt"
dir <- "../result/nov_9"
pattern <- "GxEout_strat_gender"
GxE.scan.combine(out.file,dir) 
result <- read.table("../result/nov_9/all_output_strat_gender.txt",header=T,fill=T)
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
png("../result/nov_9/UMLGWAS.png",width=8,height=6, 
    units = "in", res = 600)
manhattan(UML.result,main="UML.result",chr="CHR",bp="BP",p="P",snp="SNP")
dev.off()
png("../result/nov_9/UMLQQ.png",width=8,height=6,units="in",res=600)
qq(UML.result$P,main="UML.result")
dev.off()
result.CML <- result2[!is.na(result2$chromosome)&!is.na(result2$CML.Omnibus.Pvalue),]
CML.result <- data.frame(CHR = result.CML$chromosome,BP = result.CML$location, 
                         P = result.CML$CML.Omnibus.Pvalue,SNP =result.CML$SNP )
png("../result/nov_9/CMLGWAS.png",width=8,height=6,units="in",res=600)
manhattan(CML.result,main="CML.result",chr="CHR",bp="BP",p="P",snp="SNP")
dev.off()
png("../result/nov_9/CMLQQ.png",width=8,height=6,units="in",res=600)
qq(CML.result$P,main="CML.result")
dev.off()
result.EB <- result2[!is.na(result2$chromosome)&!is.na(result2$EB.Omnibus.Pvalue),]
EB.result <- data.frame(CHR = result.EB$chromosome,BP = result.EB$location, 
                        P = result.EB$CML.Omnibus.Pvalue,SNP =result.EB$SNP )
png("../result/nov_9/EBGWAS.png",width=8,height=6,units="in",res=600)
manhattan(EB.result,main="EB.result",chr="CHR",bp="BP",p="P",snp="SNP")
dev.off()
png("../result/nov_9/EBQQ.png",width=8,height=6,units="in",res=600)
qq(EB.result$P,main="EB.result")
dev.off()



condition <- read.csv("../conditional analysis/condition/Sheet1_Table_1.csv")
condition.add <- condition[seq(1,nrow(condition),2),]
condition.con <- condition[seq(2,nrow(condition),2),]
write.table(condition.add$SNP,"snplist.txt",quote=F,col.names=F,row.names=F)






##conditional analysis
bim <- read.table("GBC_final.bim")
colnames(bim) <- c("chrom","snp","bp","ds","alle1","alle2")

freq <- read.table("GBC_freq.frq",header=T)


library(dplyr)
topsnplist <- read.table("topsnplist",header=F)
colnames(topsnplist) <- "snp"
bim.top <- merge(topsnplist,bim,by.x="snp",sort=F)
freq.top <- merge(topsnplist,freq,by.x="snp",by.y="SNP",sort=F)
snpdata <- read.table("topsnp.raw",header=T)

pheno_new <- read.table("pheno_SNPtest.txt",header=T)

covar <- pheno_new[2:11]
covar_GCTA <- cbind(covar_GCTA,pheno_new$agegroup_18_29,pheno_new$agegroup_30_39,pheno_new$agegroup_40_49,pheno_new$agegroup_60)
write.table(covar_GCTA,"covar_GCTA.txt",quote = F,row.names = F,col.names = F)

idx_female = which(fam$V5==2)
covar_GCTA_female = covar_GCTA[idx_female,]
covar_GCTA_female = covar_GCTA_female[,-3]
qcovar_GCTA_female = qcovar_GCTA[idx_female,]
list_female = fam[idx_female,1:2]
pheno_GCTA_female = pheno_GCTA[idx_female,]

write.table(covar_GCTA_female,file="/Users/haoyuzhang/Dropbox/project/Nilanjan/Gallballadar_Cancer/code/heritability_analysis/covar_GCTA_female.txt",col.names = F,quote=F,row.names = F)
write.table(qcovar_GCTA_female,file="/Users/haoyuzhang/Dropbox/project/Nilanjan/Gallballadar_Cancer/code/heritability_analysis/qcovar_GCTA_female.txt",col.names = F,quote=F,row.names = F)
write.table(list_female,file="/Users/haoyuzhang/Dropbox/project/Nilanjan/Gallballadar_Cancer/code/heritability_analysis/list_female.txt",col.names = F,quote=F,row.names = F)
write.table(pheno_GCTA_female,file="/Users/haoyuzhang/Dropbox/project/Nilanjan/Gallballadar_Cancer/code/heritability_analysis/pheno_GCTA_female.txt",col.names = F,quote=F,row.names = F)





topsnp <- snpdata[,11]
conditional_snp <- snpdata[,27]
topsnp.new <- topsnp
topsnp.new[topsnp==0] <-2
topsnp.new[topsnp==2] <- 0
topsnp <- topsnp.new
x <- conditional_snp
x.new <- x
x.new[x==2] <- 0
x.new[x==0] <- 2
x <- x.new
conditional_snp <- x
model <- glm(pheno_new$affected_status_Case~covar$agegroup_18_29+covar$agegroup_30_39+
               covar$agegroup_40_49+covar$agegroup_60+covar$gender_F+covar$EV3+
               covar$EV4+covar$EV5+covar$EV8+covar$EV9+topsnp,family=binomial)
result <- summary(model)
places=2
log.or <- result$coefficients[2,1]
sd <- result$coefficients[2,2]
log.or.low <- -1.96*sd+log.or
log.or.high <- 1.96*sd+log.or
p <- result$coefficients[2,4]
log.or.top <- result$coefficients[13,1]
sd.top <- result$coefficients[13,2]
log.or.low.top <- -1.96*sd.top+log.or.top
log.or.high.top <- 1.96*sd.top+log.or.top
p.top <- result$coefficients[13,4]
r2 <- round(cor(x,topsnp,use="complete.obs")^2,3)
temp <- cbind(x,topsnp,pheno_new$affected_status_Case)
temp <- temp[complete.cases(temp),]
control.case <- table(temp[,3])
or <- round(exp(log.or),places)
or.low <- round(exp(log.or.low),places)
or.high <- round(exp(log.or.high),places)
or.top <- round(exp(log.or.top),places)
or.low.top <- round(exp(log.or.low.top),places)
or.high.top <- round(exp(log.or.high.top),places) 
cov.log.or <- vcov(model)[c(2,13),c(2,13)]

replication <- read.csv("replication_data.csv",na.strings = "")
replication <- replication[,c(1,12,14,26:31)]
###exclude age with missing value
replication <- replication[c(-322,-325,-344),]
replication <- replication[complete.cases(replication),]

topsnp.r <- replication[,2]
topsnp.r.new <- rep(0,length(topsnp.r))
topsnp.r.new[topsnp.r=="CT"] <- 1
topsnp.r.new[topsnp.r=="TT"] <- 2
topsnp.r <- topsnp.r.new

conditionalsnp.r <- replication[,3]
conditionalsnp.r.new <- rep(0,length(topsnp.r))
conditionalsnp.r.new[conditionalsnp.r=="CT"] <- 1
conditionalsnp.r.new[conditionalsnp.r=="TT"] <- 2
conditionalsnp.r <- conditionalsnp.r.new
age.r <- replication[,8]
age.r.category <- as.factor(cut(age.r,c(17,29,39,49,60,100),labels = F))
age.group <- model.matrix(~age.r.category-1)
age.group <- age.group[,-5]
model.r <- glm(replication$gb~conditionalsnp.r+age.group+replication$sex+topsnp.r,family = binomial)
result.r <- summary(model.r)
places=2
log.or.r <- result.r$coefficients[2,1]
sd.r <- result.r$coefficients[2,2]
log.or.low.r <- -1.96*sd.r+log.or.r
log.or.high.r <- 1.96*sd.r+log.or.r
p.r <- result.r$coefficients[2,4]
log.or.top.r <- result.r$coefficients[8,1]
sd.top.r <- result.r$coefficients[8,2]
log.or.low.top.r <- -1.96*sd.top.r+log.or.top.r
log.or.high.top.r <- 1.96*sd.top.r+log.or.top.r
p.top.r <- result.r$coefficients[8,4]


or.r <- round(exp(log.or.r),places)
or.low.r <- round(exp(log.or.low.r),places)
or.high.r <- round(exp(log.or.high.r),places)
or.top.r <- round(exp(log.or.top.r),places)
or.low.top.r <- round(exp(log.or.low.top.r),places)
or.high.top.r <- round(exp(log.or.high.top.r),places) 
cov.log.or.r <- vcov(model.r)[c(2,8),c(2,8)]

beta.dis <- c(log.or,log.or.top)
beta.r <- c(log.or.r,log.or.top.r)
W.dis <- solve(cov.log.or)
W.r <- solve(cov.log.or.r)
beta.meta <- solve((W.dis+W.r),(W.dis%*%beta.dis+W.r%*%beta.r))
cov.meta <- solve(W.dis+W.r)
log.or.meta <- beta.meta[1,1]
sd.or.meta <- sqrt(cov.meta[1,1])
p.meta <- 2*(1-pnorm(log.or.meta/sd.or.meta))
log.or.low.meta <- log.or.meta-1.96*sd.or.meta
log.or.high.meta <- log.or.meta+1.96*sd.or.meta
or.meta <- exp(log.or.meta)
or.high.meta <- exp(log.or.high.meta)
or.low.meta <- exp(log.or.low.meta)
log.or.top.meta <- beta.meta[2,1]
sd.or.top.meta <- sqrt(cov.meta[2,2])
p.top.meta <- 2*(1-pnorm(log.or.top.meta/sd.or.top.meta))
or.top.meta <- exp(log.or.top.meta)
log.or.low.top.meta <- log.or.top.meta-1.96*sd.or.top.meta
log.or.high.top.meta <- log.or.top.meta+1.96*sd.or.top.meta
or.top.meta <- exp(log.or.top.meta)
or.high.top.meta <- exp(log.or.high.top.meta)
or.low.top.meta <- exp(log.or.low.top.meta)


final.table <- data.frame(conditionalsnp.or=c(or,or.r,or.meta),
                          conditionalsnp.or.low = c(or.low,or.low.r,or.low.meta),
                          conditionalsnp.or.high = c(or.high,or.high.r,or.high.meta),
                          conditionalsnp.p = c(p,p.r,p.meta),
                          topsnp.or = c(or.top,or.top.r,or.top.meta),
                          topsnp.or.low = c(or.low.top,or.low.top.r,or.low.top.meta),
                          topsnp.or.high = c(or.high.top,or.high.top.r,or.high.top.meta),
                          topsnp.p = c(p.top,p.top.r,p.top.meta)
                          )
rownames(final.table) <- c("discovery","replicates","meta")







topsnp.new <- topsnp
topsnp.new[topsnp==0] <-2
topsnp.new[topsnp==2] <- 0
topsnp <- topsnp.new
snpdata <- snpdata[,-12]
snpdata <- snpdata[,-11]

conditionsnp <- function(x){
# x.new <- x
# x.new[x==2] <- 0
# x.new[x==0] <- 2
# x <- x.new
model <- glm(pheno_new$affected_status_Case~x+covar$agegroup_18_29+covar$agegroup_30_39+
               covar$agegroup_40_49+covar$agegroup_60+covar$gender_F+covar$EV3+
               covar$EV4+covar$EV5+covar$EV8+covar$EV9+topsnp,family=binomial)
result <- summary(model)
places=2
log.or <- result$coefficients[2,1]
sd <- result$coefficients[2,2]
log.or.low <- -1.96*sd+log.or
log.or.high <- 1.96*sd+log.or
p <- result$coefficients[2,4]
log.or.top <- result$coefficients[13,1]
sd.top <- result$coefficients[13,2]
log.or.low.top <- -1.96*sd.top+log.or.top
log.or.high.top <- 1.96*sd.top+log.or.top
p.top <- result$coefficients[13,4]
r2 <- round(cor(x,topsnp,use="complete.obs")^2,3)
temp <- cbind(x,topsnp,pheno_new$affected_status_Case)
temp <- temp[complete.cases(temp),]
control.case <- table(temp[,3])
or <- round(exp(log.or),places)
or.low <- round(exp(log.or.low),places)
or.high <- round(exp(log.or.high),places)
or.top <- round(exp(log.or.top),places)
or.low.top <- round(exp(log.or.low.top),places)
or.high.top <- round(exp(log.or.high.top),places) 
return(
  c(control.case[1],control.case[2],r2,
    paste0(or,"(",or.low,"-",
    or.high,")"),
    p,
    paste0(or.top,"(",or.low.top,
             "-",or.high.top,")"),
    p.top
        )
       )
}
condition.result <- apply(snpdata[,7:25],2,conditionsnp)


names <- colnames(condition.result)
names <- gsub("_.*","",names)
condition.result <- t(condition.result)
colnames(condition.result) <- c("control","cases","R2","OR","P-value","P-value for top","OR for top")
condition.result <- as.data.frame(condition.result)
condition.result$snp <- names
condition.result <- merge(topsnplist,condition.result,by="snp",all.x = T,all.y = T
                          ,sort = F)
write.csv(condition.result,"condition.result.csv",col.names =T,quote=F,row.names=F)






#### 2.df association test
library(lmtest)
topsnplist <- read.table("topsnplist",header=F)
colnames(topsnplist) <- "snp"
snpdata <- read.table("topsnp.raw",header=T)
pheno_new <- read.table("pheno_SNPtest.txt",header=T)
replicates <- read.csv("replication_data2.csv")
covar <- pheno_new[2:11]

x <- snpdata[,27]
dftest2 <- function(x){
 
  homo <- rep(0,length(x))
  homo[x==0] <- 1
  heter <- rep(0,length(x))
  heter[x==1] <- 1
  model2 <- glm(pheno_new$affected_status_Case~homo+heter+covar$agegroup_18_29 +covar$agegroup_30_39+
                  covar$agegroup_40_49+covar$agegroup_60+covar$gender_F+covar$EV3+
                  covar$EV4+covar$EV5+covar$EV8+covar$EV9,family=binomial)
  result2 <- summary(model2)
  model3 <- glm(pheno_new$affected_status_Case~covar$agegroup_18_29 +covar$agegroup_30_39+
                  covar$agegroup_40_49+covar$agegroup_60+covar$gender_F+covar$EV3+
                  covar$EV4+covar$EV5+covar$EV8+covar$EV9,family=binomial)
  result3 <- lrtest(model2,model3)
  log.or.homo <- result2$coefficients[2,1]
  sd.or.homo <- result2$coefficients[2,2]
  log.or.homo.low <- log.or.homo -1.96*sd.or.homo
  log.or.homo.high <- log.or.homo +1.96*sd.or.homo
 
  log.or.heter <- result2$coefficients[3,1]
  sd.or.heter <- result2$coefficients[3,2]
  log.or.heter.low <- log.or.heter -1.96*sd.or.heter
  log.or.heter.high <- log.or.heter +1.96*sd.or.heter
  p <- format(result3$`Pr(>Chisq)`[2],digits=3,scientific = T)
  places <- 2
  or.homo <- round(exp(log.or.homo),places)
  or.homo.low <- round(exp(log.or.homo.low),places)
  or.homo.high <- round(exp(log.or.homo.high),places)
  or.heter <- round(exp(log.or.heter),places)
  or.heter.low <- round(exp(log.or.heter.low),places)
  or.heter.high <- round(exp(log.or.heter.high),places)
  return(c(
    paste0(or.homo,"(",or.homo.low,"-",or.homo.high,")"),
    paste0(or.heter,"(",or.heter.low,"-",or.heter.high,")"),
    p
  )
  )
 
 
  }
 


association.test.2df <- apply(snpdata[,7:27],2,dftest2)


names <- colnames(association.test.2df)
names <- gsub("_.*","",names)
association.test.2df <- t(association.test.2df)
colnames(association.test.2df) <- c("or.homo","or.heter","p")
association.test.2df <- as.data.frame(association.test.2df)
association.test.2df$snp <- names
association.test.2df <- merge(topsnplist,association.test.2df,by="snp",all.x = T,all.y = T,sort = F)

write.csv(association.test.2df,"association.test.2df.csv",col.names =T,quote=F,row.names=F)

###dominant analysis
library(lmtest)
topsnplist <- read.table("topsnplist",header=F)
colnames(topsnplist) <- "snp"
snpdata <- read.table("topsnp.raw",header=T)
pheno_new <- read.table("pheno_SNPtest.txt",header=T)

covar <- pheno_new[2:11]

x <- snpdata[,27]
dominant <- function(x){
  x.new <- x
  x.new[x==0] <- 1
  x.new[x==1] <- 1
  x.new[x==2] <- 0
  x <- x.new
  model2 <- glm(pheno_new$affected_status_Case~x+covar$agegroup_18_29 +covar$agegroup_30_39+
                  covar$agegroup_40_49+covar$agegroup_60+covar$gender_F+covar$EV3+
                  covar$EV4+covar$EV5+covar$EV8+covar$EV9,family=binomial)
  result2 <- summary(model2)
  
  log.or <- result2$coefficients[2,1]
  sd.or<- result2$coefficients[2,2]
  log.or.low <- log.or -1.96*sd.or
  log.or.high <- log.or +1.96*sd.or
  p <- result2$coefficients[2,4]
  
  
  places <- 2
  or <- round(exp(log.or),places)
  or.low <- round(exp(log.or.low),places)
  or.high <- round(exp(log.or.high),places)
  
  return(c(
    paste0(or,"(",or.low,"-",or.high,")"),
    p
  )
  )
  
  
}




association.test.2df <- apply(snpdata[,7:27],2,dftest2)


names <- colnames(association.test.2df)
names <- gsub("_.*","",names)
association.test.2df <- t(association.test.2df)
colnames(association.test.2df) <- c("or.homo","or.heter","p")
association.test.2df <- as.data.frame(association.test.2df)
association.test.2df$snp <- names
association.test.2df <- merge(topsnplist,association.test.2df,by="snp",all.x = T,all.y = T,sort = F)

write.csv(association.test.2df,"association.test.2df.csv",col.names =T,quote=F,row.names=F)

###heter and homo/heter
library(lmtest)
topsnplist <- read.table("topsnplist",header=F)
colnames(topsnplist) <- "snp"
snpdata <- read.table("topsnp.raw",header=T)
pheno_new <- read.table("pheno_SNPtest.txt",header=T)

covar <- pheno_new[2:11]

x <- snpdata[,27]
homo.heter <- function(x){
  
  homo <- rep(0,length(x))
  homo[x==0] <- 1
  heter <- rep(0,length(x))
  heter[x==1] <- 1
  model2 <- glm(pheno_new$affected_status_Case~homo+heter+covar$agegroup_18_29 +covar$agegroup_30_39+
                  covar$agegroup_40_49+covar$agegroup_60+covar$gender_F+covar$EV3+
                  covar$EV4+covar$EV5+covar$EV8+covar$EV9,family=binomial)
  result2 <- summary(model2)
  log.or.homo.cov <- result2$cov.unscaled[2,2]
  log.or.heter.cov <- result2$cov.unscaled[3,3]
  log.or.homo.heter.cov <- result2$cov.unscaled[2,3]
  log.or.homoheter.cov <- log.or.homo.cov+log.or.heter.cov-2*log.or.homo.heter.cov
  log.or.homo <- result2$coefficients[2,1]
  sd.or.homo <- result2$coefficients[2,2]
  log.or.homo.low <- log.or.homo -1.96*sd.or.homo
  log.or.homo.high <- log.or.homo +1.96*sd.or.homo
  
  log.or.heter <- result2$coefficients[3,1]
  sd.or.heter <- result2$coefficients[3,2]
  log.or.heter.low <- log.or.heter -1.96*sd.or.heter
  log.or.heter.high <- log.or.heter +1.96*sd.or.heter
 log.or.homoheter <- log.or.homo-log.or.heter
 log.or.homoheter.sd <- sqrt(log.or.homoheter.cov)
 log.or.homoheter.low <- log.or.homoheter -1.96*log.or.homoheter.sd
 log.or.homoheter.high <- log.or.homoheter+1.96*log.or.homoheter.sd
  places <- 2
  or.homoheter <- round(exp(log.or.homoheter),places)
  or.homoheter.low <- round(exp(log.or.homoheter.low),places)
  or.homoheter.high <- round(exp(log.or.homoheter.high),places)
  or.heter <- round(exp(log.or.heter),places)
  or.heter.low <- round(exp(log.or.heter.low),places)
  or.heter.high <- round(exp(log.or.heter.high),places)
  return(c(
    paste0(or.homoheter,"(",or.homoheter.low,"-",or.homoheter.high,")"),
    paste0(or.heter,"(",or.heter.low,"-",or.heter.high,")")
  )
  )
  
  
}



homo.heter.result <- apply(snpdata[,7:27],2,homo.heter)


names <- colnames(homo.heter.result)
names <- gsub("_.*","",names)
homo.heter.result <- t(homo.heter.result)
colnames(homo.heter.result) <- c("or.homo/heter","or.heter")
homo.heter.result <- as.data.frame(homo.heter.result)
homo.heter.result$snp <- names
homo.heter.result <- merge(topsnplist,homo.heter.result,by="snp",all.x = T,all.y = T,sort = F)

write.csv(homo.heter.result,"homo.heter.result.csv",col.names =T,quote=F,row.names=F)



####marginal analysis
library(lmtest)
topsnplist <- read.table("topsnplist",header=F)
colnames(topsnplist) <- "snp"
snpdata <- read.table("topsnp.raw",header=T)
pheno_new <- read.table("pheno_SNPtest.txt",header=T)

covar <- pheno_new[2:11]
x <- snpdata[,27]
marginal <- function(x){
  
#x[x==2] <- 0
#x[x==0] <- 2
  model2 <- glm(pheno_new$affected_status_Case~x+covar$agegroup_18_29 +covar$agegroup_30_39+
                  covar$agegroup_40_49+covar$agegroup_60+covar$gender_F+covar$EV3+
                  covar$EV4+covar$EV5+covar$EV8+covar$EV9,family=binomial)
  result2 <- summary(model2)
 
  log.or <- result2$coefficients[2,1]
  sd.or<- result2$coefficients[2,2]
  log.or.low <- log.or -1.96*sd.or
  log.or.high <- log.or +1.96*sd.or
  p <- result2$coefficients[2,4]

 
  places <- 2
  or <- round(exp(log.or),places)
  or.low <- round(exp(log.or.low),places)
  or.high <- round(exp(log.or.high),places)
 
  return(c(
    paste0(or,"(",or.low,"-",or.high,")"),
    p
  )
  )
  
  
}



marginal <- apply(snpdata[,7:27],2,marginal)


names <- colnames(marginal)
names <- gsub("_.*","",names)
marginal <- t(marginal)
colnames(marginal) <- c("or","p")
marginal <- as.data.frame(marginal)
marginal$snp <- names
marginal <- merge(topsnplist,marginal,by="snp",all.x = T,all.y = T,sort = F)


####test snp

snpdata <- read.table("topsnp.raw",header=T)
idx <- which(pheno$affected==1)
pheno.case <- pheno[idx,]
snp.case <- snpdata[idx,]

gallstonecase <- function(pheno,snp){
  model1 <- glm(pheno$stoneyn~snp+pheno$age+pheno$gender+pheno$EV3+
                  pheno$EV4+pheno$EV5+
                  pheno$EV8+pheno$EV9,family=binomial())
  result <- summary(model1)
  estimate <- result$coefficients[2,1]
  std <- result$coefficients[2,2]
  odds <- exp(estimate)
  CI_lowerbound <- exp(estimate-1.96*std)
  CI_upperbound <- exp(estimate+1.96*std)
  p <- result$coefficients[2,4]
  return(list(P = p,oddsratio = odds, CI_lowerbound = CI_lowerbound,
                    CI_upperbound = CI_upperbound))
}

caseonly <- data.frame(snp=colnames(snpdata)[7:26],
                       P=apply(snp.case[7:26], 2, 
                               function(x){gallstonecase(pheno.case,x)[[1]]}),
                       oddsratio=apply(snp.case[7:26], 
                                       2, function(x){gallstonecase(pheno.case,x)[[2]]}),
                       CI_lowerbound=apply(snp.case[7:26], 2, 
                                           function(x){gallstonecase(pheno.case,x)[[3]]}),
                       CI_upperbound = apply(snp.case[7:26], 2, 
                                             function(x){gallstonecase(pheno.case,x)[[4]]}))

topsnp.order <- read.table("topsnporder.txt",header=F)
write.csv(caseonly,"caseonly.csv",quote=F,row.names=F,col.names=T)

risklogistic<- function(pheno,snp,risk){
  idx <- which(pheno$highrisk==risk)
  testpheno <- pheno[idx,]
  model1 <- glm(pheno[idx,]$affected~snp[idx]+pheno[idx,]$age+
                  pheno[idx,]$gender+
                  pheno[idx,]$EV3+pheno[idx,]$EV4+pheno[idx,]$EV5+
                  pheno[idx,]$EV8+pheno[idx,]$EV9,family = binomial())
return(-log10(summary(model1)$coefficients[2,4]))
}


highrisk <- data.frame(snp=colnames(snpdata)[7:26],
                       high=apply(snpdata[7:26], 2, function(x){risklogistic(pheno,x,1)}),
                       low=apply(snpdata[7:26], 2, function(x){risklogistic(pheno,x,0)}))

colnames(highrisk) <- c("snp","-log10(P) high risk area","-log10(P) low risk area")


write.csv(highrisk,"../result/highriskarea_
          _analysis.csv",quote=F,row.names = F)



genderlogistic<- function(pheno,snp,gender){
  idx <- which(pheno$gender==gender)
  testpheno <- pheno[idx,]
  model1 <- glm(pheno[idx,]$affected~snp[idx]+pheno[idx,]$age+
                  pheno[idx,]$EV3+pheno[idx,]$EV4+pheno[idx,]$EV5+
                  pheno[idx,]$EV8+pheno[idx,]$EV9,family = binomial())
  return((summary(model1)$coefficients[2,4]))
}

genderinteract <- function(pheno,snptry){
  model <- glm(pheno$affected~snptry+pheno$age+pheno$EV3+pheno$EV4+pheno$EV5+
                 pheno$EV8+pheno$EV9+pheno$gender+pheno$gender*snptry+pheno$gender*pheno$EV3+pheno$gender*pheno$EV4+pheno$gender*pheno$EV5+pheno$gender*pheno$EV8+pheno$gender*pheno$EV9+pheno$age*pheno$gender,family=binomial())
  result <- summary(model)
  return(result$coefficients[10,4])
}

genderinteract.result <- apply(snpdata[,7:27],2,function(x){genderinteract(pheno,x)})

pheno$highrisk <- as.numeric(pheno$highrisk)
pheno$highrisk[pheno$highrisk==3] <- NA
pheno$highrisk <- pheno$highrisk-1
highriskinteract <- function(pheno,snptry){
  model <- glm(pheno$affected~snptry+pheno$age+pheno$EV3+pheno$EV4+pheno$EV5+
                 pheno$EV8+pheno$EV9+pheno$highrisk+pheno$highrisk*snptry+pheno$highrisk*pheno$EV3+pheno$highrisk*pheno$EV4+pheno$highrisk*pheno$EV5+pheno$highrisk*pheno$EV8+pheno$highrisk*pheno$EV9+pheno$age*pheno$highrisk+pheno$gender+pheno$gender*pheno$highrisk,family=binomial())
  result <- summary(model)
  return(result$coefficients[11,4])
}
highriskinteract.result <- apply(snpdata[,7:27],2,function(x){highriskinteract(pheno,x)})

idx <- which(pheno$highrisk ==1)
model <- glm(pheno$affected[idx]~pheno$gender[idx]+snptry[idx]+pheno$age[idx]+pheno$EV3[idx]+pheno$EV4[idx]+
               pheno$EV5[idx]+pheno$EV8[idx]+pheno$EV9[idx],family = binomial())
pheno$stoneyn <- as.numeric(pheno$stoneyn)
pheno$stoneyn[pheno$stoneyn==3] <- NA
pheno$stoneyn <- pheno$stoneyn-1

stoneinteract <- function(pheno,snptry){
  model <- glm(pheno$affected~snptry+pheno$age+pheno$EV3+pheno$EV4+pheno$EV5+
                 pheno$EV8+pheno$EV9+pheno$stoneyn+pheno$stoneyn*snptry+pheno$stoneyn*pheno$EV3+pheno$stoneyn*pheno$EV4+pheno$stoneyn*pheno$EV5+pheno$stoneyn*pheno$EV8+pheno$stoneyn*pheno$EV9+pheno$age*pheno$stoneyn+pheno$gender+pheno$gender*pheno$stoneyn,family=binomial())
  result <- summary(model)
  return(result$coefficients[11,4])
}
stoneinteract.result <- apply(snpdata[,7:27],2,function(x){stoneinteract(pheno,x)})

genderstrat <- data.frame(snp=colnames(snpdata)[7:26],
                       Female=apply(snpdata[7:26], 2, function(x){genderlogistic(pheno,x,"F")}),
                       Male=apply(snpdata[7:26], 2, function(x){genderlogistic(pheno,x,"M")}))
colnames(genderstrat) <- c("snp","-log10(p) Female","-log10(P) Male")


write.csv(genderstrat,"../result/gender_
          _analysis.csv",quote=F,row.names = F)



gallstonelogistic<- function(pheno,snp,stoneyn){
  idx <- which(pheno$stoneyn==stoneyn)
  testpheno <- pheno[idx,]
  model1 <- glm(pheno[idx,]$affected~snp[idx]+pheno[idx,]$age+
                  pheno[idx,]$gender+
                  pheno[idx,]$EV3+pheno[idx,]$EV4+pheno[idx,]$EV5+
                  pheno[idx,]$EV8+pheno[idx,]$EV9,family = binomial())
  return(
    c((1/exp(summary(model1)$coefficients[2,1])),
      1/(exp(summary(model1)$coefficients[2,1]+1.96*summary(model1)$coefficients[2,2])),
      1/(exp(summary(model1)$coefficients[2,1]-1.96*summary(model1)$coefficients[2,2])),
      summary(model1)$coefficients[2,4])
    )
}


gallstone <- cbind( apply(snpdata[7:26], 2, function(x){gallstonelogistic(pheno,x,1)}),
                    apply(snpdata[7:26], 2, function(x){gallstonelogistic(pheno,x,0)}))

write.csv(t(gallstone),"gallstone.csv",quote=F,row.names = T,col.names=T)
gallstonestrat <- data.frame(snp=colnames(snpdata)[7:26],
                       gallstone=apply(snpdata[7:26], 2, function(x){gallstonelogistic(pheno,x,1)}),
                       nogallstone=apply(snpdata[7:26], 2, function(x){gallstonelogistic(pheno,x,0)}))

colnames(gallstonestrat) <- c("snp","-log10(p) gallstone","-log10(P) nogallstone")


write.csv(gallstonestrat,"../result/gallstone_
          _analysis.csv",quote=F,row.names = F)




apply(snpdata[7:26], 2, function(x){logistic(pheno,x)})



apply(snpdata[7:26], 2, function(x){risklogistic(pheno,x,0|1)})



model1 <- glm(pheno$affected~snpdata$rs1558376_A+pheno$age+pheno$gender+
                pheno$EV3+pheno$EV4+pheno$EV5+pheno$EV8+pheno$EV9,family = binomial())







#########caseonly study
## the code to run on cluster




library(ggplot2)
library(CGEN)
library(qqman)
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
pheno.list$main.vars <- ~age+gender+EV3+EV4+EV5+EV8+EV9+stoneyn
pheno.list$int.vars <- ~stoneyn
##  build op list for funciton GxE.scan
op<- list()


result <- GxE.scan(snp.list,pheno.list,op)
result <- read.table("GxE.scan.output.txt",header=T)
caseonly_result2 <- result[,c(1,7,10,11,19,20)]
caseonly_result2$oddsratio <- exp(caseonly_result2$CML.SNP.stoneyn.Beta)
caseonly_result2$lower_CI <- exp(caseonly_result2$CML.SNP.stoneyn.Beta-1.96*caseonly_result2$CML.SNP.stoneyn.SE)
caseonly_result2$upper_CI <- exp(caseonly_result2$CML.SNP.stoneyn.Beta+1.96*caseonly_result2$CML.SNP.stoneyn.SE)


write.csv(caseonly_result2,"case-case_analysis.csv",quote=F,row.names = F)



####stratifed analysis forest plot
library(forestplot)
gall_stone <- read.csv("stratified_analysis_gallstone.csv")
OR <- with(gall_stone,cbind(OR,OR.1))

low <- with(gall_stone,cbind(OR_Low,OR_Low.1))

high <- with(gall_stone,cbind(OR_high,OR_high.1))
row_names <- as.list(c(expression(paste("Cluster1")),rep(NA,11),expression(paste("Cluster2")),rep(NA,8)))
row_names <- list(categories=row_names,snp=as.list(as.character(gall_stone$snp)))
png(filename = "forest_plot.png",res=600,units="in",width=6,height = 6)
forestplot(row_names,OR,low,high,new_page = TRUE,boxsize = 0.2,
           xticks = c(0.5,1,1.5,2,2.5,3,3.5,4,4.5,5),
           hrzl_lines = gpar(lwd=15),
           clip=c(0,5),
           zero=c(1,1),
           title="
           Stratified analysis for gall stone",
           col=fpColors(box=c("royalblue","gold"),
                        line=c("darkblue","orange")), 
           legend=c("With gallstone","Without gallstone"),
           xlab="Odds Ratio",
           legend_args=fpLegend(pos=list("topright"),
                                title="Group",
                                r=unit(.1,"snpc"),
                                gp=gpar(col="#CCCCCC",lwd=1.5)
                              ))
dev.off()


forestplot(row_names,OR,low,high,new_page = TRUE,boxsize = 0.2,
           xticks = c(0.5,1,1.5,2,2.5,3,3.5,4,4.5,5),
           hrzl_lines = gpar(lwd=15),
           clip=c(0,5),
           zero=c(1,1),
           title="
           Stratified analysis for gall stone",
           col=fpColors(box=c("royalblue","gold"),
                        line=c("darkblue","orange")), 
           legend=c("With gallstone","Without gallstone"),
           xlab="Odds Ratio",
           legend_args=fpLegend(pos=list("right"),
                                title="Group",
                                r=unit(.1,"snpc"),
                                gp=gpar(col="#CCCCCC",lwd=1.4)
           ))
### 
analysis and conditional analysis based on gallstone
setwd("../data")
pheno <- read.table("pheno_ev",header=T,sep="\t")

pheno.new <- read.csv("pheno.csv")
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

topsnplist <- read.table("topsnplist",header=F)
colnames(topsnplist) <- "snp"
snpdata <- read.table("topsnp.raw",header=T)
pheno_new <- read.table("pheno_SNPtest.txt",header=T)
covar <- pheno_new[,2:11]
covar <- pheno_new[2:11]
topsnp <- snpdata[,12]
topsnp.new <- topsnp
topsnp.new[topsnp==0] <-2
topsnp.new[topsnp==2] <- 0
topsnp <- topsnp.new
snpdata <- snpdata[,-12]
snpdata <- snpdata[,-11]
x <- snpdata[,8]
condition.
<- function(x,i){
  idx <- which(pheno$stoneyn==i)
  x.new <- x
  x.new[x==2] <- 0
  x.new[x==0] <- 2
  x <- x.new
  model <- glm(pheno_new$affected_status_Case[idx]~x[idx]+covar$agegroup_18_29[idx] +covar$agegroup_30_39[idx]+
                 covar$agegroup_40_49[idx]+covar$agegroup_60[idx]+covar$gender_F[idx]+covar$EV3[idx]+
                 covar$EV4[idx]+covar$EV5[idx]+covar$EV8[idx]+covar$EV9[idx]+topsnp[idx],family=binomial)
  result <- summary(model)
  places=2
  log.or <- result$coefficients[2,1]
  sd <- result$coefficients[2,2]
  log.or.low <- -1.96*sd+log.or
  log.or.high <- 1.96*sd+log.or
  p <- result$coefficients[2,4]
  log.or.top <- result$coefficients[13,1]
  sd.top <- result$coefficients[13,2]
  log.or.low.top <- -1.96*sd.top+log.or.top
  log.or.high.top <- 1.96*sd.top+log.or.top
  p.top <- result$coefficients[13,4]
  r2 <- round(cor(x[idx],topsnp[idx],use="complete.obs")^2,3)
  temp <- cbind(x[idx],topsnp[idx],pheno_new$affected_status_Case[idx])
  temp <- temp[complete.cases(temp),]
  control.case <- table(temp[,3])
  or <- round(exp(log.or),places)
  or.low <- round(exp(log.or.low),places)
  or.high <- round(exp(log.or.high),places)
  or.top <- round(exp(log.or.top),places)
  or.low.top <- round(exp(log.or.low.top),places)
  or.high.top <- round(exp(log.or.high.top),places) 
  return(
    c(control.case[1],control.case[2],r2,
      paste0(or,"(",or.low,"-",
             or.high,")"),
      p,
      paste0(or.top,"(",or.low.top,
             "-",or.high.top,")"),
      p.top
    )
  )
}
condition.nostone <- apply(snpdata[,7:25],2,function(x)condition.
                           (x,0))


names <- colnames(condition.nostone)
names <- gsub("_.*","",names)
condition.nostone <- t(condition.nostone)
colnames(condition.nostone) <- c("control","cases","R2","OR","P-value","P-value for top","OR for top")
condition.nostone <- as.data.frame(condition.nostone)
condition.nostone$snp <- names
condition.nostone <- merge(topsnplist,condition.nostone,by="snp",all.x = T,all.y = T
                          ,sort = F)
write.csv(condition.nostone,"condition.nostone.csv",col.names =T,quote=F,row.names=F)


condition.
.for.plot <- function(x,i){
  idx <- which(pheno$stoneyn==i)
  x.new <- x
  x.new[x==2] <- 0
  x.new[x==0] <- 2
  x <- x.new
  model <- glm(pheno_new$affected_status_Case[idx]~x[idx]+covar$agegroup_18_29[idx] +covar$agegroup_30_39[idx]+
                 covar$agegroup_40_49[idx]+covar$agegroup_60[idx]+covar$gender_F[idx]+covar$EV3[idx]+
                 covar$EV4[idx]+covar$EV5[idx]+covar$EV8[idx]+covar$EV9[idx]+topsnp[idx],family=binomial)
  result <- summary(model)
  places=2
  log.or <- result$coefficients[2,1]
  sd <- result$coefficients[2,2]
  log.or.low <- -1.96*sd+log.or
  log.or.high <- 1.96*sd+log.or
  p <- result$coefficients[2,4]
  log.or.top <- result$coefficients[13,1]
  sd.top <- result$coefficients[13,2]
  log.or.low.top <- -1.96*sd.top+log.or.top
  log.or.high.top <- 1.96*sd.top+log.or.top
  p.top <- result$coefficients[13,4]
  r2 <- round(cor(x[idx],topsnp[idx],use="complete.obs")^2,3)
  temp <- cbind(x[idx],topsnp[idx],pheno_new$affected_status_Case[idx])
  temp <- temp[complete.cases(temp),]
  control.case <- table(temp[,3])
  or <- round(exp(log.or),places)
  or.low <- round(exp(log.or.low),places)
  or.high <- round(exp(log.or.high),places)
  or.top <- round(exp(log.or.top),places)
  or.low.top <- round(exp(log.or.low.top),places)
  or.high.top <- round(exp(log.or.high.top),places) 
  return(
   c(or,or.low,or.high)
    )
  
}

condition.nostone <- apply(snpdata[,7:25],2,function(x)condition.
                           .for.plot(x,0))


names <- colnames(condition.nostone)
names <- gsub("_.*","",names)
condition.nostone <- t(condition.nostone)
colnames(condition.nostone) <- c("OR","OR_low","OR_high")
condition.nostone <- as.data.frame(condition.nostone)
condition.nostone$snp <- names
condition.nostone <- merge(topsnplist,condition.nostone,by="snp",all.x = T,all.y = T
                           ,sort = F)
condition.nostone <- condition.nostone[1:19,]
condition.nostone$Set <- c(rep(1,10),rep(2,9))
write.csv(condition.nostone,"condition.nostone.plot.csv",col.names =T,quote=F,row.names=F)
library(forestplot)
gall_stone <- read.csv("condition.nostone.plot.csv")
OR <- gall_stone$OR

low <- gall_stone$OR_low

high <- gall_stone$OR_high
row_names <- as.list(c("Set1",rep(NA,9),"Set2",rep(NA,8)))
row_names <- list(categories=row_names,snp=as.list(as.character(gall_stone$snp)))
png(filename = "forest_plot.nostone.conditional.png",res=600,units="in",width=6,height = 6)
forestplot(row_names,OR,low,high,new_page = TRUE,boxsize = 0.2,
           xticks = c(0.5,1,1.5,2,2.5,3,3.5,4,4.5,5),
           hrzl_lines = gpar(lwd=15),
           clip=c(0,5),
           zero=c(1,1),
           title="conditional analysis reference to rs1558376 in no stone group",
           xlab="Odds Ratio",
           )
dev.off()


data <- read.csv("data_region.csv")
fam <- read.table("GBC_final.fam")
data_all <- merge(data,fam,by.x ="id",by.y ="V2",all=T)

idx_central <- which(data_all$Birth.state.code=="Central")
central <- data_all[idx_central,c(6,1)]
write.table(central,"./pca_freq/central_id.txt",row.names = F,col.names = F,quote = F)

idx_north <- which(data_all$Birth.state.code=="North")
north <- data_all[idx_north,c(6,1)]
write.table(north,"./pca_freq/north_id.txt",row.names = F,col.names = F,quote = F)

idx_north_east <- which(data_all$Birth.state.code=="North East")
north_east <- data_all[idx_north_east,c(6,1)]
write.table(north_east,"./pca_freq/north_east_id.txt",row.names = F,col.names = F,quote = F)

idx_south <- which(data_all$Birth.state.code=="South")
south <- data_all[idx_south,c(6,1)]
write.table(south,"./pca_freq/south_id.txt",row.names = F,col.names = F,quote = F)

idx_west <- which(data_all$Birth.state.code=="West")
west <- data_all[idx_west,c(6,1)]
write.table(west,"./pca_freq/west_id.txt",row.names = F,col.names = F,quote = F)

idx_west_bengal <- which(data_all$Birth.state.code=="West Bengal")
west_bengal <- data_all[idx_west_bengal,c(6,1)]
write.table(west_bengal,"./pca_freq/west_bengal_id.txt",row.names = F,col.names = F,quote = F)