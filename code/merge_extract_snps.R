#Goal: merge extracted SNPs and statistical analysis
setwd("/dcl01/chatterj/data/GB/NEW/KG/snpR")


# filesDir <- '/dcl01/chatterj/data/GB/NEW/KG/snpR/extract_snps/'
# files <- dir(filesDir,pattern="result_extract",full.names=T)
# result <- NULL
# result.all <- NULL
# for(i in 1:length(files)){
#   print(i)
#   load(files[i])
#   result.all <- rbind(result.all,result)
# }
# setwd('/users/hzhang1/R/GBV/result')
# save(result.all,file = "./extract_snp_all.Rdata")
SNPwiseFun <- function(casecontrol,covar,gene,stone){
  model <- glm(casecontrol~covar$agegroup_18_29+covar$agegroup_30_39+
                 covar$agegroup_40_49+covar$agegroup_60+covar$gender+covar$EV3+
                 covar$EV4+covar$EV5+covar$EV8+covar$EV9+gene,family=binomial)
  result <- summary(model)
  temp.result <-round(exp(c(result$coefficient[12,1],
                            result$coefficient[12,1]-1.96*result$coefficient[12,2],
                            result$coefficient[12,1]+1.96*result$coefficient[12,2])),2)
  
  function.result <- list()
  function.result [[1]]<-paste0(temp.result[1]," [",
                                temp.result[2],",",
                                temp.result[3],"]")
  function.result [[2]] <- result$coefficients[12,4]
  function.result[[3]] <- result$coefficient[12,1]
  model2 <- glm(stone~covar$agegroup_18_29+covar$agegroup_30_39+
                  covar$agegroup_40_49+covar$agegroup_60+covar$gender+covar$EV3+
                  covar$EV4+covar$EV5+covar$EV8+covar$EV9+gene,family=binomial)
  result <- summary(model2)
  temp.result <-round(exp(c(result$coefficient[12,1],
                            result$coefficient[12,1]-1.96*result$coefficient[12,2],
                            result$coefficient[12,1]+1.96*result$coefficient[12,2])),2)
  
  #function.result <- list()
  function.result [[4]]<-paste0(temp.result[1]," [",
                                temp.result[2],",",
                                temp.result[3],"]")
  function.result [[5]] <- result$coefficients[12,4]
  idx.control <- which(casecontrol==0)
  freq <- sum(gene[idx.control])/(2*length(idx.control))
  function.result[[6]] <- result$coefficient[12,1]
  function.result[[7]] <- freq
  return(function.result)
}

PRSFun <- function(y,PRS){
  model <- glm(y~PRS,family=binomial)
  result <- summary(model)
  temp.result <-round(exp(c(result$coefficient[2,1],
                            result$coefficient[2,1]-1.96*result$coefficient[2,2],
                            result$coefficient[2,1]+1.96*result$coefficient[2,2])),2)
  function.result <- list()
  function.result[[1]] <- result$coefficient[2,1]
  function.result[[2]] <- result$coefficient[2,2]
  function.result [[3]]<-paste0(temp.result[1]," [",
                                temp.result[2],",",
                                temp.result[3],"]")
  function.result [[4]] <- result$coefficients[2,4]
  return(function.result)
}

setwd('/users/hzhang1/R/GBV/result')
load("./extract_snp_all.Rdata")
snp.data <- read.csv("/users/hzhang1/R/GBV/data/extracted_SNPs_information_new.csv",stringsAsFactors = F)
#three SNPs can't be found in the dataset
snp.data <- snp.data[-c(28:30),]
snp.list <- snp.data$SNP
chr.list <-  snp.data$CHR
pos.list <- snp.data$Position
snp.infor <- data.frame(snp.list,chr.list,pos.list)
library(dplyr)
colnames(snp.infor)[3] <- "pos"
#head(result.all[,1:5])
colnames(result.all)[3] <- "pos"
#match by position
result.new <- left_join(snp.infor,result.all,by="pos")
result.new <- result.new[,c(1,2,3,6,7,8:ncol(result.new))]
#number of subject
n <- (ncol(result.new)-5)/3
#compute imputed genotype
genotype <- matrix(0,n,nrow(result.new))
temp = 6
for(i in 1:n){
  genotype[i,] <- 0*result.new[,temp]+1*result.new[,temp+1]+2*result.new[,temp+2]
  temp = temp+3
}
##Minor.allele reference
MA <- snp.data[,4]
##code the genotype with reference to minor allele
idx <- which(as.character(result.new$V5)!=MA)
genotype[,idx] <- 2-genotype[,idx]
#freq = apply(genotype,2,function(x){sum(x)/(2*length(x))})
#match the order of imputed file to pheno file
sample.order <- read.table("/dcl01/chatterj/data/GB/NEW/sample.lst")
genotype <- data.frame(sample.order,genotype,stringsAsFactors = F)
colnames(genotype)[1] <- "GENO_PID"
head(result.new[,1:6])
library(data.table)
pheno <-as.data.frame(fread("/users/hzhang1/R/GBC/data/pheno.txt",header = T))
n <- nrow(pheno)
agegroup_18_29 <- rep(0,n)
idx <- which(pheno$age>=18&pheno$age<=29)
agegroup_18_29[idx] <- 1
#sum(agegroup_18_29)
agegroup_30_39 <- rep(0,n)
idx <- which(pheno$age>=30&pheno$age<=39)
agegroup_30_39[idx] <- 1
agegroup_40_49 <- rep(0,n)
idx <- which(pheno$age>=40&pheno$age<=49)
agegroup_40_49[idx] <- 1
agegroup_50_59 <- rep(0,n)
idx <- which(pheno$age>=50&pheno$age<=59)
agegroup_50_59[idx] <- 1
agegroup_60 <- rep(0,n)
idx <- which(pheno$age>=60)
agegroup_60[idx] <- 1
pheno_new <- cbind(pheno[,c(2)],
                   agegroup_18_29,
                   agegroup_30_39,
                   agegroup_40_49,
                   agegroup_50_59,
                   agegroup_60,
                   pheno[,c(13,19,20,21,24,25,10,30)])
pheno_order <- pheno_new[,1,drop=F]
colnames(pheno_order) <- "GENO_PID"
genotype_match <- left_join(pheno_order,genotype,by="GENO_PID")
colnames(pheno_new)[1] <- "GENO_PID"
pheno_all <- left_join(pheno_new,genotype,by="GENO_PID")
colnames(pheno_all)[15:43] <- as.character(result.new[,1])

covar <- pheno_all[,2:12]
gene <- pheno_all[,15]
casecontrol <- rep(0,n)
casecontrol[pheno_all$affected_status=="Case"] <- 1
stone <- as.numeric(pheno_all[,14])
n.snp <- ncol(genotype)-1
OR_1 <- rep("c",n.snp)
p1 <- rep(0,n.snp)
logOR_1 <- rep(0,n.snp)
OR_2 <- rep("c",n.snp)
p2 <- rep(0,n.snp)
logOR_2 <- rep(0,n.snp)
freq <- rep(0,n.snp)
for(i in 1:n.snp){
  gene <- pheno_all[,14+i]
  function.result <- SNPwiseFun(casecontrol,covar,gene,stone)
  OR_1[i] <- function.result[[1]]
  p1[i] <- function.result[[2]]
  logOR_1[i] <- function.result[[3]]
  OR_2[i] <- function.result[[4]]
  p2[i] <- function.result[[5]]
  logOR_2[i] <- function.result[[6]]
  freq[i] <- function.result[[7]]
}
model.result <- data.frame(freq,OR_1,p1,logOR_1,OR_2,p2,logOR_2,stringsAsFactors=F)
colnames(model.result) <- c("MAF IND",
                            "OR (95% CI) GallBladderCancer (IND)",
                            "P GallBladderCancer",
                            "log OR GallBladderCancer (IND)",
                            "OR (95% CI) GallStone (IND)",
                            "P GallStone (IND)",
                            "log OR GallStone (IND)")

write.csv(model.result,file = "./Gallbladder_analysis_29SNPs.csv")

#prs
idx <- which(freq>=0.005)
all.gene <- as.matrix(pheno_all[,14+idx])
#PRS1 <- all.gene%*%logOR_1
logOR <- logOR_2[idx]
PRS <- all.gene%*%logOR
idx.control <- which(casecontrol==0)
scale.factor <- sd(PRS[idx.control])
#PRS <- all.gene%*%logOR/sd(PRS[idx.control])
PRS <- all.gene%*%logOR/scale.factor
temp.result <- PRSFun(casecontrol,PRS)
prs.result <- data.frame(temp.result[[1]],
                         temp.result[[2]],
                         temp.result[[3]],
                         temp.result[[4]])
colnames(prs.result) <- c("logORpersd",
                          "sdoflogOR",
                          "ORpersd",
                          "P")


idx.stone <- which(pheno_all$stoneyn==1)
temp.result <- PRSFun(casecontrol[idx.stone],PRS[idx.stone])
prs.result.temp <- data.frame(temp.result[[1]],
                              temp.result[[2]],
                              temp.result[[3]],
                              temp.result[[4]])
colnames(prs.result.temp) <- c("logORpersd",
                               "sdoflogOR",
                               "ORpersd",
                               "P")

prs.result <- rbind(prs.result,
                    prs.result.temp)

idx.nostone <- which(pheno_all$stoneyn==0)
temp.result <- PRSFun(casecontrol[idx.nostone],PRS[idx.nostone])
prs.result.temp <- data.frame(temp.result[[1]],
                              temp.result[[2]],
                              temp.result[[3]],
                              temp.result[[4]])
colnames(prs.result.temp) <- c("logORpersd",
                               "sdoflogOR",
                               "ORpersd",
                               "P")
prs.result <- rbind(prs.result,
                    prs.result.temp)
#use the log OR from the meta-analysis of UKBB and Iceland
log.or.UK <- rep(0,n.snp)
for(i in 1:n.snp){
  log.or.UK[i] <- log(as.numeric(strsplit(snp.data[i,8]," ")[[1]][1]))  
}
logOR <- log.or.UK[idx]
PRS <- all.gene%*%logOR
idx.control <- which(casecontrol==0)
scale.factor.uk <- sd(PRS[idx.control])
#PRS <- all.gene%*%logOR/sd(PRS[idx.control])
PRS <- all.gene%*%logOR/scale.factor.uk
temp.result <- PRSFun(casecontrol,PRS)
prs.result.new <- data.frame(temp.result[[1]],
                         temp.result[[2]],
                         temp.result[[3]],
                         temp.result[[4]])
colnames(prs.result.new) <- c("logORpersd",
                          "sdoflogOR",
                          "ORpersd",
                          "P")
idx.stone <- which(pheno_all$stoneyn==1)
temp.result <- PRSFun(casecontrol[idx.stone],PRS[idx.stone])
prs.result.temp <- data.frame(temp.result[[1]],
                              temp.result[[2]],
                              temp.result[[3]],
                              temp.result[[4]])
colnames(prs.result.temp) <- c("logORpersd",
                               "sdoflogOR",
                               "ORpersd",
                               "P")

prs.result.new <- rbind(prs.result.new,
                    prs.result.temp)

idx.nostone <- which(pheno_all$stoneyn==0)
temp.result <- PRSFun(casecontrol[idx.nostone],PRS[idx.nostone])
prs.result.temp <- data.frame(temp.result[[1]],
                              temp.result[[2]],
                              temp.result[[3]],
                              temp.result[[4]])
colnames(prs.result.temp) <- c("logORpersd",
                               "sdoflogOR",
                               "ORpersd",
                               "P")
prs.result.new <- rbind(prs.result.new,
                    prs.result.temp)

#stratified analysis
#gallstone group
pheno_all_stone <- pheno_all[(pheno_all$stoneyn==1),]
covar <- pheno_all_stone[,2:12]
gene <- pheno_all_stone[,15]
n <- nrow(pheno_all_stone)
casecontrol <- rep(0,n)
casecontrol[pheno_all_stone$affected_status=="Case"] <- 1
stone <- as.numeric(pheno_all_stone[,14])
n.snp <- ncol(genotype)-1
OR_1 <- rep("c",n.snp)
p1 <- rep(0,n.snp)
logOR_1 <- rep(0,n.snp)
freq <- rep(0,n.snp)
for(i in 1:n.snp){
  gene <- pheno_all_stone[,14+i]
  function.result <- SNPwiseFun(casecontrol,covar,gene,stone)
  OR_1[i] <- function.result[[1]]
  p1[i] <- function.result[[2]]
  logOR_1[i] <- function.result[[3]]
  freq[i] <- function.result[[7]]
}
model.result.stone <- data.frame(freq,OR_1,p1,logOR_1,freq,stringsAsFactors=F)
colnames(model.result.stone) <- c("MAF IND",
                            "OR (95% CI) GallBladderCancer (IND)",
                            "P GallBladderCancer",
                            "log OR GallBladderCancer (IND)")
write.csv(model.result.stone,file = "./Gallbladderstone_analysis_29SNPs.csv")
idx <- which(freq>=0.005)


#nogallstone group
pheno_all_stone <- pheno_all[(pheno_all$stoneyn==0),]
covar <- pheno_all_stone[,2:12]
gene <- pheno_all_stone[,15]
n <- nrow(pheno_all_stone)
casecontrol <- rep(0,n)
casecontrol[pheno_all_stone$affected_status=="Case"] <- 1
stone <- as.numeric(pheno_all_stone[,14])
n.snp <- ncol(genotype)-1
OR_1 <- rep("c",n.snp)
p1 <- rep(0,n.snp)
logOR_1 <- rep(0,n.snp)
freq <- rep(0,n.snp)
for(i in 1:n.snp){
  gene <- pheno_all_stone[,14+i]
  function.result <- SNPwiseFun(casecontrol,covar,gene,stone)
  OR_1[i] <- function.result[[1]]
  p1[i] <- function.result[[2]]
  logOR_1[i] <- function.result[[3]]
  freq[i] <- function.result[[7]]
}
model.result.nostone <- data.frame(freq,OR_1,p1,logOR_1,stringsAsFactors=F)
colnames(model.result.nostone) <- c("MAF IND",
                            "OR (95% CI) GallBladderCancer (IND)",
                            "P GallBladderCancer",
                            "log OR GallBladderCancer (IND)")
write.csv(model.result.nostone,file = "./Gallbladdernostone_analysis_29SNPs.csv")
write.csv(pheno_all,file = "pheno_all_GBC.csv")
write.csv(prs.result,file = "PRS_result_GBC.csv")
write.csv(prs.result.new,file = "PRS_result_GBC_UK.csv")



##statistical function



# pheno_SNPtest.txtx
# file <- files[i1]

#load("/dcl01/chatterj/data/GB/NEW/KG/mergesnp.Rdata")
# load(file)
# chr.temp <- gsub("/dcl01/chatterj/data/GB/NEW/KG/snpR/chr","",file)
# chr.temp <- gsub("","",chr.temp)
# chr <- as.numeric(strsplit(chr.temp,"_")[[1]][1])
# chr.vec <- rep(chr,nrow(data))
# pos.vec <- data[,3]
# chr.pos <- paste0(chr.vec,":",pos.vec)
# chr.pos.target <- paste0(chr.list,":",pos.list)
# # snpid <- data$V2
# # snpid <- as.character(snpid)
# # idx <- grep("rs",snpid)
# # data <- data[idx,]
# # snpid <- data$V2
# # snpid <- as.character(snpid)
# # snpid <- gsub(":(.*)","",snpid)
# try <- chr.pos%in%chr.pos.target
# if(sum(try)!=0){
#   result <- data[try,]
#   save(result,file=paste0("/dcl01/chatterj/data/GB/NEW/KG/snpR/extract_snps/result_extract",i1,".rda"))
# }
# print(sum(try))
# print("finished")
