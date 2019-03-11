bim <- read.table("GBC_final.bim")
colnames(bim) <- c("chrom","snp","bp","ds","alle1","alle2")

freq <- read.table("GBC_freq.frq",header=T)


library(dplyr)
topsnplist <- read.table("topsnplist",header=F)
colnames(topsnplist) <- "snp"
bim.top <- merge(topsnplist,bim,by.x="snp",sort=F)
freq.top <- merge(topsnplist,freq,by.x="snp",by.y="SNP",sort=F)
snpdata <- read.table("topsnp.raw",header=T)

pheno_new <- read.table("../data/pheno_SNPtest.txt",header=T)
covar <- pheno_new[,2:11]
load("predexp_all.rda")
gene <- predexp.all[,2]
model <- glm(pheno_new$affected_status_Case~covar$agegroup_18_29+covar$agegroup_30_39+
               covar$agegroup_40_49+covar$agegroup_60+covar$gender_F+covar$EV3+
               covar$EV4+covar$EV5+covar$EV8+covar$EV9+gene,family=binomial)


twas <- function(gene){
  model <- glm(pheno_new$affected_status_Case~covar$agegroup_18_29+covar$agegroup_30_39+
                 covar$agegroup_40_49+covar$agegroup_60+covar$gender_F+covar$EV3+
                 covar$EV4+covar$EV5+covar$EV8+covar$EV9+gene,family=binomial)
  result <- summary(model)
  return(c(result$coefficients[12,1],result$coefficients[12,4]))
  
}

twas_result <- matrix(0,(ncol(predexp.all)-1),2)
for(i in 1:(ncol(predexp.all)-1)){
  print(i)
  geno <- predexp.all[,i+1]
  if(is.na(geno)){
    twas_result[i,] <- c(NA,NA)
  }else{
    twas_result[i,] <- twas(predexp.all[,i+1])
  }
 }


covar <- pheno_new[2:11]


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
replication <- replication[,c(1,12,14,18,26:31)]
###exclude age with missing value
replication <- replication[c(-322,-325,-344),]
replication <- replication[complete.cases(replication),]










#############for snp rs1558375 

topsnp.r <- replication[,2]
topsnp.r.new <- rep(0,length(topsnp.r))
topsnp.r.new[topsnp.r=="CT"] <- 1
topsnp.r.new[topsnp.r=="TT"] <- 2
topsnp.r <- topsnp.r.new


age.r <- replication[,9]
age.r.category <- as.factor(cut(age.r,c(17,29,39,49,60,100),labels = F))
age.group <- model.matrix(~age.r.category-1)
age.group <- age.group[,-5]
region <- model.matrix(~as.factor(replication$statename)-1)
region <- region[,-5]
sex <- replication$sex
model.r <- glm(replication$gb~age.group+sex+region+topsnp.r,family = binomial)
result.r <- summary(model.r)
places=2
log.or.r <- result.r$coefficients[11,1]
sd.r <- result.r$coefficients[11,2]
log.or.low.r <- -1.96*sd.r+log.or.r
log.or.high.r <- 1.96*sd.r+log.or.r
p.r <- result.r$coefficients[11,4]


or.r <- round(exp(log.or.r),places)
or.low.r <- round(exp(log.or.low.r),places)
or.high.r <- round(exp(log.or.high.r),places)

sd.trace <- function(OR,p){
  logor <- log(OR)
  sd <- logor/qnorm(1-0.5*p)
  return(sd)
}

or.d <- 1.51
p.d <- 3.81*1e-09
sd.d <- sd.trace(or.d,p.d)
log.or.d <- log(or.d)
log.or.low.d <- -1.96*sd.d+log.or.d
log.or.high.d <- 1.96*sd.d+log.or.d
or.low.d <- round(exp(log.or.low.d),places)
or.high.d <- round(exp(log.or.high.d),places)



beta.dis <- log.or.d
beta.r <- log.or.r
W.dis <- solve(sd.d^2)
W.r <- solve(sd.r^2)
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



final.table <- data.frame(snp.or=c(or.d,or.r,or.meta),
                          snp.or.low = c(or.low.d,or.low.r,or.low.meta),
                          snp.or.high = c(or.high.d,or.high.r,or.high.meta),
                          snp.p = c(p.d,p.r,p.meta)
                        
)
rownames(final.table) <- c("discovery","replicates","meta")

final.table.1 <- final.table


#############for snp rs17209837

topsnp.r <- replication[,3]
topsnp.r.new <- rep(0,length(topsnp.r))
topsnp.r.new[topsnp.r=="CT"] <- 1
topsnp.r.new[topsnp.r=="TT"] <- 2
topsnp.r <- topsnp.r.new


age.r <- replication[,9]
age.r.category <- as.factor(cut(age.r,c(17,29,39,49,60,100),labels = F))
age.group <- model.matrix(~age.r.category-1)
age.group <- age.group[,-5]
region <- model.matrix(~as.factor(replication$statename)-1)
region <- region[,-5]
sex <- replication$sex
model.r <- glm(replication$gb~age.group+sex+region+topsnp.r,family = binomial)
result.r <- summary(model.r)
places=2
log.or.r <- result.r$coefficients[11,1]
sd.r <- result.r$coefficients[11,2]
log.or.low.r <- -1.96*sd.r+log.or.r
log.or.high.r <- 1.96*sd.r+log.or.r
p.r <- result.r$coefficients[11,4]


or.r <- round(exp(log.or.r),places)
or.low.r <- round(exp(log.or.low.r),places)
or.high.r <- round(exp(log.or.high.r),places)

sd.trace <- function(OR,p){
  logor <- log(OR)
  sd <- logor/qnorm(1-0.5*p)
  return(sd)
}

or.d <- 1.66
p.d <- 1.99*1e-08
sd.d <- sd.trace(or.d,p.d)
log.or.d <- log(or.d)
log.or.low.d <- -1.96*sd.d+log.or.d
log.or.high.d <- 1.96*sd.d+log.or.d
or.low.d <- round(exp(log.or.low.d),places)
or.high.d <- round(exp(log.or.high.d),places)



beta.dis <- log.or.d
beta.r <- log.or.r
W.dis <- solve(sd.d^2)
W.r <- solve(sd.r^2)
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



final.table <- data.frame(snp.or=c(or.d,or.r,or.meta),
                          snp.or.low = c(or.low.d,or.low.r,or.low.meta),
                          snp.or.high = c(or.high.d,or.high.r,or.high.meta),
                          snp.p = c(p.d,p.r,p.meta)
                          
)
rownames(final.table) <- c("discovery","replicates","meta")

final.table.2 <- final.table





#############for snp rs4148808

topsnp.r <- replication[,4]
topsnp.r.new <- rep(0,length(topsnp.r))
topsnp.r.new[topsnp.r=="AG"] <- 1
topsnp.r.new[topsnp.r=="AA"] <- 2
topsnp.r <- topsnp.r.new


age.r <- replication[,9]
age.r.category <- as.factor(cut(age.r,c(17,29,39,49,60,100),labels = F))
age.group <- model.matrix(~age.r.category-1)
age.group <- age.group[,-5]
region <- model.matrix(~as.factor(replication$statename)-1)
region <- region[,-5]
sex <- replication$sex
model.r <- glm(replication$gb~age.group+sex+region+topsnp.r,family = binomial)
result.r <- summary(model.r)
places=2
log.or.r <- result.r$coefficients[11,1]
sd.r <- result.r$coefficients[11,2]
log.or.low.r <- -1.96*sd.r+log.or.r
log.or.high.r <- 1.96*sd.r+log.or.r
p.r <- result.r$coefficients[11,4]


or.r <- round(exp(log.or.r),places)
or.low.r <- round(exp(log.or.low.r),places)
or.high.r <- round(exp(log.or.high.r),places)

sd.trace <- function(OR,p){
  logor <- log(OR)
  sd <- logor/qnorm(1-0.5*p)
  return(sd)
}

or.d <- 1.71
p.d <- 2.37*1e-08
sd.d <- sd.trace(or.d,p.d)
log.or.d <- log(or.d)
log.or.low.d <- -1.96*sd.d+log.or.d
log.or.high.d <- 1.96*sd.d+log.or.d
or.low.d <- round(exp(log.or.low.d),places)
or.high.d <- round(exp(log.or.high.d),places)



beta.dis <- log.or.d
beta.r <- log.or.r
W.dis <- solve(sd.d^2)
W.r <- solve(sd.r^2)
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



final.table <- data.frame(snp.or=c(or.d,or.r,or.meta),
                          snp.or.low = c(or.low.d,or.low.r,or.low.meta),
                          snp.or.high = c(or.high.d,or.high.r,or.high.meta),
                          snp.p = c(p.d,p.r,p.meta)
                          
)
rownames(final.table) <- c("discovery","replicates","meta")

final.table.3 <- final.table

