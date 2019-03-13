setwd('/Users/zhangh24/GoogleDrive/project/Nilanjan/Gallballadar_Cancer/')
snp.data <- read.csv("data/extracted_SNPs_information_new.csv",stringsAsFactors = F)
#three SNPs can't be found in the dataset
snp.data <- snp.data[-c(28:30),]
snp.list <- snp.data$SNP
chr.list <-  snp.data$CHR
pos.list <- snp.data$Position
PRS.result <- read.csv("./result/PRS_result_GBC.csv")[,-1]
model.result <- read.csv("./result/Gallbladder_analysis_29SNPs.csv")[,-1]
OR_str <- as.character(model.result[,2])
n.snp <- nrow(snp.data)
low <- high <- OR <- rep(0,n.snp+1)
for(i in 1:n.snp){
  OR[i] <- as.numeric(strsplit(OR_str[i]," ") [[1]][1])
  temp <- strsplit(OR_str[i]," ")[[1]][2]
  low[i] <- as.numeric(gsub("\\[","",strsplit(temp,",")[[1]][1]))
  high[i] <- as.numeric(gsub("\\]","",strsplit(temp,",")[[1]][2]))
}
OR[n.snp+1] <- exp(PRS.result[1,1])
high[n.snp+1] <- exp(PRS.result[1,1]+1.96*PRS.result[1,2])
low[n.snp+1] <- exp(PRS.result[1,1]-1.96*PRS.result[1,2])
