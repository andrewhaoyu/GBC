# filename = "/dcl01/chatterj/data/GB/NEW/KG/snpR"
# files = dir(filename,pattern="result_new")
# n = 1000
# result_all = data.frame(a=rep("c",n),snpid=rep("c",n),position=rep(0,n),ALLELE1=rep("c",n),ALLELE1=rep("c",n),
# 	matrix(0,n,8253),    stringsAsFactors=F)
# total = 0
# for(i in 1:length(files)){
# 	print(i)
# 	load(files[i])
# 	temp = nrow(result)
# 	result_all[(1+total):(temp+total),]= result
# 	result_all[(1+total):(temp+total),2]= as.character(result[,2])
# 	total = temp+total
# }
# result_all =result_all[1:total,]
# 



commanarg <- commandArgs(trailingOnly = T)
print(commanarg)
i1 <- as.numeric(commanarg[[1]])
i2 <- as.numeric(commanarg[[2]])
print(i1)

#snplist.gtex = c("rs10004195","rs368433")
#snp.data <- read.csv("/users/hzhang1/R/GBV/data/Snps_pos_chr_extract_041120.csv",stringsAsFactors=F)
# chr.list <-  as.numeric(gsub("chr","",snp.data$Chromosome))
# pos.list <- snp.data$Position
trait = c("bmi","whr")

snp.data <- read.csv(paste0("/users/hzhang1/R/GBV/data/",trait[i2],"_snps.csv"),stringsAsFactors=F)
chr.list <-  as.numeric(snp.data$chr)
pos.list <- snp.data$position

#snp.data <- read.csv("/users/hzhang1/R/GBV/data/extracted_SNPs_information_new.csv",stringsAsFactors = F)
#three SNPs (rs1260326,rs756082276,rs756935975 ) can't be found in the dataset
# snp.data <- snp.data[-c(28:30),]
# snp.list <- snp.data$SNP
# chr.list <-  snp.data$CHR
# pos.list <- snp.data$Position
#names <- paste(snplist.gtex,collapse="|")

# pos.list <- c(87105795,87060844)
# 
# 
#try2 <- which(data$position%in%pos.list)
# length(try2)
# data[try2,]



setwd("/dcl01/chatterj/data/GB/NEW/KG/snpR")
filesDir <- getwd()
files <- dir(filesDir,pattern="imputed",full.names=T)
result <- NULL


file <- files[i1]
#load("/dcl01/chatterj/data/GB/NEW/KG/mergesnp.Rdata")
load(file)
chr.temp <- gsub("/dcl01/chatterj/data/GB/NEW/KG/snpR/chr","",file)
chr.temp <- gsub(".imputed.txt.gz.rda","",chr.temp)
chr <- as.numeric(strsplit(chr.temp,"_")[[1]][1])

chr.vec <- rep(chr,nrow(data))
pos.vec <- data[,3]
chr.pos <- paste0(chr.vec,":",pos.vec)
chr.pos.target <- paste0(chr.list,":",pos.list)
# snpid <- data$V2
# snpid <- as.character(snpid)
# idx <- grep("rs",snpid)
# data <- data[idx,]
# snpid <- data$V2
# snpid <- as.character(snpid)
# snpid <- gsub(":(.*)","",snpid)
try <- chr.pos%in%chr.pos.target
if(sum(try)!=0){
	result <- data[try,]
  	save(result,file=paste0("/dcl01/chatterj/data/GB/NEW/KG/snpR/extract_snps/",trait[i2],"_snps_result_extract",i1,".rda"))
}
print(sum(try))
print("finished")



