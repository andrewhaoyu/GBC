filename = "/dcl01/chatterj/data/GB/NEW/KG/snpR"
files = dir(filename,pattern="result_new")
n = 1000
result_all = data.frame(a=rep("c",n),snpid=rep("c",n),position=rep(0,n),ALLELE1=rep("c",n),ALLELE1=rep("c",n),
	matrix(0,n,8253),    stringsAsFactors=F)
total = 0
for(i in 1:length(files)){
	print(i)
	load(files[i])
	temp = nrow(result)
	result_all[(1+total):(temp+total),]= result
	result_all[(1+total):(temp+total),2]= as.character(result[,2])
	total = temp+total
}
result_all =result_all[1:total,]




commanarg <- commandArgs(trailingOnly = T)
i1 <- as.numeric(commanarg[1])
print(i1)

snplist.gtex = c("rs10004195","rs368433")
#names <- paste(snplist.gtex,collapse="|")

setwd("/dcl01/chatterj/data/GB/NEW/KG/snpR")
filesDir <- getwd()
files <- dir(filesDir,pattern="imputed",full.names=T)
result <- NULL


file <- files[i1]
load(file)
snpid <- data$V2
snpid <- as.character(snpid)
idx <- grep("rs",snpid)
data <- data[idx,]
snpid <- data$V2
snpid <- as.character(snpid)
snpid <- gsub(":(.*)","",snpid)
try <- snpid%in%snplist.gtex
if(sum(try)!=0){
	result <- data[try,]
	save(result,file=paste0("/dcl01/chatterj/data/GB/NEW/KG/snpR/result_new",i1,".rda"))
}

