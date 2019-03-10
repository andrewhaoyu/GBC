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
print(i1)

#snplist.gtex = c("rs10004195","rs368433")
snp.list <- c("rs2291428","rs2290846","rs1800961","rs601338","rs708686",
              "rs28929474","rs34851490","rs1169288","rs13280055","rs56398830","rs174567","rs11012737","rs2469991","rs1935","rs17240268","rs12004","rs55971546","rs11641445","rs17138478","rs2292553","rs12968116","rs11887534","rs212100","rs12633863","rs4148808","rs6471717","rs686030","rs1260326","rs2070959","rs45575636")
chr.list <-  c(10,4,20,19,19,14,19,12,8,13,11,10,8,10,15,22,13,16,17,2,18
             ,2,19,3,7,8,9,2,2,7)
pos.list <- c(45958856,
              151199080,
              43042364,
              49206674,
              5840619,
              94844947,
              46384554,
              121416650,
              11522353,
              103701690,
              61593005,
              21849769,
              120347476,
              64927823,
              90347814,
              38877461,
              103718308,
              11639484,
              36073320,
              219146803,
              55322502,
              44066247,
              48376995,
              149211512,
              87105795,
              59377357,
              15304782,
              27730940,
              234602191,
              87060844)
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
	save(result,file=paste0("/dcl01/chatterj/data/GB/NEW/KG/snpR/extract_snps/result_extract",i1,".rda"))
}
print(sum(try))
print("finished")
