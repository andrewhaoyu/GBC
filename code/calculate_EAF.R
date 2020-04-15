library(dplyr)
library(tidyr)
data <- read.csv("/Users/zhangh24/GoogleDrive/project/Nilanjan/Gallballadar_Cancer/result/extracted_snp_genotype_041120.csv",stringsAsFactors = F)
data = data %>% separate(Subject_id,c("famid","subid"),sep=":") %>% 
  mutate(subid=as.numeric(subid))
SNP <- colnames(data)[2:ncol(data)]
EAF <- colSums(data[,2:ncol(data)])/(2*nrow(data))

result <- data.frame(SNP,EAF,stringsAsFactors = F)
write.csv(result,"/Users/zhangh24/GoogleDrive/project/Nilanjan/Gallballadar_Cancer/result/extracted_snp_EAF_041320.csv",row.names = F)


post.menu <- read.csv("/Users/zhangh24/GoogleDrive/project/Nilanjan/Gallballadar_Cancer/data/post_menual.csv",stringsAsFactors = F)
library(dplyr)
colnames(post.menu)[1] <- "subid"

colnames(data)[2] <- "subid"

new.data <- left_join(post.menu,data,
                      by="subid")
#control only
idx <- which(new.data$status=="Control")

data.temp <- new.data[idx,]
SNP <- colnames(data.temp)[2:ncol(data)]
EAF.control <- colSums(data[,2:ncol(data)])/(2*nrow(data))


