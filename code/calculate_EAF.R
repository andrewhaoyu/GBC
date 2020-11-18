library(dplyr)
library(tidyr)
data <- read.csv("/Users/zhangh24/GoogleDrive/project/Nilanjan/Gallballadar_Cancer/result/extracted_snp_genotype_041420.csv",stringsAsFactors = F)
data = data %>% 
  separate(subject.id,c("famid","subid"),sep=":") %>% 
  mutate(subid=as.numeric(subid))
SNP <- colnames(data)[3:ncol(data)]
EAF <- colSums(data[,3:ncol(data)])/(2*nrow(data))




post.menu <- read.csv("/Users/zhangh24/GoogleDrive/project/Nilanjan/Gallballadar_Cancer/data/post_menual.csv",stringsAsFactors = F)
library(dplyr)
colnames(post.menu)[1] <- "subid"

colnames(data)[2] <- "subid"

new.data <- inner_join(post.menu,data,
                      by="subid")
#control only
idx <- which(new.data$status=="Control")

data.temp <- new.data[idx,]
SNP <- colnames(data.temp)[5:ncol(data)]
EAF.control <- colSums(as.matrix(data.temp[,5:ncol(data.temp)]))/(2*nrow(data.temp))
#control premenonly
idx <- which(new.data$status=="Control"&
               new.data$Menopause.status=="Pre-menopause")
EAF.control.premeno <- colSums(as.matrix(data.temp[,5:ncol(data.temp)]))/(2*nrow(data.temp))

#control post
idx <- which(new.data$status=="Control"&
               new.data$Menopause.status=="Post-menopause")
EAF.control.postmeno <- colSums(as.matrix(data.temp[,5:ncol(data.temp)]))/(2*nrow(data.temp))


#case only
idx <- which(new.data$status=="Case")

EAF.case <- colSums(as.matrix(data.temp[,5:ncol(data.temp)]))/(2*nrow(data.temp))
#case premenonly
idx <- which(new.data$status=="Case"&
               new.data$Menopause.status=="Pre-menopause")
EAF.case.premeno <- colSums(as.matrix(data.temp[,5:ncol(data.temp)]))/(2*nrow(data.temp))

#control post
idx <- which(new.data$status=="Case"&
               new.data$Menopause.status=="Post-menopause")
EAF.case.postmeno <- colSums(as.matrix(data.temp[,5:ncol(data.temp)]))/(2*nrow(data.temp))




result <- data.frame(SNP,EAF,
                     EAF.control,
                     EAF.control.premeno,
                     EAF.control.postmeno,
                     EAF.case,
                     EAF.case.premeno,
                     EAF.case.postmeno,
                     stringsAsFactors = F)
colnames(result) <- c("SNP",
                      "EAF_overall",
                      "EAF_control",
                      "EAF_control_premenopause",
                      "EAF_control_postmenopause",
                      "EAF_cases",
                      "EAF_cases_premenopause",
                      "EAF_cases_postmenopause")
 write.csv(result,"/Users/zhangh24/GoogleDrive/project/Nilanjan/Gallballadar_Cancer/result/extracted_snp_EAF_041420.csv",row.names = F)
