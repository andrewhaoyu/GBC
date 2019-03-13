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
snp.list[n.snp+1] <- "PRS"
high[high>=8] <- 8
snp.list_f <- factor(snp.list,levels=rev(snp.list))
Type <- c(rep("ind",n.snp),"sum")
new.data <- data.frame(Type,snp.list_f,OR,high,low,stringsAsFactors = F)

library(ggplot2)
png("./result/gallbladder_cancer_result.png",height=20,width = 15,res=300,units="cm")
ggplot(new.data,aes(x=snp.list_f,y=OR,ymin=low,ymax=high,shape=Type,colour=Type))+
  geom_pointrange()+
  scale_colour_manual(values=c("#386cb0","#fdb462"))+
  #geom_line(yintercept=1,lty=2)+
  coord_flip()+
  theme_bw()+
  geom_hline(yintercept=1,size=1,lty=2)+
  theme(legend.position="none")+
  ylab("Odds ratio")+
  scale_y_continuous(breaks=c(0,1,2.5,5,8))+
  xlab("SNP")+
 # facet_grid(.~method)+
  ggtitle(paste0("Forest plot for odds ratio of gallbladder cancer")
  )+
  theme(plot.title = element_text(hjust=0.5,face="bold"),
        axis.text=element_text(face="bold"),
        axis.title.x = element_text(face="bold"),
        axis.title.y = element_text(face="bold"))
dev.off()


##stratified analysis
model.result <- read.csv("./result/Gallbladderstone_analysis_29SNPs.csv")[,-1]
OR_str <- as.character(model.result[,2])
n.snp <- nrow(snp.data)
low <- high <- OR <- rep(0,n.snp+1)
for(i in 1:n.snp){
  OR[i] <- as.numeric(strsplit(OR_str[i]," ") [[1]][1])
  temp <- strsplit(OR_str[i]," ")[[1]][2]
  low[i] <- as.numeric(gsub("\\[","",strsplit(temp,",")[[1]][1]))
  high[i] <- as.numeric(gsub("\\]","",strsplit(temp,",")[[1]][2]))
}
OR[n.snp+1] <- exp(PRS.result[2,1])
high[n.snp+1] <- exp(PRS.result[2,1]+1.96*PRS.result[2,2])
low[n.snp+1] <- exp(PRS.result[2,1]-1.96*PRS.result[2,2])
OR[c(6,10,29)] <- NA
high[c(6,10,29)] <- NA
low[c(6,10,29)] <- NA
#snp.list[n.snp+1] <- "PRS"
# high[high>=8] <- 8
# snp.list_f <- factor(snp.list,levels=rev(snp.list))
new.data1 <- data.frame(Type,snp.list_f,OR,high,low,stringsAsFactors = F)

model.result <- read.csv("./result/Gallbladdernostone_analysis_29SNPs.csv")[,-1]
OR_str <- as.character(model.result[,2])
n.snp <- nrow(snp.data)
low <- high <- OR <- rep(0,n.snp+1)
for(i in 1:n.snp){
  OR[i] <- as.numeric(strsplit(OR_str[i]," ") [[1]][1])
  temp <- strsplit(OR_str[i]," ")[[1]][2]
  low[i] <- as.numeric(gsub("\\[","",strsplit(temp,",")[[1]][1]))
  high[i] <- as.numeric(gsub("\\]","",strsplit(temp,",")[[1]][2]))
}
OR[n.snp+1] <- exp(PRS.result[3,1])
high[n.snp+1] <- exp(PRS.result[3,1]+1.96*PRS.result[3,2])
low[n.snp+1] <- exp(PRS.result[3,1]-1.96*PRS.result[3,2])
#snp.list[n.snp+1] <- "PRS"
high[high>=8] <- 8
# snp.list_f <- factor(snp.list,levels=rev(snp.list))
new.data2 <- data.frame(Type,snp.list_f,OR,high,low,stringsAsFactors = F)
new.data.final <- rbind(new.data1,new.data2)
label <- c(rep("stone",n.snp+1),rep("no stone",n.snp+1))
new.data.final <- cbind(new.data.final,label)
png("./result/gallbladder_cancer_strat_result.png",height=20,width = 20,res=300,units="cm")
ggplot(new.data.final,aes(x=snp.list_f,y=OR,ymin=low,ymax=high,shape=Type,colour=Type))+
  geom_pointrange()+
  scale_colour_manual(values=c("#386cb0","#fdb462"))+
  #geom_line(yintercept=1,lty=2)+
  coord_flip()+
  theme_bw()+
  geom_hline(yintercept=1,size=1,lty=2)+
  theme(legend.position="none")+
  ylab("Odds ratio")+
  scale_y_continuous(breaks=c(0,1,2.5,5,8))+
  xlab("SNP")+
  facet_grid(.~label)+
  ggtitle(paste0("Forest plot for odds ratio of gallbladder cancer")
  )+
  theme(plot.title = element_text(hjust=0.5,face="bold"),
        axis.text=element_text(face="bold"),
        axis.title.x = element_text(face="bold"),
        axis.title.y = element_text(face="bold"))
dev.off()
