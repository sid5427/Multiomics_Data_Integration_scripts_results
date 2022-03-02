#install.packages("doRNG")
#install.packages("biocLite")
#BiocManager::install("GENIE3")
library(dplyr)
library(tidyverse)
library(tibble)
library(stringr)
#load('SS_vs_WW_lab_FR697_combined.RData')


FPKM_flower <- read.csv('annotated_new_flower_results_fpkm_10_18_2021_withkegg.txt',header=TRUE,sep='\t')
TF_list<- read.csv('TF_list.txt',header=TRUE,sep='\t')

#merge TF ids with main matrix
joined_FPKM_flower_TF<-left_join(FPKM_flower,TF_list,by="gene_id")
#write.table(joined_FPKM_flower_TF, file = "joined_FPKM_flower_TF.txt",sep="\t")

colnames(joined_FPKM_flower_TF)

## EPIC!!!
genes_of_interest<-filter(joined_FPKM_flower_TF,str_detect(arabi.symbol, "CYP707A") | str_detect(arabi.symbol, "STOMAGEN"))
###

write.csv(genes_of_interest, file = "genes_of_interest.csv")

#Filter out significant genes for C_vs_D comparison
C_vs_D_Significant_FPKM<-filter(joined_FPKM_flower_TF,C_vs_D_Significant == "YES")
#Extract only FPKM matrix with gene ids and TF ids for significant genes for C_vs_D comparison
C_vs_D_Significant_FPKM_with_TF<-C_vs_D_Significant_FPKM %>%
  arrange(TF_ID)%>%
            select(gene_id,TF_ID,gene.names.from.assembly,arabi.symbol,C1,C2,C3,D1,D2,D3,H1,H2,H3,HD1,HD2,HD3)

#Filter out significant genes for C_vs_H comparison
C_vs_H_Significant_FPKM<-filter(joined_FPKM_flower_TF,C_vs_H_Significant == "YES")
#Extract only FPKM matrix with gene ids and TF ids for significant genes for C_vs_D comparison
C_vs_H_Significant_FPKM_with_TF<-C_vs_H_Significant_FPKM %>%
  arrange(TF_ID) %>%
  select(gene_id,TF_ID,gene.names.from.assembly,arabi.symbol,C1,C2,C3,D1,D2,D3,H1,H2,H3,HD1,HD2,HD3)

#Filter out significant genes for C_vs_HD comparison
C_vs_HD_Significant_FPKM<-filter(joined_FPKM_flower_TF,C_vs_HD_Significant == "YES")
#Extract only FPKM matrix with gene ids and TF ids for significant genes for C_vs_D comparison
C_vs_HD_Significant_FPKM_with_TF<-C_vs_HD_Significant_FPKM %>%
  arrange(TF_ID) %>%
  select(gene_id,TF_ID,gene.names.from.assembly,arabi.symbol,C1,C2,C3,D1,D2,D3,H1,H2,H3,HD1,HD2,HD3)


#Filter out All groups significant genes
ALL_Significant_FPKM<-joined_FPKM_flower_TF %>% filter(C_vs_D_Significant == "YES"|
                             C_vs_H_Significant == "YES"| 
                             C_vs_HD_Significant == "YES")

#Extract only FPKM matrix with gene ids and TF ids for significant genes for C_vs_D comparison
ALL_Significant_FPKM_with_TF<-ALL_Significant_FPKM %>%
  arrange(TF_ID) %>%
  select(gene_id,TF_ID,gene.names.from.assembly,arabi.symbol,C1,C2,C3,D1,D2,D3,H1,H2,H3,HD1,HD2,HD3)




library(GENIE3)
set.seed(123) # For reproducibility of results

#C_vs_D_Significant_FPKM_genie_mt<-as.matrix(C_vs_D_Significant_FPKM_with_TF %>% 
ALL_Significant_FPKM_genie_mt<-as.matrix(ALL_Significant_FPKM_with_TF %>% 
  remove_rownames %>% 
  column_to_rownames(var="gene_id") %>%
  #select(C1,C2,C3,HD1,HD2,HD3))
  select(C1,C2,C3,D1,D2,D3,H1,H2,H3,HD1,HD2,HD3))

# Genes that are used as candidate regulators
#regulators <- c(2, 4, 7)
##change to the number of TFs in the fpkm filtered matrix##
regulators <- c(1:1316)
# Or alternatively:
#regulators <- c("Gene2", "Gene4", "Gene7")
weightMat <- GENIE3(ALL_Significant_FPKM_genie_mt, regulators=regulators, nCores=4, verbose=TRUE)
dim(weightMat)
weightMat[1:5,1:5]

##Get all the regulatory links
linkList <- getLinkList(weightMat)
dim(linkList)
head(linkList)

##Get only the top-ranked links
#linkList_max5 <- getLinkList(weightMat, reportMax=5)
#linkList_max200 <- getLinkList(weightMat, reportMax=200)
linkList_max3k <- getLinkList(weightMat, reportMax=3000)

write.csv(linkList_max3k, file = "ALL_max3k.csv")
##weight threshold
linkList_threshold <- getLinkList(weightMat, threshold=0.01)
write.csv(linkList_threshold, file = "ALL_thres_0.01.csv")

## remember me?
colnames(genes_of_interest)

##yes I know - I forgot to save them in their own data.frames....
C_vs_D_max3k<-read.csv('C_vs_D_max3k.csv',header=TRUE,sep=',')
C_vs_H_max3k<-read.csv('C_vs_H_max3k.csv',header=TRUE,sep=',')
C_vs_HD_max3k<-read.csv('C_vs_HD_max3k.csv',header=TRUE,sep=',')

C_vs_D_thres_0.01<-read.csv('C_vs_D_thres_0.01.csv',header=TRUE,sep=',')
C_vs_H_thres_0.01<-read.csv('C_vs_H_thres_0.01.csv',header=TRUE,sep=',')
C_vs_HD_thres_0.01<-read.csv('C_vs_HD_thres_0.01.csv',header=TRUE,sep=',')


colnames(C_vs_D_max3k)

C_vs_D_goi_TF<-inner_join(genes_of_interest,C_vs_D_max3k,by = c("gene_id"="targetGene"))
C_vs_H_goi_TF<-inner_join(genes_of_interest,C_vs_H_max3k,by = c("gene_id"="targetGene"))
C_vs_HD_goi_TF<-inner_join(genes_of_interest,C_vs_HD_max3k,by = c("gene_id"="targetGene"))
ALL_goi_TF<-inner_join(genes_of_interest,linkList_max3k,by = c("gene_id"="targetGene"))
ALL_goi_TF<-inner_join(ALL_goi_TF,joined_FPKM_flower_TF,by = c("regulatoryGene"="gene_id"))
write.csv(ALL_goi_TF, file = "results/ALL_goi_TF_max3K.csv")

C_vs_D_goi_TF_thres_0.01<-inner_join(genes_of_interest,C_vs_D_thres_0.01,by = c("gene_id"="targetGene"))
C_vs_H_goi_TF_thres_0.01<-inner_join(genes_of_interest,C_vs_H_thres_0.01,by = c("gene_id"="targetGene"))
C_vs_HD_goi_TF_thres_0.01<-inner_join(genes_of_interest,C_vs_HD_thres_0.01,by = c("gene_id"="targetGene"))
ALL_goi_TF_thres_0.01<-inner_join(genes_of_interest,linkList_threshold,by = c("gene_id"="targetGene"))

write.csv(C_vs_D_goi_TF_thres_0.01, file = "results/C_vs_D_goi_TF_thres_0.01.csv")
write.csv(C_vs_H_goi_TF_thres_0.01, file = "results/C_vs_H_goi_TF_thres_0.01.csv")
write.csv(C_vs_HD_goi_TF_thres_0.01, file = "results/C_vs_HD_goi_TF_thres_0.01.csv")
colnames(ALL_goi_TF_thres_0.01)
colnames(joined_FPKM_flower_TF)

C_vs_D_goi_TF_thres_0.01<-inner_join(C_vs_D_goi_TF_thres_0.01,joined_FPKM_flower_TF,by = c("regulatoryGene"="gene_id"))
write.csv(C_vs_D_goi_TF_thres_0.01, file = "results/C_vs_D_goi_TF_thres_0.01.csv")

ALL_goi_TF_thres_0.01<-inner_join(ALL_goi_TF_thres_0.01,joined_FPKM_flower_TF,by = c("regulatoryGene"="gene_id"))
write.csv(ALL_goi_TF_thres_0.01, file = "results/ALL_goi_TF_thres_0.01.csv")


