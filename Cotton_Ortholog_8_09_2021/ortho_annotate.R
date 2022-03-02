

library(dplyr)
library(tidyverse)

mapping_table = read.table("mapping.txt",header = TRUE,sep = "\t")

DAT_table = read.table("DATs_in_Ortholog_Group_Tables.txt",header = TRUE,sep = "\t")

##DAT_table$Gh_IDs,DAT_table$ZM_ids

##mapping_table$NBI_Gossypium_hirsutum_v1.1, mapping_table$Zea_mays.B73_RefGen_v4_orthofinder

##mapping_table$NBI_Gossypium_hirsutum_v1.1, mapping_table$Zea_mays.B73_RefGen_v4_OMA

#A$C %in% B$C

#if(DAT_table$Gh_IDs %in% mapping_table$NBI_Gossypium_hirsutum_v1.1)

mylist <- vector()
mylist_2 <- vector()

for (row_ids in 1:nrow(DAT_table)){
  #print(DAT_table[row_ids,3])
  #print(DAT_table[row_ids,6])
  if(DAT_table[row_ids,9] %in% mapping_table$GH_Zm_Orthofinder){
    mylist<-append(mylist,"Orthofinder")
  } else {mylist<-append(mylist,"NA")}
  
  if(DAT_table[row_ids,9] %in% mapping_table$GH_Zm_OMA){
    mylist_2<-append(mylist_2,"OMA")
  } else {mylist_2<-append(mylist_2,"NA")}
}

DAT_table$orthfinder_match <- mylist
DAT_table$OMA_match <- mylist_2

write.table(DAT_table,"DATs_in_Ortholog_Group_Tables_WITHGROUPMATCH.txt",sep = "\t")

#   if(DAT_table[row_ids,9] %in% mapping_table$GH_Zm_OMA){
#     print ("Orthofinder,OMA")
#   }
#   # if((DAT_table[row_ids,9] %in% mapping_table$GH_Zm_Orthofinder) && 
#   #    DAT_table[row_ids,9] %in% mapping_table$GH_Zm_OMA){
#   #     print(row_ids)
#   #   }
#   # if((DAT_table[row_ids,3] %in% mapping_table$NBI_Gossypium_hirsutum_v1.1) && 
#   #    (DAT_table[row_ids,6] %in% mapping_table$Zea_mays.B73_RefGen_v4_orthofinder) &&
#   #    (DAT_table[row_ids,6] %in% mapping_table$Zea_mays.B73_RefGen_v4_OMA)){
#   #   print(DAT_table[row_ids,3])
#   # }
# }









newflower_C_vs_H = read.table("leaf/C_vs_H_gene_exp.diff",header = TRUE,sep = "\t")
newflower_C_vs_HD = read.table("leaf/C_vs_HD_gene_exp.diff",header = TRUE,sep = "\t")

combined_table <- full_join(newflower_fpkm_table,newflower_C_vs_D,by = "gene_id")
combined_table <- full_join(combined_table,newflower_C_vs_H,by = "gene_id")
combined_table <- full_join(combined_table,newflower_C_vs_HD,by = "gene_id")


#write.table(annotated_counts,"combined_new_leaf_results.txt",sep = "\t",row.names = TRUE,col.names = TRUE)
#comb_table = read.csv("combined_new_leaf_results.txt",header = TRUE,sep = "\t")

anno_table = read.csv("Gmax_275_Wm82.a2.v1.annotation_info.txt",header = TRUE,sep = "\t")

#colnames(anno_table)
#reduced_anno_table<-select(anno_table, c(Gene_ID,Best.hit.arabi.name,arabi.symbol,arabi.defline))
#colnames(reduced_anno_table)

#distinct_anno_table<-distinct(reduced_anno_table)
distinct_anno_table_geneID<-distinct(anno_table,gene_id, .keep_all= TRUE)

length(unique(anno_table$gene_id))

annotated_counts_exp <- left_join(combined_table,distinct_anno_table_geneID,by = "gene_id")

write.table(annotated_counts_exp,"annotated_leaf_results_fpkm.txt",sep = "\t",row.names = TRUE,col.names = TRUE)

TF_leaf_HD_uniq = read.table("TF_leaf_HD_uniq.txt",header = TRUE,sep = "\t")
TF_annotated_exp_leaf <- left_join(TF_leaf_HD_uniq,annotated_counts_exp,by = "gene_id")
write.table(TF_annotated_exp_leaf,"TF_annotated_leaf.txt",sep = "\t",row.names = TRUE,col.names = TRUE)
