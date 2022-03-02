
library(dplyr)
library(tidyverse)

orthofinder_results = as.data.frame(read.table("G_hirs_vs_Zmays4_flatfile.txt",header=T))

OMA_results = as.data.frame(read.table("OMA_results_zm_gh.txt",header=T))

joined_datasets<-full_join(orthofinder_results,OMA_results,by = "NBI_Gossypium_hirsutum_v1.1")

write.table(joined_datasets, file = "orthofinder_with_OMA_results.txt", sep = "\t",
            col.names = TRUE,row.names = FALSE)
