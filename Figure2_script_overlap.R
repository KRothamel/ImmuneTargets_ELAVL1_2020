### We are testing push/pull 
library(tidyverse)
library(dplyr)
library(reshape2)
library(ggplot2)
library(data.table)
library(mltools)

### Clear all previous data. Load in filtered counts_RIP_par_genelevel data
data<-read_csv("all_counts_RIP_par_genelevel1")

#Breaking up data into Naive and IRF3 
naive_data<-data %>% 
  dplyr::select(gene_short_name,naive_mean_counts, normalized_Flag_naive, naive_5UTR, naive_Intron, naive_Exon, naive_3UTR)

IRF3_data<-data %>% 
  dplyr::select(gene_short_name,IRF3_mean_counts, normalized_Flag_IRF3, IRF3_5UTR, IRF3_Intron, IRF3_Exon, IRF3_3UTR)

IntronSpecific_naive<-naive_data %>% 
  filter(naive_3UTR == 0 & naive_Intron > 0 & naive_5UTR == 0 & naive_Exon== 0 & naive_mean_counts > 8); 

ThreeUTRSpecific_naive<-naive_data %>% 
  filter(naive_Intron == 0 & naive_3UTR > 0 & naive_5UTR==0 & naive_Exon==0, naive_mean_counts> 8); 

Intron_and_ThreeUTR_naive<- naive_data %>% 
  filter(naive_Intron !=0 & naive_3UTR !=0, naive_mean_counts > 8) 

IntronSpecific_IRF3<-filter(IRF3_data, IRF3_3UTR == 0 & IRF3_Intron !=0 & IRF3_5UTR == 0 & IRF3_Exon== 0, 
                            IRF3_mean_counts > 8)

ThreeUTRSpecific_IRF3<-filter(IRF3_data, IRF3_Intron ==0 & IRF3_3UTR !=0 & IRF3_5UTR==0 & IRF3_Exon==0,  
                              IRF3_mean_counts > 8); 

Intron_and_ThreeUTR_IRF3<- IRF3_data %>% 
  filter(IRF3_Intron !=0 & IRF3_3UTR !=0, IRF3_mean_counts > 8 )

# Overlap and find the intersecting and specific for each transcript group
`%notin%` = Negate(`%in%`)
x<- IntronSpecific_naive$gene_short_name %in% IntronSpecific_IRF3$gene_short_name
intron_naive_spec <-IntronSpecific_naive[x,]


x<- ThreeUTRSpecific_IRF3$gene_short_name %notin% ThreeUTRSpecific_naive$gene_short_name
intron_IRF3_spec <-ThreeUTRSpecific_IRF3[x,]
summary(intron_IRF3_spec)
write_csv(IRF3_intron_shared, "IRF3_intron_shared")

