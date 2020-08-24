### Figure 2 ELAVL1 cluster overlap 
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

#Breaking up into three categories mentioned in paper:intron exclusive, 3'UTR exclusive or both (naive)
IntronSpecific_naive<-naive_data %>% 
  filter(naive_3UTR == 0 & naive_Intron > 0 & naive_5UTR == 0 & naive_Exon== 0 & naive_mean_counts > 8); 

ThreeUTRSpecific_naive<-naive_data %>% 
  filter(naive_Intron == 0 & naive_3UTR > 0 & naive_5UTR==0 & naive_Exon==0, naive_mean_counts> 8); 

Intron_and_ThreeUTR_naive<- naive_data %>% 
  filter(naive_Intron !=0 & naive_3UTR !=0, naive_mean_counts > 8) 

#Breaking up into three categories mentioned in paper:intron exclusive, 3'UTR exclusive or both (stimulated)
IntronSpecific_IRF3<-filter(IRF3_data, IRF3_3UTR == 0 & IRF3_Intron !=0 & IRF3_5UTR == 0 & IRF3_Exon== 0, 
                            IRF3_mean_counts > 8)

ThreeUTRSpecific_IRF3<-filter(IRF3_data, IRF3_Intron ==0 & IRF3_3UTR !=0 & IRF3_5UTR==0 & IRF3_Exon==0,  
                              IRF3_mean_counts > 8); 

Intron_and_ThreeUTR_IRF3<- IRF3_data %>% 
  filter(IRF3_Intron !=0 & IRF3_3UTR !=0, IRF3_mean_counts > 8 )

# Overlap and find the intersecting and specific for each transcript group
`%notin%` = Negate(`%in%`)

x<- IntronSpecific_IRF3$gene_short_name %notin% IntronSpecific_naive$gene_short_name
intron_shared <-IntronSpecific_naive[x,]
intron_naive <-IntronSpecific_naive[x,]
intron_IRF3<- IntronSpecific_IRF3[x,]

y<- Intron_and_ThreeUTR_IRF3$gene_short_name %notin% Intron_and_ThreeUTR_naive$gene_short_name
Intron_and_ThreeUTR_shared <-Intron_and_ThreeUTR_naive[y,]
In_and_3_naive <-Intron_and_ThreeUTR_naive[y,]
In_and_3_IRF3<- Intron_and_ThreeUTR_IRF3[y,]

z<- ThreeUTRSpecific_IRF3$gene_short_name %notin% ThreeUTRSpecific_naive$gene_short_name
ThreeUTR_shared <-ThreeUTRSpecific_naive[z,]
Three_naive <-ThreeUTRSpecific_naive[z,]
Three_IRF3<- ThreeUTRSpecific_IRF3[z,]

