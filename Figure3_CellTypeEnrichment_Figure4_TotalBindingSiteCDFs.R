#Figure 3 
library(tidyverse)
library(dplyr)
library(reshape2)
library(ggplot2)
library(data.table)
library(mltools)

master_data<-read_csv("all_counts_RIP_par_genelevel1")

naive_mRNA_transcriptbound<-master_data %>%
  filter(naive_Exon > 0 | 
           naive_3UTR > 0 | 
           naive_Intron > 0 | 
           naive_5UTR > 0)

IRF3_mRNA_transcriptbound<-master_data %>% 
  filter(IRF3_Exon > 0 | 
           IRF3_3UTR > 0 | 
           IRF3_Intron > 0 | 
           IRF3_5UTR > 0 )

#write.csv(clusterdataIRF3, "clusterdataIRF3")
all_naive<- naive_mRNA_transcriptbound #%>% 
  #dplyr::select(gene_short_name, naive_mean_counts, normalized_Flag_naive, naive_5UTR, naive_Intron, naive_Exon, naive_3UTR)

all_IRF3<-IRF3_mRNA_transcriptbound #%>% 
  #dplyr::select(gene_short_name, IRF3_mean_counts, normalized_Flag_IRF3, IRF3_5UTR, IRF3_Intron, IRF3_Exon, IRF3_3UTR)

#Filtering by enrichment over IgG background
all_naive_enrich <- all_naive %>% filter(normalized_Flag_naive > 0)
all_IRF3_enrich <- all_IRF3 %>% filter(normalized_Flag_IRF3 > 0)

`%notin%` = Negate(`%in%`)
spec_naive_enrich <- all_naive_enrich$gene_short_name %notin% all_IRF3_enrich$gene_short_name
spec_naive_enrich <- all_naive_enrich[spec_naive_enrich,]

shared_enriched<- all_naive_enrich$gene_short_name %in% all_IRF3_enrich$gene_short_name
shared_enriched <- all_naive_enrich[shared_enriched,]

spec_IRF3_enrich <- all_IRF3_enrich$gene_short_name %notin% all_naive_enrich$gene_short_name
spec_IRF3_enrich <- all_IRF3_enrich[spec_IRF3_enrich,]

spec_IRF3_enrich_static<- spec_IRF3_enrich %>% filter(IRF3_mean_counts - naive_mean_counts < 1) %>% 
  filter(naive_mean_counts - IRF3_mean_counts < 1)

x<-spec_IRF3_enrich_static$gene_short_name %in% functional_targtets$gene
spec_IRF3_enrich_static_fun <- spec_IRF3_enrich_static[x,]


############################################################################
#plotting enrichment based of of a cumulative number of "total" binding sites 
#Figure 4 CDF plots on number of sites
master_data_cumulative<- master_data %>% 
  mutate(total_number_naive = rowSums(master_data[6:9])) %>%
  mutate(total_number_IRF3 = rowSums(master_data[10:13]))

master_data_cumulative[,"naive_total_quantiles"] <- bin_data(master_data_cumulative$total_number_naive, 
                                        bins=c(0, 1, 2, 5, 10,25, 415),
                                        binType = "explicit")

master_data_cumulative[,"IRF3_total_quantiles"] <- bin_data(master_data_cumulative$total_number_IRF3, 
                                      bins=c(0, 1, 2, 5, 10,25, 220), 
                                      binType = "explicit")


ggplot(naive_total, aes(color= naive_total, x = normalized_Flag_naive)) +
  stat_ecdf() +
  theme_minimal() +
  xlim(-2, 2)

ggplot(IRF3_total, aes(x = normalized_Flag_IRF3, color = IRF3_total)) +
  stat_ecdf() +
  xlim(-4,4)+
  theme_minimal()

ggplot(IRF3_total, aes(color= IRF3_total, x = normalized_Flag_IRF3)) +
  stat_ecdf() +
  theme_minimal() +
  xlim(-4, 4)

nontargets <- master_data_cumulative %>% filter(naive_total_quantiles == "[0, 1)")
greater_2 <- master_data_cumulative %>% filter(naive_total_quantiles != "[0, 1)" & naive_total_quantiles != "[1, 2)")

mean(greater_2$normalized_Flag_IRF3)/ mean(nontargets$normalized_Flag_IRF3)

2^mean(nontargets$normalized_Flag_naive)
2^mean(greater_2$normalized_Flag_naive)

all_IRF3 <- all_IRF3 %>% 
  mutate(regulation = (IRF3_mean_counts - naive_mean_counts > 1))

ggplot(all_IRF3, aes(x= normalized_Flag_IRF3, color = as.factor(regulation))) + stat_ecdf() #+ 
  stat_ecdf(aes( x= normalized_Flag_IRF3, color = as.factor(IRF3_3UTR )))

ggplot(all_IRF3, aes(x= naive_mean_counts, y = IRF3_mean_counts, color = normalized_Flag_IRF3))  + 
  geom_point() +
  theme_minimal()


#Read in Mukherjee's 2011 HEK293 PAR-CLIP data and compare to THP-1 data
HEK<- read_csv("molcel_3918_mmc3.csv")

HEK <- HEK %>% filter(Binding_Site_Bins != "No Sites")

#Filtering based of Mukherjee's enrichmend score "LOD"  
HEK_enrich<- HEK %>% filter(LOD_HuR > 0) 
xx<-HEK_enrich$Symbol %in% spec_IRF3_enrich$gene_short_name
xy<- HEK_enrich[xx,]

zz<-HEK$Symbol %in% spec_IRF3_enrich$gene_short_name
zy<- HEK[zz,]

write_csv(spec_IRF3_enrich_static, "xy")

plot(spec_IRF3_enrich$naive_mean_counts, spec_IRF3_enrich$IRF3_mean_counts, pch=20)
abline(a=0, b=1)
abline(a=1, b=1)
abline(a=-1, b=1)

###############################
# What portion out of all the binding sites are in the 3'UTR 
RIP_par_data<-read_csv("all_counts_RIP_par_genelevel1")

RIP_par_data %>% mutate(total_sites_naive = rowSums(RIP_par_data[6:9]), 
                        total_sites_IRF3= rowSums(RIP_par_data[10:13])) %>%
  ggplot(aes(y= IRF3_3UTR/total_sites_IRF3, x= naive_3UTR/total_sites_naive)) +
  geom_point() +
  theme_minimal() +
  geom_abline(x= 0, b= 1, color = "grey")

IRF3_total<-RIP_par_data %>% 
  mutate(total_number = rowSums(RIP_par_data[10:13])) 

summary(as.data.frame(naive_total$naive_total))

##### 
#plotting enrichment based of of a cumulative number of "total" binding sites 
naive_total<- RIP_par_data %>% 
  mutate(total_number = rowSums(RIP_par_data[6:9]))

naive_total[,"naive_total"] <- bin_data(naive_total$total_number, 
                                        bins=c(0, 1, 2, 5, 10,25, 415),
                                        binType = "explicit")

IRF3_total[,"IRF3_total"] <- bin_data(IRF3_total$total_number, 
                                      bins=c(0, 1, 2, 5, 10,25, 220), 
                                      binType = "explicit")

ggplot(naive_total, aes(x= naive_total, y = normalized_Flag_naive)) +
  geom_boxplot(width=.5) +
  theme_minimal()

ggplot(IRF3_total, aes(x= IRF3_total, y = normalized_Flag_IRF3)) +
  geom_boxplot(width=.5) +
  ylim(-4,4) +
  theme_minimal()

ggplot(naive_total, aes(color= naive_total, x = normalized_Flag_naive)) +
  stat_ecdf() +
  theme_minimal() +
  xlim(-2, 2)

ggplot(IRF3_total, aes(x = normalized_Flag_IRF3, color = IRF3_total)) +
  stat_ecdf() +
  xlim(-4,4)+
  theme_minimal()

ggplot(IRF3_total, aes(color= IRF3_total, x = normalized_Flag_IRF3)) +
  stat_ecdf() +
  theme_minimal() +
  xlim(-4, 4)

#Testing increasing number on introns vs. increasing number of 3'UTR in the transcripts that have both types of binding-sites
#############################################
ThreeUTRSpecific[,"quantiles"] <- bin_data(ThreeUTRSpecific$IRF3_3UTR, bins=c(0, 1, 2, 5, 10, 20, 56), binType = "explicit")
ThreeUTRSpecific <- ThreeUTRSpecific %>% mutate(source = "ThreeSpecific")

Intron_and_ThreeUTR[,"quantiles"] <- bin_data(Intron_and_ThreeUTR$IRF3_3UTR, bins=c(0, 1, 2, 5, 10, 20, 56), binType = "explicit")

Intron_and_ThreeUTR <- Intron_and_ThreeUTR %>% mutate(source = "Three and Intron")

All_3UTRs <- rbind(ThreeUTRSpecific, Intron_and_ThreeUTR)

ggplot(All_3UTRs, aes(x= halflife.x, color = as.factor(source)))+
  stat_ecdf() +
  theme_minimal() +  
  facet_wrap(~quantiles) 