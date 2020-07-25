library(tidyverse)
library(dplyr)
library(reshape2)
library(ggplot2)
library(data.table)
library(mltools)

### Clear all previous data. Load in filtered counts_RIP_par_genelevel data
data<-read_csv("all_counts_RIP_par_genelevel")

#Breaking up data into Naive and IRF3 
naive_data<-data %>% 
  dplyr::select(gene_short_name,naive_mean_counts, normalized_Flag_naive, naive_5UTR, naive_Intron, naive_Exon, naive_3UTR)

IRF3_data<-data %>% 
  dplyr::select(gene_short_name,IRF3_mean_counts, normalized_Flag_IRF3, IRF3_5UTR, IRF3_Intron, IRF3_Exon, IRF3_3UTR)

summary(IRF3_data)

#### Run one dataset at a time
IRF3_data_3UTR<-IRF3_data #%>% filter(naive_3UTR >=1 )
quantile(IRF3_data_3UTR$IRF3_3UTR)

IRF3_data_3UTR[,"quantiles"] <- bin_data(IRF3_data_3UTR$IRF3_3UTR, bins=c(0, 1, 2, 5, 10, 20, 56), binType = "explicit")
### Look at the breakdown of the data based on number of 3'utr binding sites
plot(IRF3_data_3UTR$quantiles)

HighBinders_3UTR_Base <-IRF3_data_3UTR %>%
  filter(quantiles == "[20, 56]" ) ; 

HighMediumBinders_3UTR_Base <-IRF3_data_3UTR %>%
  filter(quantiles == "[10, 20)"); 

MediumBinders_3UTR_Base <- IRF3_data_3UTR%>%
  filter(quantiles =="[5, 10)")

MediumLowBinders_3UTR_Base<-IRF3_data_3UTR%>%
  filter(quantiles == "[2, 5)")

LowBinders_3UTR_Base<-IRF3_data_3UTR%>%
  filter(quantiles == "[1, 2)")

No_3UTR<- IRF3_data_3UTR %>% 
  filter(quantiles == "[0, 1)")

###################################
HighBinders_3UTR_Express<-HighBinders_3UTR_Base$normalized_Flag_IRF3
HighMediumBinders_3UTR_Express<-HighMediumBinders_3UTR_Base$normalized_Flag_IRF3
MidBinders_3UTR_Express<-MediumBinders_3UTR_Base$normalized_Flag_IRF3
MidLowBinders_3UTR_Express<-MediumLowBinders_3UTR_Base$normalized_Flag_IRF3
LowBinders_3UTR_Express<-LowBinders_3UTR_Base$normalized_Flag_IRF3
NoBinders_3UTR_Express<-No_3UTR$normalized_Flag_IRF3
AllTranscripts<-IRF3_data_3UTR$normalized_Flag_IRF3

#KS test or Mann Whitney Wilcoxon test 
x<-wilcox.test(NoBinders_3UTR_Express,HighMediumBinders_3UTR_Express, alternative= "less")
p.adjust(x$p.value, method= "bonferroni", n=2)
#### TRYING with ggplot for better looking graph
library("reshape2")
library("plyr")
library("ggplot2")

df <- data.frame(x = c(AllTranscripts,NoBinders_3UTR_Express,LowBinders_3UTR_Express,
                       MidLowBinders_3UTR_Express,MidBinders_3UTR_Express,HighMediumBinders_3UTR_Express,
                       HighBinders_3UTR_Express),ggg = factor(rep(1:7, c(13582,9127,1472,1820,842,286, 35))))

#Boxplot
ggplot(df, aes(x=as.factor(ggg), y=x, fill=ggg)) +
  geom_boxplot(width=.5) +
  ylim(-4,4) +
  theme_minimal()

#CDF plot
ggplot(df, aes(x, colour = ggg)) +
  coord_cartesian(xlim = c(-4, 4))+
  stat_ecdf()+
  scale_colour_hue(name="my legend", labels=c('AllTranscripts', 'No Binders', 
                                              'Low (1-2)','Medium_Low (2-5)', 'Medium (5-10)', 
                                              'HighMedium (10-20 )',  'High > 20' ))

#################### Filtering data based on specific transcript feature 
IntronSpecific_IRF3<-filter(IRF3_data, IRF3_data$IRF3_3UTR == 0 & IRF3_data$IRF3_Intron !=0 &
                         IRF3_data$IRF3_5UTR == 0 & IRF3_data$IRF3_Exon== 0, 
                       IRF3_data$IRF3_mean_counts > 8)

ThreeUTRSpecific_IRF3<-filter(IRF3_data, IRF3_data$IRF3_Intron ==0 & IRF3_data$IRF3_3UTR !=0 
                         & IRF3_data$IRF3_5UTR==0 & IRF3_data$IRF3_Exon==0,  
                         IRF3_data$IRF3_mean_counts > 8); 
dim(ThreeUTRSpecific)

FiveUTRSpecific<- filter(IRF3_data, IRF3_data$IRF3_3UTR ==0 & IRF3_data$IRF3_Intron == 0 & 
                           IRF3_data$IRF3_Exon ==0 & IRF3_data$IRF3_5UTR !=0, 
                         IRF3_data$IRF3_mean_counts > 8); dim(FiveUTRSpecific)

CodingSpecific<- filter(IRF3_data, IRF3_data$IRF3_3UTR ==0 & IRF3_data$IRF3_Intron == 0 & 
                          IRF3_data$IRF3_5UTR ==0 & IRF3_data$IRF3_Exon !=0, 
                        IRF3_data$IRF3_mean_counts > 8); dim(CodingSpecific)

Intron_and_ThreeUTR_IRF3<- IRF3_data %>% 
  filter(IRF3_data$IRF3_Intron !=0 & IRF3_data$IRF3_3UTR !=0) #& IRF3_data$IRF3_5UTR==0 & 
           IRF3_data$IRF3_Exon==0 & IRF3_data$IRF3_mean_counts > 8 ); 

NonTargets<-filter(IRF3_data, IRF3_data$IRF3_3UTR ==0 & IRF3_data$IRF3_Intron ==0 &
                     IRF3_data$IRF3_5UTR==0 & IRF3_data$IRF3_Exon==0 | 
                     IRF3_data$IRF3_mean_counts < 8); 
dim(NonTargets)

Multi_Targets<- IRF3_data %>% filter(IRF3_3UTR > 0 & IRF3_Intron > 0, 
                                     IRF3_data$IRF3_mean_counts > 8) 
Multi_Targets <- Multi_Targets %>% filter( IRF3_5UTR > 0| IRF3_Exon > 0)

dim(Multi_Targets)

AllmRNAdata<-filter(IRF3_data, IRF3_data$IRF3_5UTR !=0 | IRF3_data$IRF3_Intron !=0 |IRF3_data$IRF3_Exon !=0 | IRF3_data$IRF3_3UTR !=0); dim(AllmRNAdata)

Intron_Only_Express<-(IntronSpecific_IRF3$normalized_Flag_IRF3)
ThreeUTR_Only_Express<-(ThreeUTRSpecific_IRF3$normalized_Flag_IRF3)
FiveUTR_Only_Express<-(FiveUTRSpecific$normalized_Flag_IRF3)
Coding_Only_Express<-(CodingSpecific$normalized_Flag_IRF3)
ThreeIntron_Express<-(Intron_and_ThreeUTR_IRF3$normalized_Flag_IRF3)
All_RNAtargets_Express_0H<-(IRF3_data$normalized_Flag_IRF3)
NonTargets_Express<-(NonTargets$normalized_Flag_IRF3)
Multi_Targets_Express<- Multi_Targets$normalized_Flag_IRF3
AllmRNAtargets_Express<-(AllmRNAdata$normalized_Flag_IRF3)

x<-wilcox.test(NonTargets_Express,Multi_Targets_Express)
p.adjust(x$p.value, method= "bonferroni", n=2)

df <- data.frame(x = c(All_RNAtargets_Express_0H, NonTargets_Express,FiveUTR_Only_Express,Coding_Only_Express, ThreeUTR_Only_Express,Intron_Only_Express, ThreeIntron_Express, Multi_Targets_Express) 
                 ,ggg = factor(rep(1:8, c(13582,8087,19,42,1791, 871,1992,527))))
ggplot(df, aes(x=as.factor(ggg), y=x, fill=ggg)) +
  geom_boxplot(width= .5) +
  ylim(-3,4) +
  theme_minimal()

ggplot(df, aes(x, colour = ggg)) +
  coord_cartesian(xlim = c(-3, 3))+
  stat_ecdf()+
  theme_minimal() +
  scale_colour_hue(name="my legend", labels=c('AllTranscripts', 'NonTargets','FiveUTR_Only', 
                                              'Coding_Only',  'ThreeUTR_Only', 'Intron_Only',
                                              'Three and Intron Only', 'Multi'))
ggplot(df, aes(x, colour = ggg)) +
  geom_density()+
  theme_minimal() +
  scale_colour_hue(name="my legend", labels=c('AllTranscripts', 'NonTargets','FiveUTR_Only', 
                                              'Coding_Only','ThreeUTR_Only', 'Intron_Only', 
                                              'Three and Intron Only', 'Multi'))

### Compare to Intron (collecting all that are intronic)
Introns<-IRF3_data %>% filter(IRF3_data$IRF3_Intron >= 0 & IRF3_data$IRF3_3UTR == 0 & IRF3_data$IRF3_5UTR ==0 & IRF3_Exon == 0)

Introns[, "Intron_quantiles"] <- bin_data(Introns$IRF3_Intron, bins=c( 0, 1, 5,10,20, 100), binType = "explicit")
quantile(Introns$IRF3_Intron)

######
### Look at the breakdown of the data based on number of 3'utr binding sites
plot(Introns$Intron_quantiles)

HighBinders_Intron_Base <-Introns %>%
  filter(Intron_quantiles == "[20, 100]"); 

MidBinders_Intron_Base<- Introns %>% 
  filter(Intron_quantiles == "[10, 20)")

LowMediumBinders_Intron_Base <- Introns%>%
  filter(Intron_quantiles =="[5, 10)")

LowerBinders_Intron_Base<- Introns %>%
  filter(Intron_quantiles == "[1, 5)")

No_Intron<- Introns %>% 
  filter(Intron_quantiles == "[0, 1)")

HighBinders_Intron_Express<-(HighBinders_Intron_Base$normalized_Flag_IRF3)
MidBinders_Intron_Express<- MidBinders_Intron_Base$normalized_Flag_IRF3
LowMediumBinders_Intron_Express<-(LowMediumBinders_Intron_Base$normalized_Flag_IRF3)
LowerBinders_Intron_Express<-(LowerBinders_Intron_Base$normalized_Flag_IRF3)
NoBinders_Intron_Express<-No_Intron$normalized_Flag_IRF3
All_RNAs_Express<-Introns$normalized_Flag_IRF3

x<-wilcox.test(NoBinders_Intron_Express,HighBinders_Intron_Express, alternative= "less")
p.adjust(x$p.value, method= "bonferroni", n=2)
x

df <- data.frame(x = c(All_RNAs_Express,NoBinders_Intron_Express,LowerBinders_Intron_Express,
                       LowMediumBinders_Intron_Express, MidBinders_Intron_Express, HighBinders_Intron_Express) 
                 ,ggg = factor(rep(1:6, c(8958,8087,704,108,44, 15))))

ggplot(df, aes(x, colour = ggg)) +
  coord_cartesian(xlim = c(-4, 4))+
  stat_ecdf()+
  theme_minimal() +
  scale_colour_hue(name="my legend", labels=c('AllTranscripts', 'No_binders', 'Low', 
                                              'LowMedium', 'Mid', 'MidHigh', 'High'))

ggplot(df, aes(x=as.factor(ggg), y=x, fill=ggg)) +
  geom_boxplot(width=.5) +
  ylim(-4,4) +
  theme_minimal()

#Boxplot of number of binding sites vs. enrichment
ggplot(data, aes(x = as.factor(naive_3UTR), y = normalized_Flag_naive)) +
  #geom_violin(scale="area", fill="steelblue", alpha=.5) +
  geom_boxplot(width=.5, fill="gold")
  theme_minimal() +
  labs(x= "Number of 3' UTR binding sites") +
  labs(y= "RIP-seq enrichment")

ggplot(Introns, aes(x= as.factor(Intron_quantiles), y = normalized_Flag_IRF3)) +
  geom_violin(scale="area", fill="orange", alpha=.5) +
  geom_boxplot(alpha = 0.5, fill = "red", width=.1) +
  theme_minimal() +
  labs(x= "Number of binding sites") +
  labs(y= "RIP-seq enrichment")

alldata2 <- alldata2 %>% filter(IRF3_3UTR == 0)
lines(density(alldata2$FC.halflife), col="red")


### Now making CDF plots comparing exclusively bound 3'UTR or intronic transcripts to intron & 3'UTR.Comparing them based off number of 3'UTR sites 

ThreeUTRSpecific[,"quantiles"] <- bin_data(ThreeUTRSpecific$IRF3_3UTR, bins=c(0, 1, 2, 5, 10, 20, 56), binType = "explicit")
ThreeUTRSpecific <- ThreeUTRSpecific %>% mutate(source = "ThreeSpecific")

Intron_and_ThreeUTR[,"quantiles"] <- bin_data(Intron_and_ThreeUTR$IRF3_Intron, bins=c(0, 1, 5,10,20, 121), binType = "explicit")
Intron_and_ThreeUTR <- Intron_and_ThreeUTR %>% mutate(source = "Three and Intron")
All_3UTRs <- rbind(ThreeUTRSpecific, Intron_and_ThreeUTR)

glimpse(All_3UTRS_df)

ggplot(All_3UTRs, aes(x= normalized_Flag_IRF3, color = as.factor(source)))+
  stat_ecdf() +
  theme_minimal() + 
  facet_wrap(~quantiles)


Intron_and_ThreeUTR[,"quantiles"] <- bin_data(Intron_and_ThreeUTR$IRF3_Intron, bins=c(0, 1, 2, 5, 10, 20, 400), binType = "explicit")

I_Intron_and_ThreeUTR <- Intron_and_ThreeUTR %>% mutate(source = "I_Intron_and_ThreeUTR")

Introns<-IRF3_data %>% filter(IRF3_data$IRF3_Intron > 0 & IRF3_data$IRF3_3UTR == 0)
Introns[, "quantiles"] <- bin_data(Introns$IRF3_Intron, bins=c( 0, 1, 2, 5, 10, 20, 400)), binType = "explicit")
Introns <- Introns %>% mutate(source = "Introns")
All_Intron <- rbind(Intron_and_ThreeUTR, Introns)

ggplot(All_Intron, aes(x= normalized_Flag_IRF3, color = as.factor(source), na.rm = TRUE))+
  stat_ecdf() +
  theme_minimal() + 
  facet_wrap(~quantiles)


all_bound %>% 
  gather(RIPseq, source, c(normalized_Flag_naive,  normalized_Flag_IRF3)) %>%
                       ggplot(aes(x=source, color=RIPseq)) + geom_density() +
  theme_minimal() + xlim(-1,5)

all_data %>% mutate()


xx %>% group_by()
