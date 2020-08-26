
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

#### Run one dataset at a time
naive_data_3UTR<-naive_data %>% filter(naive_Intron == 0 &  naive_5UTR == 0 & naive_Exon == 0)

naive_data_3UTR[,"quantiles"] <- bin_data(naive_data_3UTR$naive_3UTR, bins=c(0, 1, 2, 5, 10, 20, 56), binType = "explicit")

#plot(naive_data_3UTR$quantiles)

IRF3_data[,"quantiles"] <- bin_data(IRF3_data$IRF3_3UTR, bins=c(0, 1, 2, 5, 10, 20, 56), binType = "explicit")

### Look at the breakdown of the data based on number of 3'utr binding sites
HighBinders_3UTR_Base <-naive_data_3UTR %>%
  filter(quantiles == "[20, 56]" ) ; 

HighMediumBinders_3UTR_Base <-naive_data_3UTR %>%
  filter(quantiles == "[10, 20)"); 

MediumBinders_3UTR_Base <- naive_data_3UTR%>%
  filter(quantiles =="[5, 10)")

MediumLowBinders_3UTR_Base<-naive_data_3UTR%>%
  filter(quantiles == "[2, 5)")

LowBinders_3UTR_Base<-naive_data_3UTR%>%
  filter(quantiles == "[1, 2)")

No_3UTR<- naive_data_3UTR %>% 
  filter(quantiles == "[0, 1)")

###################################
HighBinders_3UTR_Express<-HighBinders_3UTR_Base$normalized_Flag_naive
HighMediumBinders_3UTR_Express<-HighMediumBinders_3UTR_Base$normalized_Flag_naive
MidBinders_3UTR_Express<-MediumBinders_3UTR_Base$normalized_Flag_naive
MidLowBinders_3UTR_Express<-MediumLowBinders_3UTR_Base$normalized_Flag_naive
LowBinders_3UTR_Express<-LowBinders_3UTR_Base$normalized_Flag_naive
NoBinders_3UTR_Express_naive<-No_3UTR$normalized_Flag_naive
AllTranscripts_naive<-naive_data_3UTR$normalized_Flag_naive

# Wilcox test (or KS test) significance between two groups 
x<-wilcox.test(NonTargets_Express,Multi_Targets_Express)
p.adjust(x$p.value, method= "bonferroni", n=2)
x$p.value

#### TRYING with ggplot for better looking graph
library("reshape2")
library("plyr")
library("ggplot2")

df <- data.frame(x = c(AllTranscripts_naive,NoBinders_3UTR_Express_naive, LowBinders_3UTR_Express,
                       MidLowBinders_3UTR_Express, MidBinders_3UTR_Express,
                       HighMediumBinders_3UTR_Express,
                       HighBinders_3UTR_Express ),
                 ggg = factor(rep(1:7, c(8441,7242,556,463,150,28, 2))))

# CDF plot 
ggplot(df, aes(x, colour = ggg)) +
  coord_cartesian(xlim = c(-2, 2))+
  stat_ecdf()+
  theme_minimal() +
  scale_colour_hue(name="my legend", labels=c('AllTranscripts', 'No Binders', 
                                              'Low (1-2)','Medium_Low (2-5)', 'Medium (5-10)', 
                                              'HighMedium (10-20 )',  'High > 20'))
#Boxplot 
ggplot(df, aes(x=as.factor(ggg), y=x, fill=ggg)) +
  geom_boxplot(width= .5) +
  ylim(-5,5) +
  theme_minimal()

#historgram 
ggplot(naive_data_3UTR, aes(x= normalized_Flag_naive, fill= quantiles)) + 
  geom_histogram() + 
  theme_minimal()

#################### Filtering data based on specific transcript feature 
IntronSpecific_naive<-filter(naive_data, naive_data$naive_3UTR == 0 & naive_data$naive_Intron  > 0 & naive_data$naive_5UTR == 0 & naive_data$naive_Exon== 0 & naive_data$naive_mean_counts > 8); 

ThreeUTRSpecific_naive<-filter(naive_data, naive_data$naive_Intron == 0 & naive_data$naive_3UTR > 0 & naive_data$naive_5UTR==0 & naive_data$naive_Exon==0, naive_data$naive_mean_counts> 8); 

FiveUTRSpecific<- filter(naive_data, naive_data$naive_3UTR ==0 & naive_data$naive_Intron == 0 & naive_data$naive_Exon ==0 & naive_data$naive_5UTR !=0); dim(FiveUTRSpecific)

CodingSpecific<- filter(naive_data, naive_data$naive_3UTR ==0 & naive_data$naive_Intron == 0 & naive_data$naive_5UTR ==0 & naive_data$naive_Exon !=0, naive_data$naive_mean_counts > 8); dim(CodingSpecific)

Intron_and_ThreeUTR_naive<- naive_data %>% 
  filter(naive_data$naive_Intron !=0 & naive_data$naive_3UTR !=0 & naive_data$naive_5UTR==0 & naive_data$naive_Exon==0, naive_data$naive_mean_counts > 8 ); dim(Intron_and_ThreeUTR)

NonTargets<-filter(naive_data, naive_data$naive_3UTR ==0 & naive_data$naive_Intron ==0 & naive_data$naive_5UTR==0 & naive_data$naive_Exon==0 | naive_data$naive_mean_counts < 8); 

Multi_Targets<- naive_data %>% filter(naive_3UTR > 0 & naive_Intron > 0, naive_data$naive_mean_counts > 8 )
Multi_Targets<- Multi_Targets %>% filter(naive_5UTR > 0 | naive_Exon > 0)

AllmRNAdata<-filter(naive_data, naive_data$naive_5UTR !=0 | naive_data$naive_Intron !=0 |naive_data$naive_Exon !=0 | naive_data$naive_3UTR !=0); dim(AllmRNAdata)

Intron_Only_Express<-(IntronSpecific_naive$normalized_Flag_naive) 
ThreeUTR_Only_Express<-(ThreeUTRSpecific_naive$normalized_Flag_naive)
FiveUTR_Only_Express<-(FiveUTRSpecific$normalized_Flag_naive)
Coding_Only_Express<-(CodingSpecific$normalized_Flag_naive)
ThreeIntron_Express<-(Intron_and_ThreeUTR_naive$normalized_Flag_naive)
All_RNAtargets_Express_0H<-(naive_data$normalized_Flag_naive)
NonTargets_Express<-(NonTargets$normalized_Flag_naive)
Multi_Targets_Express<- Multi_Targets$normalized_Flag_naive
AllmRNAtargets_Express<-(AllmRNAdata$normalized_Flag_naive)

#Wilcox test of difference in enrichment between groups of transcripts 
x<-wilcox.test(NonTargets_Express,ThreeIntron_Express)
p.adjust(x$p.value, method= "bonferroni", n=2)

df <- data.frame(x = c(All_RNAtargets_Express_0H,NonTargets_Express,FiveUTR_Only_Express,
                       Coding_Only_Express,ThreeUTR_Only_Express, Intron_Only_Express,
                       ThreeIntron_Express,Multi_Targets_Express ) 
                 ,ggg = factor(rep(1:8, c(13582,7230,15,27,1201,1236,2721,889))))

#boxplot
ggplot(df, aes(x=as.factor(ggg), y=x, fill=ggg)) +
  geom_boxplot(width= .5) +
  ylim(-3,4) +
  theme_minimal()

#CDF plot
ggplot(df, aes(x, color=ggg)) +
  stat_ecdf() +
  xlim(-3,4) +
  theme_minimal()
  
df_intron <- data.frame(x = c(All_RNAtargets_Express_0H,NonTargets_Express,
                        Intron_Only_Express) ,ggg = factor(rep(1:3, c(13582,7230,1236))))

ggplot(df_intron, aes(x, color=ggg)) +
  geom_density() +
  xlim(0,2)+
  theme_minimal()

x<-wilcox.test(NonTargets_Express, Intron_Only_Expres, alternative= "less")
p.adjust(x$p.value, method= "bonferroni", n=2)

### Compare to Intron 
naive_data_intron<-naive_data %>% filter(naive_Intron  >= 0 & naive_3UTR == 0 & naive_Exon == 0 & naive_5UTR == 0)
IRF3_data_intron<-IRF3_data %>% filter(IRF3_Intron  >= 0 & IRF3_3UTR == 0 & IRF3_Exon == 0 & IRF3_5UTR == 0)

naive_data_intron[, "Intron_quantiles"] <- bin_data(naive_data_intron$naive_Intron, bins=c( 0, 1, 5, 10, 20, 406), binType = "explicit")

IRF3_data_intron[, "Intron_quantiles"] <- bin_data(IRF3_data_intron$IRF3_Intron, bins=c( 0, 1, 5, 10, 50, 100, 406), binType = "explicit")

### Look at the breakdown of the data based on number of Intronic binding sites
plot(IRF3_data_intron$Intron_quantiles)
glimpse(IRF3_data_intron)

HighBinders_Intron_Base <-naive_data_intron %>%
  filter(Intron_quantiles == "[20, 406]" ) ; 

MediumBinder_Intron_Base<-naive_data_intron %>%
  filter(Intron_quantiles == "[10, 20)")

LowMediumBinders_Intron_Base <- naive_data_intron%>%
  filter(Intron_quantiles =="[5, 10)")

LowerBinders_Intron_Base<- naive_data_intron %>%
  filter(Intron_quantiles == "[1, 5)")

No_Intron<- naive_data_intron %>% 
  filter(Intron_quantiles == "[0, 1)")

HighBinders_Intron_Express<-(HighBinders_Intron_Base$normalized_Flag_naive )
MediumBinder_Intron_Express<-(MediumBinder_Intron_Base$normalized_Flag_naive)
LowMediumBinders_Intron_Express<-(LowMediumBinders_Intron_Base$normalized_Flag_naive)
LowerBinders_Intron_Express<-(LowerBinders_Intron_Base$normalized_Flag_naive)
NoBinders_Intron_Express<-(No_Intron$normalized_Flag_naive)
All_RNAs_Express<-naive_data_intron$normalized_Flag_naive

x<-wilcox.test(NoBinders_Intron_Express, HighBinders_Intron_Express, alternative= "less")
p.adjust(x$p.value, method= "bonferroni", n=2)
x

df <- data.frame(x = c(All_RNAs_Express, NoBinders_Intron_Express,LowerBinders_Intron_Express,
                       LowMediumBinders_Intron_Express,MediumBinder_Intron_Express,
                       HighBinders_Intron_Express) 
                 ,ggg = factor(rep(1:6, c(8466,7230, 811, 195, 117, 113))))

#CDF plots of number of introns and level of enrichment
ggplot(df, aes(x, colour = ggg)) +
  coord_cartesian(xlim = c(-4, 4))+
  stat_ecdf()+
  theme_minimal()+
  scale_colour_hue(name="my legend", labels=c('AllTranscripts','No_binders','Lower', 
                                              'Medium Low','Medium', 'High_Intron'))

#boxplot plots of number of introns and level of enrichment
ggplot(df, aes(x=as.factor(ggg), y=x, fill=ggg)) +
  geom_boxplot(width=.5) +
  ylim(-5,5) +
  theme_minimal()

#Boxplot of number of binding sites vs. enrichment
ggplot(naive_data_3UTR, aes(x = as.factor(quantiles), y = normalized_Flag_naive)) +
  geom_boxplot(width= .5,alpha = 0.7, fill = "steelblue") +
  theme_minimal() +
  labs(x= "Number of 3' UTR binding sites") +
  labs(y= "RIP-seq enrichment")

ggplot(naive_data, aes(x = as.factor(Intron_quantiles), y = normalized_Flag_naive)) +
  geom_boxplot(width= .5,alpha = 0.7, fill = "steelblue") +
  theme_minimal() +
  labs(x= "Number of 3' UTR binding sites") +
  labs(y= "RIP-seq enrichment")

ggplot(naive_data, aes(x= normalized_Flag_naive, color = as.factor(Intron_quantiles)))+
  stat_ecdf()

###comparison between naive and IRF3
ggplot(IRF3_data, aes(x = as.factor(IRF3_3UTR ), y = normalized_Flag_IRF3)) +
  geom_boxplot(alpha = 0.7, fill = "maroon") +
  theme_minimal() +
  labs(x= "Number of 3' UTR binding sites") +
  labs(y= "RIP-seq enrichment")

### plotted in bins
ggplot(naive_data_3UTR, aes(x = as.factor(quantiles), y = normalized_Flag_naive)) +
  geom_boxplot(alpha = 0.7, fill = "steelblue", width = 0.5) +
  theme_minimal() +
  labs(x= "Number of 3' UTR binding sites") +
  labs(y= "RIP-seq enrichment")

ggplot(IRF3_data, aes(x = normalized_Flag_IRF3, color = quantiles)) +
  stat_ecdf() +
  theme_minimal() +
  labs(x= "Number of 3' UTR binding sites") +
  labs(y= "RIP-seq enrichment")

#### compare to Introns
naive_data_intron<-naive_data %>% filter(naive_Intron  >= 0 & naive_3UTR == 0)
glimpse(naive_data_intron)

ggplot(naive_data_intron, aes(x = as.factor(Intron_quantiles), y = normalized_Flag_naive)) +
  geom_boxplot(alpha = 0.7, fill = "steelblue", width= .5) +
  ylim(-5,5) +
  theme_minimal() +
  labs(x= "Number of Intronic binding sites") +
  labs(y= "RIP-seq enrichment")
###comparison between naive and IRF3
ggplot(IRF3_data_intron, aes(x = as.factor(Intron_quantiles), y = normalized_Flag_IRF3)) +
  geom_boxplot(alpha = 0.7, fill = "maroon", width= .5) +
  theme_minimal() +
  labs(x= "Number of Intronic sites") +
  labs(y= "RIP-seq enrichment")

###################################
#plotting difference in intron/3'UTR fraction between Naive and IRF3
glimpse(data)

plot((data$IRF3_Intron/data$IRF3_3UTR), (data$naive_Intron/data$naive_3UTR), pch=20)
abline(a= 0, b =1)

ggplot(data, aes(x= (IRF3_Intron/IRF3_3UTR), y = (naive_Intron/naive_3UTR))) + 
  geom_point(color = "steelblue", alpha = .2) +
  theme_minimal() +
  geom_abline(slope = 1, intercept = 0)

glimpse(data %>% 
  mutate(naive_fraction = (naive_Intron/naive_3UTR), IRF3_fraction = (IRF3_Intron/IRF3_3UTR)) %>%
  gather(source,fraction, c(naive_fraction, IRF3_fraction))) %>%
  ggplot(aes(x= log2(fraction),  fill = source)) +
  geom_histogram() +
  theme_minimal() 

data %>% 
  gather(source, trans_feat, c(naive_3UTR, naive_Intron, IRF3_3UTR, IRF3_Intron)) %>%
  ggplot(aes(y = log2(trans_feat), x= source, fill= source)) +
  geom_boxplot(width= .5) +
  theme_minimal() 

shared <- data %>% filter(naive_Intron > 0 & naive_3UTR > 0 & IRF3_Intron > 0 & IRF3_3UTR > 0)

glimpse(shared %>% gather(source, fraction, c(naive_Intron, naive_3UTR, IRF3_Intron, IRF3_3UTR))) %>%
  ggplot(aes(y = log(fraction), x = as.factor(source))) +
  geom_boxplot()


#plotting enrichment based of of a cumulative number of "total" binding sites
naive_total<- all_counts_RIP_par %>% 
  mutate(total_number = rowSums(all_counts_RIP_par[7:10]))

naive_total[,"naive_total"] <- bin_data(naive_total$total_number, 
                                        bins=c(0, 1, 2, 5, 10,25, 415),
                                        binType = "explicit")
IRF3_total <- all_counts_RIP_par %>% 
  mutate(total_number = rowSums(all_counts_RIP_par[11:14]))

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
  stat_ecdf(aes(color= naive_total), color = "black") +
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

glimpse(alldata %>% mutate(category = as.numeric(IRF3_3UTR/sum(IRF3_5UTR:IRF3_3UTR)))) %>%
  ggplot(aes(x = halflife.x, y =category)) + geom_point() +  geom_smooth(method='lm', formula= y~x)

glimpse(IRF3_total)
zero <- naive_total %>% filter(naive_total == "[0, 1)")
test<- naive_total %>% filter(naive_total == "[25, 415]")

x<-wilcox.test(zero$normalized_Flag_naive, test$normalized_Flag_naive)
p.adjust(x$p.value, method = "bonferroni", n = 2) 


