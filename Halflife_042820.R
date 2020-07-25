BiocManager::install("GenomicFeatures")
library(RCAS)
library(ChIPseeker)
library(clusterProfiler)
library(tidyverse)
library(Biostrings)
library(biomaRt)
library(GenomicFeatures)
library(rtracklayer)

############# annotate bed file. 
gff <- import.gff("gencode.v19.chr_patch_hapl_scaff.annotation.gtf") #gtf from PAR-CLIP accessories

txdb <- GenomicFeatures::makeTxDbFromGRanges(gff)
txdbFeatures <- getTxdbFeaturesFromGRanges(gff)
transcriptCoords <- GenomicFeatures::transcripts(txdb)

#Importing the output file from alleyoop merge. File contains all UTR and T-to-C fractions per sample
data<-read_tsv("/Volumes/BackupPlus/alleyoop/merge_all", skip=3)
glimpse(data)
bed<-data %>% 
  dplyr::select(Chromosome, Start, End, Strand)
write_tsv(bed, "SLAM.bed")

# From RCAS, annotating the coordinates. Also where you get exon ID for ReactomePA GO visualizations later
library(org.Hs.eg.db)
mybed<-import("SLAM.bed")
peakAnno_SLAM <- annotatePeak(mybed, TxDb = txdb, annoDb = "org.Hs.eg.db")
peakAnno_SLAM_df <- as.data.frame(peakAnno_SLAM)
glimpse(peakAnno_SLAM_df)

my_anno<-peakAnno_SLAM_df %>% 
  dplyr::select(Chromosome= seqnames, Start= start, End= end, geneId, ENTREZID, SYMBOL,GENENAME)

glimpse(my_anno)

anno_SLAM<-left_join(data, my_anno, by="End" )
anno_SLAM<-na.omit(anno_SLAM)

full_Anno_SLAM<- anno_SLAM[!duplicated(anno_SLAM[c('SYMBOL')]),]
write_csv(full_Anno_SLAM, "full_Anno_SLAM")

##################################
library(tidyverse)
library(purrr)

#Reading in gene annotated SLAM-DUNK output
data<-read_csv("Full_Anno_SLAM")
glimpse(data)
#renaming columns 
data <- data %>% dplyr::rename(`0H.KO.1` = sample_0, `0H.TA.1`= sample_1, `0H.TA.2` = sample_2, `1H.KO.1` = sample_3 , 
                               `1H.TA.1`= sample_4, `3H.KO.1`=sample_5, `3H.KO.2`= sample_6, `3H.TA.1`= sample_7, 
                               `3H.TA.2`= sample_8, `6H.KO.1` = sample_9, `6H.KO.2`= sample_10, `6H.TA.1`= sample_11, 
                               `6H.TA.2`= sample_12, `8H.KO.1`= sample_13, `8H.TA.1`= sample_14, CTRL.1= sample_15, 
                               `1H.KO.2`= sample_16, CTRL.2= sample_17, `8H.TA.2`= sample_18)

tc.data<-data %>% 
  dplyr::select(SYMBOL, ENTREZID, geneId, "0H.KO.1","1H.KO.1","1H.KO.2", "3H.KO.1", 
                "3H.KO.2", "6H.KO.1", "6H.KO.2", "8H.KO.1","0H.TA.1","0H.TA.2","1H.TA.1",
                "3H.TA.1", "3H.TA.2","6H.TA.1", "6H.TA.2","8H.TA.1","8H.TA.2",
                CTRL.1, CTRL.2)

TA_raw_data <-tc.data %>%
  dplyr::select(SYMBOL, geneId, contains("TA")) %>% 
  mutate_each(funs(as.numeric), contains("TA"))

#Rearranging data and normalizing data to chase onset (0H)
TA_norm <- TA_raw_data %>%
  gather(key = Condition, value = T2C, -geneId, -`0H.TA.2`, -SYMBOL) %>%
  mutate(norm_T2C = T2C/`0H.TA.2`) %>%
  dplyr::select(-T2C) %>%
  spread(key = Condition, value = norm_T2C) %>%
  gather(key = Condition, value = norm_T2C, -geneId, -SYMBOL) %>%
  mutate(timepoints = as.numeric(str_sub(Condition, start = 1, end = 1))) %>%
  filter(!is.na(norm_T2C) & norm_T2C > 0)

KO_raw_data <-tc.data %>%
  dplyr::select(SYMBOL, geneId, contains("KO")) %>% 
  mutate_each(funs(as.numeric), contains("KO"))

KO_norm <- KO_raw_data %>%
  gather(key = Condition, value = T2C, -geneId, -`0H.KO.1`, -SYMBOL) %>%
  mutate(norm_T2C = T2C/`0H.KO.1`) %>%
  dplyr::select(-T2C) %>%
  spread(key = Condition, value = norm_T2C) %>%
  gather(key = Condition, value = norm_T2C, -geneId, -SYMBOL) %>%
  mutate(timepoints = as.numeric(str_sub(Condition, start = 1, end = 1))) %>%
  filter(!is.na(norm_T2C) & norm_T2C > 0)

####### Non-linear decay model equation from Herzog and tidied up by Sarah Arcos 
library(minpack.lm)

calculate_model2 <- function(gene_data){
  tryCatch({
    model = nlsLM(norm_T2C~Plat + (y0 - Plat) * exp(-k * (timepoints)),
                  data = gene_data,
                  start=list(y0= 1,
                             Plat = 0,
                             k= 0.5),
                  upper = c(1,0,Inf),
                  lower = c(1,0,0),
                  control = nls.lm.control(maxiter = 1000),
                  na.action = na.omit)
    return(broom::augment(model, gene_data) %>%
             mutate(halflife = log(2)/coef(model)["k"]))
  }, error = function(err){gene_data %>%
      mutate(.fitted = NA, .resid = NA, halflife = NA)})
}

result_TA <- TA_norm %>%
  group_by(geneId) %>%
  group_modify(~ calculate_model2(.)) %>%
  bind_rows()

result_TA_halflife <- result_TA %>%
  dplyr::select(SYMBOL, halflife) %>%
  unique()

result_KO <- KO_norm %>%
  group_by(geneId) %>%
  group_modify(~ calculate_model2(.)) %>%
  bind_rows()

result_KO_halflife <- result_KO %>%
  dplyr::select(SYMBOL, halflife) %>%
  unique()

# Herzog et al. 2017 does this in the SLAM-seq paper. Max half-life is 24 hours. 
result_TA_halflife$halflife[result_TA_halflife$halflife > 24 ] <- 24 
result_KO_halflife$halflife[result_KO_halflife$halflife > 24 ] <- 24 

result_TA_halflife<- na.omit(result_TA_halflife)
result_KO_halflife <- na.omit(result_KO_halflife)

#write.csv(result_TA_halflife, "result_TA_halflife")
#write.csv(result_KO_halflife, "result_KO_halflife")
############################

#######################
#Combining half-life data to PARCLIP, RIPseq and RNAseq data 
result_TA_halflife<- read_csv("result_TA_halflife")
result_KO_halflife<- read_csv("result_KO_halflife")

PARdata<-read_csv("all_counts_RIP_par_genelevel")
glimpse(PARdata)
colnames(PARdata)[2]<-"SYMBOL"

all_halflife<- full_join(result_TA_halflife, result_KO_halflife, by = "SYMBOL")
alldata<-full_join(PARdata, all_halflife, by = "SYMBOL")
alldata<-na.omit(alldata)
glimpse(alldata)

#Merging PAR-CLIP data and SLAM-seq RNA halflife data
IRF3_all_halflife<-full_join(PARdata, result_TA_halflife, by = "SYMBOL")
IRF3_all_halflife<-na.omit(IRF3_all_halflife)

#genes with a change in halflife by atleast 1.5 foldchange
functional_targets<- alldata %>% 
  filter(IRF3_3UTR > 0) %>% 
  filter(halflife.x/halflife.y > 1.5)


##### Looking at metadata of the functional targets. What is the expression and enrichment level of functional targets and how do they compare to non-functional targets
func<-func %>% filter(IRF3_mean_counts/naive_mean_counts > 0.6)
alldata_df<-cbind(alldata, func_df)
func_df<- alldata$SYMBOL %in% func$SYMBOL

alldata_df %>% filter(IRF3_3UTR > 0) %>% 
ggplot(aes(x= normalized_Flag_IRF3, y = IRF3_mean_counts/naive_mean_counts, color= as.factor(func_df))) + geom_point() + 
  theme_minimal() +
  scale_color_manual(values=c("#999999", "#E69F00", "#56B4E9")) +
  geom_hline(yintercept= 0) +
  geom_vline(xintercept= 0)

#write_tsv(test, "test")

# Boxplot of halflife from parental and KO (3 bins, all transcripts, targets and nontargets)
#halflife.x = parental, halflife.y = KO
alldata <- alldata %>% gather(source, halflife, c(halflife.x, halflife.y))

alldata %>%
    mutate(Category = case_when(
    IRF3_3UTR > 0  ~ " > 0 3'UTR sites", 
    IRF3_3UTR == 0 ~ " 0 3'UTR sites"
  )) %>%
    ggplot(aes(x = Category, y = halflife)) +
    geom_jitter(aes(x = "mean", y = halflife), alpha = .5 , width = .25, size = .5, color = "steelblue") +
    geom_boxplot(aes(x = "mean", y = halflife ), width = .5, size = .5, alpha= 0, color = "black") +
    geom_jitter(width = .25, alpha = .5, size = .5, color = "steelblue") +
    geom_boxplot(aes(x = Category, y = halflife ), width = .5, alpha = 0 , color= "black") +
    theme_minimal() 
  
# Simple boxplot and  violin (playing around with aesthitics)
alldata %>% filter (source == "halflife.y") %>%
    mutate(Category = case_when(
      IRF3_3UTR > 0  ~ " > 0 3'UTR sites", 
      IRF3_3UTR == 0 ~ " 0 3'UTR sites"
    )) %>%
    ggplot(aes(x = Category, y = halflife)) +
    geom_boxplot(aes(x = "mean", y = halflife ), width = .5, color = "black") +
    #geom_violin(aes(x = "mean", y = halflife ),  color = "black", alpha = 0, width = .4) +
    geom_boxplot(aes(x = Category, y = halflife ), width = .5, alpha = 0 , color= "black") +
    #geom_violin(aes(x = Category, y = halflife ), alpha = 0, color = "black", width = .4) +
    theme_minimal() +
  ylim(2,16)
  
#Testing the signiifcance Mann-Whitney wilcoxon of the differences in halflife between bins 
#Change halflife.x or y depending on parental or KO 
halflife_nosites_TA <- alldata %>% 
    filter(source == "halflife.y") %>% 
    filter(IRF3_3UTR == 0)

halflife_onesite_atleast_TA <- alldata %>% 
    filter(source == "halflife.y") %>% 
    filter(IRF3_3UTR > 0)
  
halflife_mean_TA <- alldata %>% 
    filter(source == "halflife.y")
  
#ks.test(halflife_nosites_TA$halflife, halflife_onesite_atleast_TA$halflife)

x<-wilcox.test(halflife_nosites_TA$halflife, halflife_onesite_atleast_TA$halflife,
               alternative = "less", conf.int = TRUE)

p.adjust(x$p.value, method = "bonferroni", n = 2) 

# CDF plot of 3 bins (all transcripts, non-targets and  3'UTR targets)
alldata %>% filter(source == "halflife.x") %>%
  mutate(Category = case_when(
    IRF3_3UTR > 0  ~ " > 0 3'UTR sites", 
    IRF3_3UTR == 0 ~ " 0 3'UTR sites"
  )) %>%
    ggplot(aes(x = halflife,  color = Category)) +
    stat_ecdf() +
    stat_ecdf(aes(x= halflife), color = "black") +
    theme_minimal() +
    xlim(2,16) +
    geom_line(aes(y = .5), linetype = "dotted")

#############################################
#Binning cimulative number of 3'UTR binding sites for CDF plots
library(mltools)
alldata[,"quantiles_Naive"] <- bin_data(alldata$naive_3UTR, bins=c(0, 1, 2, 5, 10, 56), binType = "explicit")

alldata[,"quantiles_IRF3"] <- bin_data(alldata$IRF3_3UTR, bins=c(0, 1, 2, 5, 10, 56), binType = "explicit")

### Look at the breakdown of the data based on number of 3'utr binding sites. Log2 FC of naive/KO
alldata %>% 
  spread(source, halflife) %>%
  ggplot(aes(x = log2(halflife.x / halflife.y), color = as.factor(quantiles_IRF3))) +
  stat_ecdf() +
  stat_ecdf(aes(x= log2(halflife.x / halflife.y)), color = "black") +
  xlim(-3, 3) +
  theme_minimal()

alldata %>%
  spread(source, halflife) %>%
  ggplot(aes(x= as.factor(quantiles_IRF3),  y= log2(halflife.x/halflife.y))) +
  geom_boxplot( width = .5) +
  ylim(-2, 2) +
  theme_minimal() 

# Look at KO data (halflife.y = ELAVL1 KO RNA halflives)
alldata %>% filter(source == "halflife.y") %>%
  ggplot(aes(x=as.factor(quantiles_IRF3), y = halflife))+
  geom_boxplot(width=.5) +
  ylim(4,16)+
  theme_minimal()

alldata %>% spread(source, halflife) %>%
  ggplot(aes(x= log2(halflife.x/halflife.y), color = as.factor(quantiles_IRF3)))+
  stat_ecdf() +
  xlim(-1.5,1.5) +
  theme_minimal()

alldata %>% spread(source, halflife) %>%
  ggplot(aes(y= log2(halflife.x/halflife.y), x= as.factor(quantiles_IRF3)))+
  geom_boxplot(width=.5) +
  ylim(-2, 2) +
  theme_minimal()

# splitting up transcripts into bins to run Mann-Whitney wilcoxon test 
q0 <- alldata %>% filter(quantiles_IRF3 == "[0, 1)" ) %>%
  spread(source, halflife)

q1 <- alldata %>% filter(quantiles_IRF3 == "[1, 2)" ) %>%
spread(source, halflife)

q2<- alldata %>% filter(quantiles_IRF3 == "[2, 5)" ) %>%
spread(source, halflife)

q3<- alldata %>% filter(quantiles_IRF3 == "[5, 10)") %>%
spread(source, halflife)

q4<- alldata %>% filter(quantiles_IRF3 == "[10, 56]") %>%
spread(source, halflife)

x<-wilcox.test((q0$halflife.x/q0$halflife.y), (q4$halflife.x/q4$halflife.y), paired = FALSE, alternatie= "less")
p.adjust(x$p.value, method = "bonferroni", n = 2) 

#############################################################
#Making CDFs
library(mltools)

threeutr_only<-alldata %>%filter(IRF3_3UTR > 0 & IRF3_Intron == 0)

alldata[,"quantiles_Naive"] <- bin_data(alldata$naive_3UTR, bins=c(0, 1, 2, 5, 10, 20, 56), binType = "explicit")

alldata[,"quantiles_IRF3"] <- bin_data(alldata$IRF3_3UTR, bins=c(0, 1, 2, 5, 10, 20, 56), binType = "explicit")

alldata <- alldata %>%
  spread(source, halflife)

HighBinders_3UTR_Base <-alldata %>%
  filter(quantiles_IRF3 == "[20, 56]" ) ; 

HighMediumBinders_3UTR_Base <-alldata %>%
  filter(quantiles_IRF3 == "[10, 20)"); 

MediumBinders_3UTR_Base <- alldata %>%
  filter(quantiles_IRF3 =="[5, 10)")

MediumLowBinders_3UTR_Base<-alldata %>%
  filter(quantiles_IRF3 == "[2, 5)")

LowBinders_3UTR_Base<-alldata %>%
  filter(quantiles_IRF3 == "[1, 2)")

No_3UTR<- alldata %>% 
  filter(quantiles_IRF3 == "[0, 1)")

###################################
HighBinders_3UTR_Express<-HighBinders_3UTR_Base$halflife.x
HighMediumBinders_3UTR_Express<-HighMediumBinders_3UTR_Base$halflife.x
MidBinders_3UTR_Express<-MediumBinders_3UTR_Base$halflife.x
MidLowBinders_3UTR_Express<-MediumLowBinders_3UTR_Base$halflife.x
LowBinders_3UTR_Express<-LowBinders_3UTR_Base$halflife.x
NoBinders_3UTR_Express<-No_3UTR$halflife.x
AllTranscripts<-alldata$halflife.x

#### TRYING with ggplot for better looking graph
library("reshape2")
library("plyr")
library("ggplot2")

df <- data.frame(x = c(AllTranscripts, HighBinders_3UTR_Express, 
                       HighMediumBinders_3UTR_Express,MidBinders_3UTR_Express, 
                       MidLowBinders_3UTR_Express, LowBinders_3UTR_Express, 
                       NoBinders_3UTR_Express),
                 ggg = factor(rep(1:7, c(4058,13,106,338,746, 652,2203))))


ggplot(df, aes(x, colour = ggg)) +
  stat_ecdf() +
  xlim(2, 18) +
  scale_colour_hue(name="my legend", labels=c('AllTranscripts', 'High > 20', 
                                              'HighMedium (10-20 )', 'Medium (5-10)',
                                              'Medium_Low (2-5)', 'Low (1-2)','No Binders'))

##############################################################
alldata %>% 
  gather(source, halflife, c(halflife.x, halflife.y)) %>%
  filter(source == "halflife.x") %>%
  ggplot(aes(x = as.factor(quantiles_IRF3), y = halflife)) +
  geom_boxplot(width = .5) +
  #geom_boxplot(aes(x = "mean", y = halflife), width = .5) +
  ylim(2,18) +
  theme_minimal()

ggplot(alldata, aes(x = halflife.x, color = quantiles_IRF3)) +
  stat_ecdf() +
  xlim(2,18) +
  theme_minimal()

alldata %>%  
gather(source, halflife, c(halflife.x, halflife.y)) %>% 
filter(IRF3_3UTR > 0) %>%
ggplot(aes(x= halflife)) +
  geom_density(aes(color = source), size = 1) +
  geom_density(aes(x= halflife), color = "black", size = 1) +
  theme_minimal() 

###############################################3
#Making fitted timecourse visualizations 
fitted_results<-full_join(result_TA, result_KO, by= c("SYMBOL", "timepoints"))
fitted_results<- na.omit(fitted_results)

fitted_results %>% 
  ggplot(aes(x= timepoints, y = .fitted.y, group= SYMBOL)) +
  geom_line(alpha = .05, color = "maroon") +
  theme_minimal()

fitted_results %>% filter(SYMBOL ==	"IFNB1") %>%
  gather(source, fitted, c(.fitted.x, .fitted.y)) %>%
  gather(resid_source, resid, c(.resid.x, .resid.y)) %>%
  ggplot(aes(x= timepoints, y = fitted, color = source)) +
  geom_line() +
  geom_point(aes(timepoints, y= fitted))+
  ylim(.2,1)+
  theme_minimal()

test_alldata<- alldata %>% 
  spread(source, halflife) %>%
  #spread(source, halflife) %>% 
  mutate(FC = halflife.x/halflife.y) %>%
  filter(IRF3_3UTR > 0) %>%
  filter(FC > 1.5)

non_functional<- alldata %>% 
  spread(source, halflife) %>%
  #spread(source, halflife) %>% 
  mutate(FC = halflife.x/halflife.y) %>%
  filter(IRF3_3UTR > 0) %>%
  filter(FC < 1.5)

#write.csv(non_functional,"non_functional")
test<-fitted_results %>% 
  mutate(FC_halflife= halflife.x/halflife.y)
  
ggplot(test, aes(x= timepoints, y= timepoint_mean.x)) +
  geom_line() +
  geom_line(aes(x= timepoints, y = timepoint_mean.y), color = "red") +
  theme_minimal() +
  ylim(0,1)

####################
#plotting CDF plot with stability and transcript feature
IntronSpecific<-filter(alldata, alldata$IRF3_3UTR == 0 & alldata$IRF3_Intron != 0 & 
                         alldata$IRF3_5UTR == 0 & alldata$IRF3_Exon == 0)

ThreeUTRSpecific<-filter(alldata, alldata$IRF3_Intron ==0 & alldata$IRF3_3UTR !=0 & alldata$IRF3_5UTR==0 & alldata$IRF3_Exon==0); 
dim(ThreeUTRSpecific)

FiveUTRSpecific<- filter(alldata, alldata$IRF3_3UTR ==0 & alldata$IRF3_Intron == 0 & alldata$IRF3_Exon ==0 & alldata$IRF3_5UTR !=0); 
dim(FiveUTRSpecific)

CodingSpecific<- filter(alldata, alldata$IRF3_3UTR ==0 & alldata$IRF3_Intron == 0 & alldata$IRF3_5UTR ==0 & alldata$IRF3_Exon !=0); 
dim(CodingSpecific)

Intron_and_ThreeUTR<- alldata %>% 
  filter(alldata$IRF3_Intron !=0 & alldata$IRF3_3UTR !=0 & alldata$IRF3_5UTR==0 & alldata$IRF3_Exon==0 ); dim(Intron_and_ThreeUTR)

NonTargets<-filter(alldata, alldata$IRF3_3UTR ==0 & alldata$IRF3_Intron ==0 & alldata$IRF3_5UTR==0 & alldata$IRF3_Exon==0); 
dim(NonTargets)

AllmRNAdata<-filter(alldata, alldata$IRF3_5UTR !=0 | alldata$IRF3_Intron !=0 |alldata$IRF3_Exon !=0 | alldata$IRF3_3UTR !=0); 
dim(AllmRNAdata)

Multi_Targets<-filter(alldata, alldata$IRF3_3UTR == !0 & alldata$IRF3_Intron == !0 & alldata$IRF3_5UTR== !0 | alldata$IRF3_Exon==!0); 
dim(Multi_Targets)

Intron_Only_Express<-(IntronSpecific$halflife.x)
ThreeUTR_Only_Express<-(ThreeUTRSpecific$halflife.x) 
FiveUTR_Only_Express<-(FiveUTRSpecific$halflife.x) 
Coding_Only_Express<-(CodingSpecific$halflife.x) 
ThreeIntron_Express<-(Intron_and_ThreeUTR$halflife.x) 
All_RNAtargets_Express_0H<-(alldata$halflife.x) 
NonTargets_Express<-(NonTargets$halflife.x)
Multi_Targets_Express<-Multi_Targets$halflife.x
AllmRNAtargets_Express<-(AllmRNAdata$halflife.x)

x<-wilcox.test(Multi_Targets_Express,ThreeIntron_Express)
p.adjust(x$p.value, method= "bonferroni", n=2)

df <- data.frame(x = c(All_RNAtargets_Express_0H, NonTargets_Express,
                       FiveUTR_Only_Express,Coding_Only_Express, 
                       ThreeUTR_Only_Express,Intron_Only_Express, 
                       ThreeIntron_Express, Multi_Targets_Express) 
                 ,ggg = factor(rep(1:8, c(4058,1757,10,22,791, 372,769, 247))))

ggplot(df, aes(x=as.factor(ggg), y=x, fill=ggg)) +
  geom_boxplot(width= .5) +
  theme_minimal() 

ggplot(df, aes(x, colour = ggg)) +
  stat_ecdf()+
  theme_minimal()+
  xlim(2,18) +
  scale_colour_hue(name="my legend", labels=c('AllTranscripts','No_binders','5UTR', 
                                              'coding','3UTR', 'intron', '3UTR_intron','multi'))

ggplot(df, aes(x, colour = ggg, width=.5)) +
  coord_flip() +
  geom_boxplot(width=.5)+
  theme_minimal()

#########################
#Stability based on number of introns and NO 3'UTR contribution 
library(mltools)
intron_only <- alldata %>% filter(IRF3_3UTR == 0 & IRF3_5UTR == 0 & IRF3_Exon == 0)
quantile(intron_only$IRF3_Intron)

intron_only[,"quantiles_Intron"] <- bin_data(intron_only$IRF3_Intron, bins=c(0, 1, 2, 5,65), binType = "explicit")

plot(intron_only$quantiles_Intron)

#HighBinders_intron_Base <-intron_only %>%
  #filter(quantiles_Intron == "[20, 65]" ) ; 

#HighMediumBinders_intron_Base <-intron_only %>%
 # filter(quantiles_Intron == "[10, 65]"); 

MediumBinders_intron_Base <- intron_only %>%
  filter(quantiles_Intron =="[5, 65]")

MediumLowBinders_intron_Base<-intron_only %>%
  filter(quantiles_Intron == "[2, 5)")

LowBinders_intron_Base<-intron_only %>%
  filter(quantiles_Intron == "[1, 2)")

No_intron<- intron_only %>% 
  filter(quantiles_Intron == "[0, 1)")

#HighBinders_intron_Express<-HighBinders_intron_Base$halflife.x
#HighMediumBinders_intron_Express<-HighMediumBinders_intron_Base$halflife.x
MidBinders_intron_Express<-MediumBinders_intron_Base$halflife.x
MidLowBinders_intron_Express<-MediumLowBinders_intron_Base$halflife.x
LowBinders_intron_Express<-LowBinders_intron_Base$halflife.x
NoBinders_intron_Express<-No_intron$halflife.x
AllTranscripts<-alldata$halflife.x

#### TRYING with ggplot for better looking graph

df <- data.frame(x = c(AllTranscripts, 
                       MidBinders_intron_Express, MidLowBinders_intron_Express, 
                       LowBinders_intron_Express, NoBinders_intron_Express
),ggg = factor(rep(1:5, c(4058,60, 119,193, 1757))))


ggplot(intron_only, aes(x = halflife.x, color = quantiles_Intron)) +
  stat_ecdf() +
  xlim(2,18) +
  theme_minimal()

ggplot(intron_only, aes(x= as.factor(IRF3_Intron), y = halflife.x)) + 
  geom_boxplot()

ggplot(intron_only, aes(x= as.factor(quantiles_Intron), y = halflife.x, fill= as.factor(quantiles_Intron))) + 
  geom_boxplot(width= .5) +
  ylim(2,18) +
  theme_minimal()

### testing to see if there is significance between number of introns and stability
x<-wilcox.test(NoBinders_intron_Express, LowBinders_intron_Express, paired = FALSE)
p.adjust(x$p.value, method = "bonferroni", n = 2) 

#### What is going on with ratio of 3'UTR/Intron 
alldata %>% filter(IRF3_3UTR > 0 & IRF3_Intron > 0) %>%
ggplot(aes(y= normalized_Flag_IRF3, x = log(IRF3_3UTR/IRF3_Intron))) +
  geom_point() +
  geom_smooth(method='lm', formula= y~x) +
  theme_minimal()

RIP_par_data<- read.csv("all_counts_RIP_par_genelevel")
IRF3_3UTR_data<- RIP_par_data %>% filter(IRF3_3UTR > 0)

IRF3_low <- alldata %>% filter(IRF3_3UTR >= 1 & IRF3_3UTR <= 2)

IRF3_low[,"quantiles_Intron"] <- bin_data(IRF3_low$IRF3_Intron, bins=c(0, 1, 2, 5, 10, 20, 150), binType = "explicit")

summary(as.data.frame(IRF3_low$quantiles_Intron))

ggplot(IRF3_low, aes(x = halflife.x, color = as.factor(quantiles_Intron))) + 
  stat_ecdf() +
  theme_minimal() 

naive_3UTR_data<- RIP_par_data %>% filter(naive_3UTR > 0)
naive_low <- naive_3UTR_data %>% filter(naive_3UTR >= 2 & naive_3UTR <= 5)

naive_low[,"quantiles_Intron"] <- bin_data(naive_low$naive_Intron, bins=c(0, 1, 2, 5, 10, 20, 450), binType = "explicit")
plot(intron_only$quantiles_Intron)

ggplot(naive_low, aes(x = normalized_Flag_naive, color = as.factor(quantiles_Intron))) + 
  stat_ecdf() +
  theme_minimal() +
  xlim(-4,4)

ggplot(naive_low, aes(y = normalized_Flag_naive, x = as.factor(quantiles_Intron))) + 
  geom_boxplot(width=.75) +
  theme_minimal() +
  ylim(-3,3)

#############################
# What portion out of all the binding sites are in the 3'UTR 
RIP_par_data %>% mutate(total_sites_naive = rowSums(RIP_par_data[7:10]), 
                        total_sites_IRF3= rowSums(RIP_par_data[11:14])) %>%
  ggplot(aes(y= IRF3_3UTR/total_sites_IRF3, x= naive_3UTR/total_sites_naive)) +
  geom_point() +
  theme_minimal() +
  geom_abline(x= 0, b= 1, color = "grey")

IRF3_total<-RIP_par_data %>% 
  mutate(total_number = rowSums(RIP_par_data[11:14])) 

summary(as.data.frame(naive_total$naive_total))

##### 
#plotting enrichment based ofd of a cumulative number of "total" binding sites 
naive_total<- RIP_par_data %>% 
  mutate(total_number = rowSums(RIP_par_data[7:10]))

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

glimpse(alldata %>% mutate(category = as.numeric(IRF3_3UTR/sum(IRF3_5UTR:IRF3_3UTR)))) %>%
ggplot(aes(x = halflife.x, y =category)) + geom_point() +  geom_smooth(method='lm', formula= y~x)

glimpse(IRF3_total)
zero <- naive_total %>% filter(naive_total == "[0, 1)")
test<- naive_total %>% filter(naive_total == "[25, 415]")

x<-wilcox.test(zero$normalized_Flag_naive, test$normalized_Flag_naive)
p.adjust(x$p.value, method = "bonferroni", n = 2) 





#############################################
ThreeUTRSpecific[,"quantiles"] <- bin_data(ThreeUTRSpecific$IRF3_3UTR, bins=c(0, 1, 2, 5, 10, 20, 56), binType = "explicit")
ThreeUTRSpecific <- ThreeUTRSpecific %>% mutate(source = "ThreeSpecific")

Intron_and_ThreeUTR[,"quantiles"] <- bin_data(Intron_and_ThreeUTR$IRF3_3UTR, bins=c(0, 1, 2, 5, 10, 20, 56), binType = "explicit")
Intron_and_ThreeUTR <- Intron_and_ThreeUTR %>% mutate(source = "Three and Intron")

All_3UTRs <- rbind(ThreeUTRSpecific, Intron_and_ThreeUTR)

ggplot(All_3UTRs, aes(x= halflife.x, color = as.factor(source)))+
  stat_ecdf() +
  theme_minimal() +  
  xlim(2,20) +
  facet_wrap(~quantiles) 
  

Intron_and_ThreeUTR[,"quantiles"] <- bin_data(Intron_and_ThreeUTR$IRF3_Intron, bins=c(0, 1, 2, 5, 10, 20, 400), binType = "explicit")

I_Intron_and_ThreeUTR <- Intron_and_ThreeUTR %>% mutate(source = "I_Intron_and_ThreeUTR")

IntronSpecific[, "quantiles"] <- bin_data(IntronSpecific$IRF3_Intron, bins=c( 0, 1, 2, 5, 10, 20, 400), binType = "explicit")
IntronSpecific <- IntronSpecific %>% mutate(source = "IntronSpecific")


All_Intron <- rbind(I_Intron_and_ThreeUTR, IntronSpecific)

ggplot(All_Intron, aes(x= normalized_Flag_IRF3, color = as.factor(source), na.rm = TRUE))+
  stat_ecdf() +
  theme_minimal() + 
  facet_wrap(~quantiles)
