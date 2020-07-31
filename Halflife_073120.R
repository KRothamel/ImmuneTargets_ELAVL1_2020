library(RCAS)
library(ChIPseeker)
library(clusterProfiler)
library(tidyverse)
library(Biostrings)
library(biomaRt)
library(GenomicFeatures)
library(rtracklayer)
library(dplyr)

############# annotate bed file. 
gff <- import.gff("gencode.v19.chr_patch_hapl_scaff.annotation.gtf") #gtf from PAR-CLIP accessories
txdb <- GenomicFeatures::makeTxDbFromGRanges(gff)
txdbFeatures <- getTxdbFeaturesFromGRanges(gff)
transcriptCoords <- GenomicFeatures::transcripts(txdb)

#Importing the output file from alleyoop merge (SLAM-DUNK output). File contains all T-to-C fractions per 3'UTR per sample. 
data<-read_tsv("/Volumes/BackupPlus/alleyoop/merge_all", skip=3)
glimpse(data)
bed<-data %>% 
  dplyr::select(Chromosome, Start, End, Strand)
#write_tsv(bed, "SLAM.bed")

# From RCAS, annotating the SLAM-Seq coordinates. Also where you get exon ID for ReactomePA GO visualizations later
library(org.Hs.eg.db)
mybed<-import("SLAM.bed")
peakAnno_SLAM <- annotatePeak(mybed, TxDb = txdb, annoDb = "org.Hs.eg.db")
peakAnno_SLAM_df <- as.data.frame(peakAnno_SLAM)
dplyr::glimpse(peakAnno_SLAM_df)

my_anno<-peakAnno_SLAM_df %>% 
  dplyr::select(Chromosome= seqnames, Start= start, End= end, geneId, ENTREZID, SYMBOL,GENENAME)

dplyr::glimpse(my_anno)

anno_SLAM<-left_join(data, my_anno, by="End" )
anno_SLAM<-na.omit(anno_SLAM)

full_Anno_SLAM<- anno_SLAM[!duplicated(anno_SLAM[c('SYMBOL')]),]
#write_csv(full_Anno_SLAM, "full_Anno_SLAM")

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
result_TA_halflife<- read_csv("result_TA_halflife")
result_KO_halflife<- read_csv("result_KO_halflife")

PARdata<-read_csv("all_counts_RIP_par_genelevel1")
glimpse(PARdata)
colnames(PARdata)[1]<-"SYMBOL" 
PARdata <- PARdata %>% mutate("SYMBOL" = as.character(SYMBOL))

all_halflife<- full_join(result_TA_halflife, result_KO_halflife, by = "SYMBOL")
alldata<-full_join(PARdata, all_halflife, by = "SYMBOL")
alldata<-na.omit(alldata)
glimpse(alldata)

#plotting CDF plot with stability and transcript feature
IntronSpecific<-filter(alldata, alldata$IRF3_3UTR == 0 & alldata$IRF3_Intron > 0 & alldata$IRF3_5UTR == 0 & alldata$IRF3_Exon == 0)

ThreeUTRSpecific<-filter(alldata, alldata$IRF3_Intron == 0 & alldata$IRF3_3UTR !=0 & alldata$IRF3_5UTR==0 & alldata$IRF3_Exon==0); 
dim(ThreeUTRSpecific)

FiveUTRSpecific<- filter(alldata, alldata$IRF3_3UTR == 0 & alldata$IRF3_Intron == 0 & alldata$IRF3_Exon == 0 & alldata$IRF3_5UTR !=0); 
dim(FiveUTRSpecific)

CodingSpecific<- filter(alldata, alldata$IRF3_3UTR == 0 & alldata$IRF3_Intron == 0 & alldata$IRF3_5UTR ==0 & alldata$IRF3_Exon> 0); 
dim(CodingSpecific)

Intron_and_ThreeUTR<- alldata %>% 
  filter(alldata$IRF3_Intron> 0 & alldata$IRF3_3UTR !=0 & alldata$IRF3_5UTR==0 & alldata$IRF3_Exon==0 ); dim(Intron_and_ThreeUTR)

NonTargets<-filter(alldata, alldata$IRF3_3UTR == 0 & alldata$IRF3_Intron == 0 & alldata$IRF3_5UTR== 0 & alldata$IRF3_Exon== 0); 
dim(NonTargets)

Multi_Targets<-filter(alldata, alldata$IRF3_3UTR > 0 & alldata$IRF3_Intron > 0  & alldata$IRF3_5UTR> 0 | alldata$IRF3_Exon > 0); 
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

x<-wilcox.test(Multi_Targets_Express,NonTargets_Express)
p.adjust(x$p.value, method= "bonferroni", n= 335)

df <- data.frame(x = c(All_RNAtargets_Express_0H, NonTargets_Express,
                       FiveUTR_Only_Express,Coding_Only_Express, 
                       ThreeUTR_Only_Express,Intron_Only_Express, 
                       ThreeIntron_Express, Multi_Targets_Express) 
                 ,ggg = factor(rep(1:8, c(4058,1918,9,18,703,330,747, 335))))

ggplot(df, aes(x, colour = ggg)) +
  stat_ecdf()+
  theme_minimal()+
  xlim(2,18) +
  scale_colour_hue(name="my legend", labels=c('AllTranscripts','No_binders','5UTR', 
                                              'coding','3UTR', 'intron', '3UTR_intron','multi'))


####################################
library(mltools)

intron_only<-filter(alldata, alldata$IRF3_3UTR == 0 & alldata$IRF3_Intron > 0 & alldata$IRF3_5UTR == 0 & alldata$IRF3_Exon == 0)

intron_only[,"quantiles_Intron"] <- bin_data(intron_only$IRF3_Intron, bins=c(0, 1, 2, 5,65), binType = "explicit")

plot(intron_only$quantiles_Intron)

MediumBinders_intron_Base <- intron_only %>%
  filter(quantiles_Intron =="[5, 65]")

MediumLowBinders_intron_Base<-intron_only %>%
  filter(quantiles_Intron == "[2, 5)")

LowBinders_intron_Base<-intron_only %>%
  filter(quantiles_Intron == "[1, 2)")


#HighBinders_intron_Express<-HighBinders_intron_Base$halflife.x
#HighMediumBinders_intron_Express<-HighMediumBinders_intron_Base$halflife.x
MidBinders_intron_Express<-MediumBinders_intron_Base$halflife.x
MidLowBinders_intron_Express<-MediumLowBinders_intron_Base$halflife.x
LowBinders_intron_Express<-LowBinders_intron_Base$halflife.x
NonTargets_Express<-(NonTargets$halflife.x)
AllTranscripts<-alldata$halflife.x

zero_introns<- alldata %>% filter(IRF3_Intron == 0 & IRF3_3UTR > 0 | IRF3_5UTR > 0 | IRF3_Exon > 0) 
zero_introns_df <- zero_introns$halflife.x

x<-wilcox.test(zero_introns_df, NonTargets_Express)
p.adjust(x$p.value, method = "bonferroni", n = 1603) 

#### TRYING with ggplot for better looking graph

df <- data.frame(x = c(AllTranscripts, 
                       MidBinders_intron_Express, MidLowBinders_intron_Express, 
                       LowBinders_intron_Express, NonTargets_Express, zero_introns_df
),ggg = factor(rep(1:6, c(4058,56, 100,174, 1918, 1063))))


ggplot(df, aes(x, colour = ggg)) +
  stat_ecdf()+
  theme_minimal()+
  xlim(2,18) +
  scale_colour_hue(name="my legend", labels=c('AllTranscripts','Mid', 'Mid-Low', 'Low', 'Nontargets', 'Zero_Express'))


##############################################################################
#Making CDFs
library(mltools)

threeutr_only<-alldata %>%filter(IRF3_3UTR > 0 & IRF3_5UTR==0 & IRF3_Exon==0 & IRF3_Intron == 0); 
threeutr_only[,"quantiles_IRF3"] <- bin_data(threeutr_only$IRF3_3UTR, bins=c(0, 1, 2, 5, 50), binType = "explicit")

#HighBinders_3UTR_Base <-threeutr_only %>%
#filter(quantiles_IRF3 == "[20, 56]" ) ; 

#HighMediumBinders_3UTR_Base <-threeutr_only %>%
  #filter(quantiles_IRF3 == "[10, 50]"); 

MediumBinders_3UTR_Base <- threeutr_only %>%
  filter(quantiles_IRF3 =="[5, 50]")

MediumLowBinders_3UTR_Base<-threeutr_only %>%
  filter(quantiles_IRF3 == "[2, 5)")

LowBinders_3UTR_Base<-threeutr_only %>%
  filter(quantiles_IRF3 == "[1, 2)")

###################################
#HighBinders_3UTR_Express<-HighBinders_3UTR_Base$halflife.x
#HighMediumBinders_3UTR_Express<-HighMediumBinders_3UTR_Base$halflife.x
MidBinders_3UTR_Express<-MediumBinders_3UTR_Base$halflife.x
MidLowBinders_3UTR_Express<-MediumLowBinders_3UTR_Base$halflife.x
LowBinders_3UTR_Express<-LowBinders_3UTR_Base$halflife.x
NonTargets_Express<-(NonTargets$halflife.x)
AllTranscripts<-alldata$halflife.x

non3URtargets<- alldata %>% filter(IRF3_3UTR == 0 ) %>% filter( IRF3_Intron > 0 | IRF3_5UTR > 0 | IRF3_Exon > 0) 
non3URtargets_df <- non3URtargets$halflife.x

x<-wilcox.test(non3URtargets_df, NonTargets_Express)
p.adjust(x$p.value, method = "bonferroni", n = 690) 


#### TRYING with ggplot for better looking graph
library("reshape2")
library("plyr")
library("ggplot2")

#df <- data.frame(x = c(AllTranscripts,
                       #MidBinders_3UTR_Express, 
                       #MidLowBinders_3UTR_Express, LowBinders_3UTR_Express, 
                       #NonTargets_Express, non3URtargets_df),
                 #ggg = factor(rep(1:5, c(456,718,568,1918,690))))

#ggplot(df, aes(x, colour = ggg)) +
 # stat_ecdf() +
  #theme_minimal() +
  #scale_colour_hue(name="my legend", labels=c( 'Medium (5-10)',
                                              'Medium_Low (2-5)', 'Low (1-2)','No Binders', 'Non-3UTR targets'))


df <- data.frame(x = c(AllTranscripts,
  MidBinders_3UTR_Express, 
  MidLowBinders_3UTR_Express, LowBinders_3UTR_Express, 
  NonTargets_Express, non3URtargets_df),
  ggg = factor(rep(1:6, c(4058,101,273,329,1918,398))))

ggplot(df, aes(x, colour = ggg)) +
  stat_ecdf() +
  theme_minimal() +
  xlim(2,18) +
  scale_colour_hue(name="my legend", labels=c('All', 'Medium (5-10)',
                                              'Medium_Low (2-5)', 'Low (1-2)','No Binders', 'Non-3UTR targets'))





