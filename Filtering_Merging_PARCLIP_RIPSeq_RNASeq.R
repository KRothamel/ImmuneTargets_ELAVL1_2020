### This script filters PAR-CLIP data based on RNA-seq expression > 1 FPKM and then catenates RIP, PAR and RNA-seq data for further analysis 
library(tidyverse)
library(dplyr)
library(devtools)
library(reshape2)
library(ggplot2)
library(data.table)
library(mltools)

###Import PARCLIP cluster data per gene. File contains number of binding sites for each transcript
naive_parclipdata<-read.csv("untrimmed1.gene_cl.csv") %>% 
  dplyr::select(GeneName, X5.utr,Intron,Exon, X3.utr)

IRF3_parclipdata<-read.csv("untrimmed2.gene_cl.csv") %>% 
  dplyr::select(GeneName, X5.utr,Intron,Exon, X3.utr)

##renaming gene id columns
colnames(IRF3_parclipdata)[1] <- "gene_short_name"
colnames(naive_parclipdata)[1] <- "gene_short_name"

#combing PAR-CLIP data from both conditions
all_parclip<-full_join(naive_parclipdata, IRF3_parclipdata, by= "gene_short_name") 

#If there is no data for that gene adding "zero" binding sites for that gene 
all_parclip[is.na(all_parclip)] <- 0

##renaming columns
colnames(all_parclip)<- c("gene_short_name",
                          "naive_5UTR", "naive_Intron", 
                          "naive_Exon", "naive_3UTR", 
                          "IRF3_5UTR", "IRF3_Intron", 
                          "IRF3_Exon", "IRF3_3UTR" )
head(all_parclip)

#### Importing  log2 counts data from DeSeq2 Expression data
expressiondata<-read_csv("DeSeq2_log2Counts")
colnames(expressiondata)[8]<-"gene_short_name"
expressiondata<-as.data.frame(expressiondata)

### averaging across replicates for naive and IRF3 samples
naive_counts<- expressiondata %>%  
  dplyr::select(THP1_0H, THP1_0H_2, THP1_0H_3)

IRF3_counts<-expressiondata %>% 
  dplyr::select(THP1_16H, THP1_16H_2, THP1_16H_3)

naive_mean_counts = apply(naive_counts, 1, mean)
IRF3_mean_counts = apply(IRF3_counts, 1, mean)

all_mean_counts<-cbind(naive_mean_counts, IRF3_mean_counts, "gene_short_name" = expressiondata$gene_short_name) 

#### Importing average RIP-Seq data from DeSeq2
RIPdat<-read_csv("averageRIPdata")
colnames(RIPdat)[6]<-"gene_short_name"

normalized_Flag_naive <-RIPdat$Flag.Zero.mean - RIPdat$IgG.Zero.mean
normalized_Flag_IRF3 <-RIPdat$Flag.Sixteen.mean -RIPdat$IgG.Sixteen.mean

all_normalized_RIPdat<-cbind(normalized_Flag_naive,normalized_Flag_IRF3, "gene_short_name"= RIPdat$gene_short_name)

### Joining counts and RIP-Seq data 
all_counts_RIP<-full_join(as.data.frame(all_normalized_RIPdat), as.data.frame(all_mean_counts), by = "gene_short_name")

all_counts_RIP<-na.omit(all_counts_RIP) %>% 
  select(gene_short_name, naive_mean_counts, IRF3_mean_counts, 
         normalized_Flag_naive, normalized_Flag_IRF3)

#### Joining counts_RIP with PAR-CLIP data
all_counts_RIP_par<-full_join(all_counts_RIP, all_parclip, by ="gene_short_name")

#### Replacing only in the PAR-CLIP data NAs with Zeros
all_counts_RIP_par<- all_counts_RIP_par %>% mutate_at(c(6:13), funs(replace(., is.na(.), 0)))

all_counts_RIP_par<- na.omit(all_counts_RIP_par)
glimpse(all_counts_RIP_par)
##########################3
# We are filtering PAR-CLIP hits off of CPM value from RNA-seq. We are filtering the naive PAR-CLIP by the naive RNA-Seq counts and IRF3 PAR-CLIP by the IRF3 RNA-Seq. So we have to separate the data, filter by CPM and then combine and fill in all metadata for the non-targets. 

naive_mRNA_transcriptbound<-all_counts_RIP_par %>% 
  select(gene_short_name, naive_5UTR, naive_Intron, naive_Exon, naive_3UTR, normalized_Flag_naive, naive_mean_counts) %>%
  filter(naive_Exon > 0 | 
           naive_3UTR > 0 | 
           naive_Intron > 0 | 
           naive_5UTR > 0)

# number of mRNA that ELAVL1 binds in specific condition (regardless of RIP-seq score)
naive_mRNA_transcriptbound<- naive_mRNA_transcriptbound %>% 
  dplyr::filter(as.numeric(naive_mean_counts) > 8)

IRF3_mRNA_transcriptbound<-all_counts_RIP_par %>% select(gene_short_name, IRF3_5UTR, IRF3_Intron, IRF3_Exon, IRF3_3UTR, normalized_Flag_IRF3, IRF3_mean_counts) %>%
  filter(IRF3_Exon > 0 | 
           IRF3_3UTR > 0 | 
           IRF3_Intron > 0 | 
           IRF3_5UTR > 0 )

IRF3_mRNA_transcriptbound<- IRF3_mRNA_transcriptbound %>% 
  filter(as.numeric(IRF3_mean_counts) > 8)

all_filter_par<-full_join(naive_mRNA_transcriptbound, IRF3_mRNA_transcriptbound, by = "gene_short_name") %>% 
  select(gene_short_name, naive_5UTR, naive_Intron, naive_Exon, naive_3UTR, IRF3_5UTR, IRF3_Intron, IRF3_Exon, IRF3_3UTR)

all_filter_par[is.na(all_filter_par)] <- 0

all_counts_RIP_PAR_filter<- full_join(all_counts_RIP, all_filter_par, by= "gene_short_name")
all_counts_RIP_PAR_filter<- all_counts_RIP_PAR_filter %>% mutate_at(c(6:13), funs(replace(., is.na(.), 0)))
all_counts_RIP_PAR_filter<- all_counts_RIP_PAR_filter %>% mutate(naive_mean_counts = as.numeric(naive_mean_counts), IRF3_mean_counts = as.numeric(IRF3_mean_counts), normalized_Flag_naive = as.numeric(normalized_Flag_naive), normalized_Flag_IRF3= as.numeric(normalized_Flag_IRF3))

#write_csv(all_counts_RIP_PAR_filter, "all_counts_RIP_par_genelevel1")

#Master table to work with#
master_data<-read_csv("all_counts_RIP_par_genelevel1")

ggplot(master_data, aes(y= (naive_mean_counts-IRF3_mean_counts), x = rank(normalized_Flag_naive) - rank(normalized_Flag_IRF3))) + 
  geom_point(alpha = .1) +
  theme_minimal()
  
  ##########################################################################################################
### Importing specific Clustering Information for Cluster level analysis
clusterdataNaive <- read_csv("untrimmed1.clusters.csv")
clusterdataIRF3<-read_csv("untrimmed2.clusters.csv") %>% as.data.frame()

#changing column names
colnames(clusterdataNaive)[6] <- "gene_short_name"
colnames(clusterdataIRF3)[6] <- "gene_short_name"
#write.csv(clusterdataIRF3, "clusterdataIRF3")
all_naive<- naive_mRNA_transcriptbound %>% 
  dplyr::select(gene_short_name, naive_mean_counts, normalized_Flag_naive, naive_5UTR, naive_Intron, naive_Exon, naive_3UTR)

all_IRF3<-IRF3_mRNA_transcriptbound %>% 
  dplyr::select(gene_short_name, IRF3_mean_counts, normalized_Flag_IRF3, IRF3_5UTR, IRF3_Intron, IRF3_Exon, IRF3_3UTR)

full_cluster_naive<-left_join(all_naive, clusterdataNaive, by = "gene_short_name")
full_cluster_naive<-na.omit(full_cluster_naive)
glimpse(full_cluster_naive)

#write.csv(full_cluster_naive, "full_cluster_naive")
full_cluster_IRF3<-left_join(all_IRF3,clusterdataIRF3, by = "gene_short_name")
full_cluster_IRF3<-na.omit(full_cluster_IRF3)
glimpse(full_cluster_IRF3)

######################################################
#filtering clusters based off of CPM value and Aligned to. 
full_cluster_IRF3$`Aligned to` <- gsub("5'utr", "utr5", full_cluster_IRF3$`Aligned to`)
full_cluster_naive$`Aligned to` <- gsub("5'utr", "utr5", full_cluster_naive$`Aligned to`)
full_cluster_IRF3$`Aligned to` <- gsub("3'utr", "utr3", full_cluster_IRF3$`Aligned to`)
full_cluster_naive$`Aligned to` <- gsub("3'utr", "utr3", full_cluster_naive$`Aligned to`)

#write.csv(full_cluster_naive, "full_cluster_naive")

# let's look at the annotation categories that ELAVL1 binds to
naive_binding_regions_all <- full_cluster_naive %>% 
  mutate(Aligned_to = `Aligned to`) %>% 
  dplyr::count(Aligned_to) %>%
  arrange(desc(n)) %>%
  mutate(freq=n/sum(n))
print(naive_binding_regions_all)

IRF3_binding_regions_all <- full_cluster_IRF3 %>% 
  mutate(Aligned_to = `Aligned to`) %>% 
  dplyr::select(Aligned_to) %>%
  dplyr::count(Aligned_to) %>%
  arrange(desc(n)) %>%
  mutate(freq=n/sum(n))
print(IRF3_binding_regions_all)

# wrangling data for visualization purposes
binding_regions_all<-full_join(naive_binding_regions_all, IRF3_binding_regions_all, by = "Aligned_to") # combines data
binding_regions_all<-as.data.frame(binding_regions_all)
glimpse(binding_regions_all)

#renaming column names for clarity
binding_regions_all<- binding_regions_all %>%
  dplyr::rename("transcript_feature" = Aligned_to, 
                "number_clusters_naive"= n.x,
                "frequency_clusters_naive"= freq.x,
                "number_clusters_IRF3"= n.y, 
                "frequency_clusters_IRF3"= freq.y)
#write.csv(binding_regions_all, "Table_1")

# wrangling data
binding_regions_all <- binding_regions_all %>%
  tidyr::gather(sample_type, frequency, c(3, 5))

# Stacked from ggplot
library(viridis)
library(viridisLite)
library(RColorBrewer)

#visualize the data
ggplot(binding_regions_all, aes(x= sample_type, y= frequency, 
                                fill=transcript_feature)) + 
  geom_bar(position="stack", stat="identity", width= .5)+
  scale_fill_viridis(discrete = T)

### Only want to look at transcript_feature above a certain frequency
major_binding_regions_all<-binding_regions_all %>%
  filter(frequency > .003)

#visualize that
ggplot(major_binding_regions_all, aes(x=sample_type, y=frequency, fill=transcript_feature)) + 
  geom_bar(position="stack", stat="identity", width=.5 ) +
  scale_fill_viridis(discrete = T)

major_binding_regions_all <- major_binding_regions_all %>%
  gather(count_type, counts, c(2, 3)) 

#visualize the number of clusters mapped to each transcript feature
ggplot(major_binding_regions_all, 
       aes(x= transcript_feature, y= counts, fill= count_type)) + geom_col()

#plotting average number of reads based on source
ggplot(Intron_UTR3_counts, aes(x=as.factor(`Aligned to`), y=log10(ReadCount))) +
  geom_boxplot(width=.5, aes(fill= source))

####################################################################
#Making bed files for overlapping of binding-sites
mRNA_naive_clusterdat_bed<-mRNA_naive_clusterdat %>% 
  dplyr::select(Chr, Start, End,Strand) %>%
  write_tsv("mRNA_naive_clusterdat.bed", col_names = FALSE)

mRNA_IRF3_clusterdat_bed<-mRNA_IRF3_clusterdat %>% 
  dplyr::select(Chr, Start, End,Strand) %>%
  write_tsv("mRNA_IRF3_clusterdat.bed", col_names = FALSE)

IRF3_utr3_cluster_bed<- mRNA_IRF3_clusterdat %>% 
  filter(`Aligned to` == "intron") %>%
  dplyr::select(Chr, Start, End,Strand) %>%
  write_tsv("IRF3_utr3_cluster.bed", col_names = FALSE)

naive_utr3_cluster_bed<- mRNA_naive_clusterdat %>% 
  filter(`Aligned to` == "intron") %>%
  dplyr::select(Chr,Start, End,Strand) %>%
  write_tsv("naive_utr3_cluster.bed", col_names = FALSE)

naive_all_cluster_bed<- clusterdataNaive %>% 
  dplyr::select(Chr,Start, End,Strand) %>%
  write_tsv("naive_all_cluster.bed", col_names = FALSE)

IRF3_all_cluster_bed<- clusterdataIRF3 %>% 
  dplyr::select(Chr,Start, End,Strand) %>%
  write_tsv("IRF3_all_cluster.bed", col_names = FALSE)

#IRF3_intron_cluster_bed %>% select(Chr, Start, End, Strand) %>% write_tsv("IRF3_intron_cluster.bed", col_names = FALSE)

# plotting the number of sites for each transcript feature per condition
all_counts_RIP_par %>%
  filter(naive_3UTR > 0 | IRF3_3UTR > 0 | naive_Intron > 0 | IRF3_Intron > 0) %>%
  gather(source, feature, c(naive_3UTR, IRF3_3UTR, naive_Intron, IRF3_Intron)) %>%
  ggplot(aes(y = feature, x= source, color = source)) +
  geom_boxplot(width= .5) +
  theme_minimal() 

all_counts_RIP_par %>%
  filter(naive_Intron > 0 | IRF3_Intron > 0) %>%
  gather(source, intron, c(naive_Intron, IRF3_Intron)) %>%
  ggplot(aes(y = log2(intron), x= source,  fill = source)) +
  geom_boxplot(width = .5) +
  theme_minimal()

