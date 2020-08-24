# This script pulls out annotated introns and 3'UTR for every gene. We also then calculate the width of introns and 3'UTRs as well as the AU frequency and the ratio of the number of (3'UTR binding sites/intronic binding sites)
library(RCAS)
library(ChIPseeker)
library(clusterProfiler)
library(tidyverse)
library(Biostrings)
library(biomaRt)
library(GenomicFeatures)
library(rtracklayer)

#Importing gtf file and converting to GRanges and pulling out transcript genomic coordinates 
gff <- import.gff("gencode.v19.chr_patch_hapl_scaff.annotation.gtf")
txdb <- GenomicFeatures::makeTxDbFromGRanges(gff)
txdbFeatures <- getTxdbFeaturesFromGRanges(gff)
transcriptCoords <- GenomicFeatures::transcripts(txdb)

#From Genomic features-> pulling out genomic coordinate of all introns defined in gtf
introns <- intronsByTranscript(txdb)
print(introns@unlistData)

#From Genomic features-> pulling out genomic coordinate of all utr3 defined in gtf
utr3  <- threeUTRsByTranscript(txdb)
print(utr3@unlistData)

# From RCAS, annotating the coordinates. Also where you get exon ID for ReactomePA GO visualizations later
peakAnno_Introns <- annotatePeak(introns@unlistData, TxDb=txdb)
peakAnno_utr3 <- annotatePeak(utr3@unlistData, TxDb=txdb)

peakAnno_Introns <- as.data.frame(peakAnno_Introns)
peakAnno_utr3 <- as.data.frame(peakAnno_utr3)
glimpse(peakAnno_Introns)

#removing extra decimal point in the ensembl gene id 
peakAnno_Introns <- peakAnno_Introns %>% 
  mutate(short_ensembl_id = str_remove(geneId, "\\.[0-9]*"))
peakAnno_utr3 <- peakAnno_utr3 %>% 
  mutate(short_ensembl_id = str_remove(geneId, "\\.[0-9]*"))

#now getting ensembl and hgnc_gene_id 
ensmart=useMart("ensembl", dataset="hsapiens_gene_ensembl")
#genes <- getBM( attributes=c("hgnc_symbol","ensembl_gene_id","chromosome_name","start_position","end_position"), mart = ensmart)

#saveRDS(genes, "ensembl_genes")
genes<-readRDS("ensembl_genes")
colnames(genes)[2]<- "short_ensembl_id"
genes<- genes %>% 
  dplyr::select(hgnc_symbol, short_ensembl_id)

full_Anno_Introns<- left_join(peakAnno_Introns, genes, by = "short_ensembl_id")
full_Anno_Introns <- na.omit(full_Anno_Introns) 

full_Anno_utr3<- left_join(peakAnno_utr3, genes, by = "short_ensembl_id")
full_Anno_utr3 <- na.omit(full_Anno_utr3) 

#rename columns for clarity
full_Anno_utr3 <- full_Anno_utr3 %>% 
 dplyr::rename(utr3.start = start, utr3.end = end, utr3.width = width, gene_short_name = hgnc_symbol)

full_Anno_Introns <- full_Anno_Introns %>% 
  dplyr::rename(intron.start = start, intron.end = end, intron.width = width, gene_short_name = hgnc_symbol)

full_Anno_Introns<- full_Anno_Introns %>% dplyr::select(intron.start, intron.end, intron.width, short_ensembl_id, gene_short_name)

glimpse(full_Anno_Introns)

#adding intron and 3UTR length to respective tables
intron_totals <- full_Anno_Introns %>%
  group_by(short_ensembl_id) %>%
  mutate(intron_full = sum(intron.width, na.rm= TRUE))

intron_totals<-intron_totals[!duplicated(intron_totals[c("short_ensembl_id")]),]
intron_totals <- intron_totals %>% 
  dplyr::select(short_ensembl_id, gene_short_name, intron_full)

# loaded in PAR_RIP_Expression data and overlapping with annotated 3'UTR metadata
all_PARCLIP_counts<-read_csv("all_counts_RIP_par_genelevel1")

PARCLIP_3UTR_annot<- left_join(all_PARCLIP_counts, full_Anno_utr3, by = "gene_short_name")
PARCLIP_3UTR_annot <- na.omit(PARCLIP_3UTR_annot)
PARCLIP_3UTR_annot <- PARCLIP_3UTR_annot %>% 
  dplyr::select(-X1, -X1_1)

PAR_3UTR_Intron_annt<-left_join(PARCLIP_3UTR_annot, intron_totals, by="short_ensembl_id")
PAR_3UTR_Intron_annt<- PAR_3UTR_Intron_annt[!duplicated(PAR_3UTR_Intron_annt[c('short_ensembl_id', 'exon_id')]),]
glimpse(PAR_3UTR_Intron_annt)

# TRYING TO add actual full 3'UTR sequence to each gene
#ensemblid = getBM("ensembl_gene_id", mart=ensmart)
#utr3_bioMART_seq = getSequence(id = ensemblid[1:67140,1], type="ensembl_gene_id", seqType="3utr", mart=ensmart)
utr3_bioMART_seq<-read_csv("bioMart_3utr_ensembl")
colnames(utr3_bioMART_seq)[3]<- "short_ensembl_id"

# combining 3'UTR with hgnc gene names
annotated_utr3<- full_join(utr3_bioMART_seq , PAR_3UTR_Intron_annt, by = "short_ensembl_id")
glimpse(annotated_utr3)

annotated_utr3<- annotated_utr3 %>% 
  dplyr::select(-X1,-exon_name, -exon_rank, -annotation) %>%
  na.omit()

#write.csv(annotated_utr3, "annotated_utr3")

#getting rid of all 3UTR that don't have a sequence
annotated_utr3<-annotated_utr3 %>% filter(`3utr` != "Sequence unavailable")
annotated_utr3<- annotated_utr3[!duplicated(annotated_utr3[c('short_ensembl_id', 'exon_id')]),]
annotated_utr3<-left_join(annotated_utr3, genes, by="short_ensembl_id")
annotated_utr3<-annotated_utr3[!is.na(annotated_utr3$gene_short_name.x),]
glimpse(annotated_utr3)

#finding the length of each 3'UTR
annotated_utr3<-annotated_utr3 %>%
  mutate(length_3UTR = width(annotated_utr3$`3utr`))

#Use Biostrings to count the AU frequency of the 3'UTR sequences
library(Biostrings)
utr3_string<-DNAStringSet(pull(annotated_utr3, `3utr` ))

AU_frequency<-letterFrequency(utr3_string, letters= c("T", "A"))

AU_frequency<-AU_frequency %>% 
  as.data.frame(AU_frequency) %>%
  mutate(AU = rowSums(AU_frequency))

# counting the literal "A" & "U" frequency and dividing it by the length of the 3'UTR
annotated_utr3<- annotated_utr3 %>%
  mutate(AUfrequence = (AU_frequency$AU/length_3UTR))

#comparing the number of binding sites in each 3UTR that ELAVL1 binds
ggplot(annotated_utr3, aes(x=naive_3UTR, y= IRF3_3UTR)) +
  geom_jitter(alpha = 0.3, color= "steelblue") +
  geom_abline(slope = 1, intercept = 0) +
  theme_minimal()

#plotting 3'UTR length of ELAVL1 binds sites vs. all transcripts 
annotated_utr3 %>%
  mutate(Category = case_when(
    IRF3_3UTR > 1 & naive_3UTR > 1 ~ "Both", 
    naive_3UTR > 1 ~ "Naive only",
    IRF3_3UTR > 1 ~ "IRF3 only",
    TRUE ~ "non-binders"
  )) %>% 
  ggplot(aes(x= log10(length_3UTR))) + 
  geom_density(aes(color = Category, fill=Category, alpha=.1)) +
  geom_density(color = "black") +
  theme_minimal() +
  labs(y="density", 
       x="log(3UTR length)", 
       title="Length of 3UTR of targets of ELAVL1")

#plotting AU frequency of ELAVL1 binds sites vs. all transcripts 
annotated_utr3 %>%
  mutate(Category = case_when(
    IRF3_3UTR > 1 & naive_3UTR > 1~ "> 1 binding site in both", 
    naive_3UTR > 1 ~ "> 1 binding site naive only",
    IRF3_3UTR > 1 ~ "> 1 binding site IRF3 only",
    TRUE ~ "non-targets"
  )) %>% 
  ggplot(aes(x= AUfrequence)) + 
  geom_density(aes(alpha = 0.1, color = Category, fill= Category)) +
  geom_density(color = "black", fill= "gray", alpha = .1) +
  theme_minimal() +
  labs(y="density", 
       x="AU frequency", 
       title="AU frequency of the 3UTRs of mRNA targets")

#plotting  full intron length for all targets based on condition 
annotated_utr3 %>%
  mutate(Category = case_when(
    IRF3_3UTR > 0 & naive_3UTR > 0 ~ "bound in both conditions", 
    naive_Intron > 0 ~ "naive only",
    IRF3_Intron > 0 ~ "IRF3 only",
    TRUE ~ "non-binders"
  )) %>% 
  ggplot(aes(x= log10(intron_full))) + 
  geom_density(aes(color = Category, alpha=.01)) +
  geom_density(alpha=.1, color = "black") +
  theme_minimal() + 
  labs(y="density", 
       x="Total intron space", 
       title="Intron length of the targets of ELAVL1")

#plot the ratio utr3 length to intron length
annotated_utr3 <- annotated_utr3 %>%
  mutate(utr3tointron = length_3UTR/intron_full)

#plotting 3UTR/intron ratio based on expression fold-change
annotated_utr3 %>%
  mutate(Category = case_when(
    naive_mean_counts & IRF3_mean_counts > 5  ~ "expressed in both ",
    naive_mean_counts > 5 ~ "naive expressed only",
    IRF3_mean_counts > 5  ~ "IRF3 expressed only",
    TRUE ~ "< 0 in both"
  )) %>% 
  ggplot(aes(x= (log(utr3tointron)))) + 
  geom_density(aes(color = Category)) +
  geom_density(alpha=.2, color = "black") +
  theme_minimal() + 
  labs(y="density", 
       x="3UTR/Intron length ratio", 
       title="Ratio of 3UTR/Intron length of the targets of ELAVL1")

#plotting 3UTR/intron ratio based on condition value
annotated_utr3 %>%
  mutate(Category = case_when(
    IRF3_3UTR > 1 & naive_3UTR > 1 ~ "bound in 3UTR in both", 
    IRF3_3UTR > 1 ~ "naive only",
    naive_3UTR > 1 ~ "IRF3 only",
    TRUE ~ "non-3UTR binders", 
  )) %>% 
  ggplot(aes(x= (log10(utr3tointron)))) + 
  geom_density(aes(color = Category)) +
  geom_density(alpha=.2, color = "black") +
  theme_minimal() +
  xlim(-5, 2) +
  labs(y="density", 
       x="3UTR/Intron length ratio", 
       title="Ratio of 3UTR/Intron length of the targets of ELAVL1")

#comparing 3utr/intron length to # of binding sites
ggplot(annotated_utr3, aes(x= as.factor(IRF3_3UTR), y = log(utr3tointron))) + 
  geom_boxplot(width= .5, color="steelblue") 
  #geom_boxplot(aes(x=as.factor(naive_3UTR), y= log(AUfrequence)), color= "maroon")
  #geom_violin(scale="area", fill="orange", alpha=.5)


#plotting 3UTR/intron ratio against number of binding sites
library(mltools)
annotated_utr3[,"UTR3_IRF3_quantiles"] <- bin_data(annotated_utr3$IRF3_3UTR, bins=c(0, 1, 2, 5, 10, 20, 30, 40,56), binType = "explicit")
annotated_utr3[,"UTR3_naive_quantiles"] <- bin_data(annotated_utr3$naive_3UTR, bins=c(0, 1, 2, 5, 10, 20,30, 40, 56), binType = "explicit")
### Look at the breakdown of the data based on number of 3'utr binding sites
plot(annotated_utr3$UTR3_naive_quantiles)

ggplot(annotated_utr3, aes(x= as.factor(UTR3_naive_quantiles), y= log(utr3tointron)), na.rm=TRUE)+
  geom_boxplot(width=.5, fill="blue", alpha=.5) +
  theme_minimal() +
  labs(y="log(3UTR to Intron length ratio)", 
       x="Number of naive 3'UTR binding-sites", 
       title="Ratio of 3UTR/Intron length based off number of naive binding sites")

###############################################
naive_urt3_bound<- annotated_utr3 %>% 
  filter(annotated_utr3$naive_3UTR > 0 )

IRF3_urt3_bound<- annotated_utr3 %>% 
  filter(annotated_utr3$IRF3_3UTR > 0)

`%notin%` <- Negate(`%in%`)

IRF3_3UTR_speicificTranscripts<-IRF3_urt3_bound%>%
  filter(as.vector(IRF3_urt3_bound$gene_short_name.x) %notin% as.vector(naive_urt3_bound$gene_short_name.x))


###########################################
#playing around with annotating the peak overlaps

###################
annotated_utr3 %>%
  mutate(Category = case_when(
    (IRF3_mean_counts/naive_mean_counts) > 1 ~ "upregulated", 
    (naive_mean_counts/IRF3_mean_counts) > 1 ~ "downregulated",
  )) %>% 
  ggplot(aes(x= (log(intron_full)))) + 
  geom_density(aes(color = Category)) +
  geom_density(alpha=.2, color = "black") +
  theme_minimal() +
  labs(y="density", 
       x="3'UTR length", 
       title="3'UTR length of upregulated and downregulated genes")


#ANOTHER WAY of calculate motif frequency, calculating frequency of motif "TTTTT"
x<-vcountPattern("TTTTT",utr3_string,max.mismatch = 1)
summary(x)
test<-cbind(annotated_utr3, x)

# seeing if more 3'UTR sites correlated with number of motifs
ggplot(test, aes(IRF3_3UTR, x)) + 
  geom_point(alpha = .2, color = "black")+
  theme_minimal()

# are there more motifs in the upregulated genes
test %>% 
  mutate(Category = case_when(
    IRF3_mean_counts/naive_mean_counts > 1  ~ "increased levels",
    naive_mean_counts/IRF3_mean_counts > 1  ~ "decreased levels"
  )) %>%
  ggplot(aes(x= log2(x), color = Category)) + 
  geom_density(aes(color = Category, alpha=.1)) +
  geom_density(color = "black") +
  theme_minimal() 
