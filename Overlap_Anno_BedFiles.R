############### Overlapping bed files to find clusters and transcripts shared between conditions
library(rtracklayer)
library(GenomicRanges)
library(tidyverse)

#Importing bed files
bed1<- read_tsv("naive_all_cluster.bed", col_names = FALSE)
bed2<-read_tsv("IRF3_all_cluster.bed", col_names = FALSE)

gr1<-import("mRNA_naive_clusterdat.bed")
gr2<-import("mRNA_IRF3_clusterdat.bed")
gr3<-import("naive_all_cluster.bed")
gr4<-import("IRF3_utr3_cluster.bed")
gr5<-import("naive_utr3_cluster.bed")
gr6<-import("db.bed")

#Finding number of overlap clusters via a bed file import
overlap_subset_THP1_HEK<-subsetByOverlaps(gr1, gr2) # can use this to create a new GRanges that can be overlapped with another bed/GRanges file
findOverlaps(gr1, gr2) #outputs number of clusters overlapped


overlap_subset_all3<-subsetByOverlaps(overlap_utr3_subset, gr5)

overlap_intron_df<-as.data.frame(overlap_intron)

#finding the unique clusters
testdiff<-setdiff(gr3, gr5)

#write_rds(overlap_subset_all3, "overlap_subset_all3")
#write_rds(overlap_subset_THP1, "overlap_subset_THP1")

#annotating GRanges output from the bed overlaps
peak_overlap_diff<-annotatePeak(testdiff)
peak_overlap<-annotatePeak(overlap_subset_THP1)

#GO/Reactome Analysis
library(ReactomePA)
peak_overlap_geneID<-as.data.frame(peak_overlap@anno$transcriptId)

pathway1<- enrichPathway(as.data.frame(peak_overlap)$geneId)
pathway2<-enrichPathway(as.data.frame(peak_overlap_all3)$geneId, pvalueCutoff = 0.05, readable=T)
pathway_diff<-enrichPathway(as.data.frame(peak_overlap_diff)$geneId, pvalueCutoff = 0.05, readable=T)

#visualization for CHIPseeker and ReactomePA libraries
barplot(pathway_diff, showCategory = 5)
emapplot(pathway2)

#Importing HEK293 PAR-CLIP data from Mukherjee et al. 2011 Molecular Cell paper 
molcell_2011<-read.csv("molcel_3918_mmc3.csv")
glimpse(molcell_2011)
molcell_2011<-molcell_2011 %>% 
  dplyr::filter(X5UTR > 0 |
                X3UTR > 0 |
                CDS > 0 |
                Intron > 0)

#making quick function 
`%notin%` <- Negate(`%in%`)

#findind the mRNA that are specific bound by ELAVL1 in THP-1s and not HEK293
mRNA_specifictoTHP1s<- naive_mRNA_transcriptbound$gene_short_name  %notin%  IRF3_mRNA_transcriptbound$gene_short_name
naivespecific<-naive_mRNA_transcriptbound[mRNA_specifictoTHP1s,]

mRNA_specifictoTHP1s<-IRF3_mRNA_transcriptbound$gene_short_name  %notin%   naive_mRNA_transcriptbound$gene_short_name
IRF3specific<- IRF3_mRNA_transcriptbound[mRNA_specifictoTHP1s,]


mRNA_THP1s<- naive_mRNA_transcriptbound$gene_short_name  %in%  IRF3_mRNA_transcriptbound$gene_short_name
THP1specific<-naive_mRNA_transcriptbound[mRNA_THP1s,]

mRNA_HEK293<- molcell_2011$Symbol %notin% naive_mRNA_transcriptbound$gene_short_name 
mRNA_HEK293<-molcell_2011[mRNA_HEK293,]

all_mRNA<- molcell_2011$Symbol %in% naive_mRNA_transcriptbound$gene_short_name 
all_mRNA<-molcell_2011[all_mRNA,]


mRNA_specifictoHEK_df<-molcell_2011[mRNA_specifictoTHP1s,]
write_csv(mRNA_specifictoHEK_df, "mRNA_specifictoHEK_df")


#making VennDiagrams of the Overlaps
library(VennDiagram)
draw.pairwise.venn(area1= 7252, area2= 7356, cross.area =5230 )
draw.triple.venn(area1= 151519, area2=133742, area3= 50074, n12=39530, 
                 n23=35607, n13=22850, n123=19404 )

library(UpSetR)
expressionInput <- c(one = 151519, two = 133742, three = 50074, 
                     `one&two` = 39530, `one&three` = 22850,
                     `two&three` = 35607, `one&two&three` = 19404)

upset(fromExpression(expressionInput), order.by = "freq")

################

reactome_THP1s<-read_tsv("ELAVL1_THP1_specificmRNA_reactome.txt")
head(reactome_THP1s)

ggplot(quad1_dat, aes(x = `Pathway name`, y =`#Entities found`, fill = "orange")) +
  geom_bar(stat = 'identity') +
  coord_flip() +
  scale_fill_manual(values=cbPalette) +
  geom_text(aes(label = `Pathway name`, y = 0, hjust = 0), size = 2.5) +
  geom_text(aes(label = signif(`Entities pValue`, digits = 3)  , y = 220, hjust = 0), size = 2.5) +
  theme_minimal() +
  theme(axis.text.x = element_text(size = 12),
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        panel.grid.major.y = element_blank(),
        strip.text.y = element_text(size = 12)) +
  labs(y = "Number of Proteins", title = "GO analysis Enriched and Overexpressed")








