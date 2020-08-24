#GO, Reactome and KEGG analysis from bed files and for transcripts specific for each PAR-CLIP and RIP-seq condition 
library(tidyverse) # data-warngling
library(reshape2) # 
library(ggplot2) # 
library(ggthemes) #
library(scales) # 
library(GenomicFeatures)
library(rtracklayer)
library(AnnotationDbi)
library(Biostrings)
library(seqLogo)
library(BSgenome.Hsapiens.UCSC.hg19)
library(DT)
library(Hmisc)
library(RColorBrewer)

#visualizing GO terms
library(ReactomePA)
library(org.Hs.eg.db)
library(clusterProfiler)

#getting ENTRREZID for GO visualization from clusterprofiler and ReactomePA
mRNA_specifictoTHP1s_df<-read_csv("mRNA_specifictoTHP1s_df")

gene.df <- bitr(mRNA_specifictoTHP1s_df$gene_short_name, fromType = "SYMBOL",
                toType = c("ENTREZID", "ENSEMBL"),
                OrgDb = org.Hs.eg.db)

#Running GO ontology from clusterprofiler
ego <- enrichGO(gene = gene.df$ENTREZID,
                OrgDb = org.Hs.eg.db,
                ont  = "ALL",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.01,
                qvalueCutoff  = 0.05)

#Running KEGG pathway analysis from clusterprofiler
ego_KEGG<- enrichKEGG(gene.df$ENTREZID, organism = "hsa",  pvalueCutoff  = 0.01, 
                     qvalueCutoff  = 0.05, pAdjustMethod = "BH")

#Running REACTOME analysis 
ego_RA<-ReactomePA::enrichPathway(gene.df$ENTREZID, organism = "human", readable = TRUE)

#Make a bargrpah of the GO ontology analysis
barplot(ego_RA) + theme_minimal() +
  coord_flip() + 
  theme(axis.text.x = element_text(size=12), 
        axis.text.y = element_text(size=8),
        panel.grid.major = element_blank(),
        strip.text.y = element_text(size=12)) +
  labs(y= "Number of transcripts", title = "REACTOME terms for THP1 specific mRNAs")

cnetplot(ego_RA, showCategory = 3, layout= "kk", label= node, circular=TRUE, colorEdge=TRUE)

#Importing genes that are specific to HEK293T PAR-CLIP 
mRNA_specifictoHEK_df<-read_csv("mRNA_specifictoHEK_df")

HEK.df<-bitr(mRNA_specifictoHEK_df$Symbol, fromType = "SYMBOL",
             toType = c("ENTREZID", "ENSEMBL"),
             OrgDb = org.Hs.eg.db)

HEKgo <- enrichGO(gene = HEK.df$ENTREZID,
                OrgDb = org.Hs.eg.db,
                ont  = "ALL",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.01,
                qvalueCutoff  = 0.05)

HEK_KEGG<-enrichKEGG(HEK.df$ENTREZID, organism = "hsa",  pvalueCutoff  = 0.01, 
                     qvalueCutoff  = 0.05, pAdjustMethod = "BH")

HEK_RA<-enrichPathway(HEK.df$ENTREZID, organism = "human", 
                      pAdjustMethod = "BH", readable = TRUE)


clusterProfiler::dotplot(HEK_KEGG)
emapplot(HEK_RA, showCategory = 3)

HEK_RA<-enrichPathway(HEK.df$ENTREZID, organism = "human", readable = TRUE)

barplot(HEKgo) + theme_minimal() +
  coord_flip() + 
  theme(axis.text.x = element_text(size=12), 
        axis.text.y = element_text(),
        panel.grid.major = element_blank(),
        strip.text.y = element_text(size=12)) +
  labs(y= "Number of transcripts", title = "GO analysis HEK293 specific mRNAs")

barplot(ego_RA)

## Create enrichGO and pathway analysis objects
library(clusterProfiler)
library(org.Hs.eg.db)
library(ReactomePA)
library(Rgraphviz)
library(cowplot)
library(ggrepel)

all_naive<- read_csv("naive_mRNA_transcriptbound")
all_IRF3<-read_csv("IRF3_mRNA_transcriptbound")

`%notin%` = Negate(`%in%`)

#Only want the enriched targets 
naive_x<- all_naive$gene_short_name %notin% all_IRF3$gene_short_name
naive_spec<-all_naive[naive_x,]

IRF3_x<- all_IRF3$gene_short_name %notin% all_naive$gene_short_name
IRF3_spec<-all_IRF3[IRF3_x,]

shared_x<- all_naive$gene_short_name %in% all_IRF3$gene_short_name
shared_spec<-all_naive[shared_x,]

# GO visualizations of PAR-CLIP and RIP specific data 
naive_mRNA<- spec_naive_enrich$gene_short_name
naive_mRNA_df <- bitr(naive_mRNA, fromType = "SYMBOL",
                      toType = c("ENTREZID", "ENSEMBL"),
                      OrgDb = org.Hs.eg.db)

naive_RA<-enrichPathway(naive_mRNA_df$ENTREZID, organism = "human", readable = TRUE, pvalueCutoff  = 0.1)

naive_mRNA_BP <- enrichGO(gene     = naive_mRNA_df$ENTREZID,
                          OrgDb         = org.Hs.eg.db,
                          ont           = "BP",
                          pAdjustMethod = "BH",
                          pvalueCutoff  = 0.05,
                          qvalueCutoff  = 0.05,
                          readable = TRUE)

barplot(naive_RA) + theme_minimal() +
  coord_flip() + 
  theme(axis.text.x = element_text(size=12), 
        axis.text.y = element_text(),
        panel.grid.major = element_blank(),
        strip.text.y = element_text(size=12)) +
  labs(y= "Number of transcripts", title = "GO analysis naive specific")


IRF3_mRNA <- spec_IRF3_enrich$gene_short_name
IRF3_mRNA_df <- bitr(IRF3_mRNA, fromType = "SYMBOL",
                toType = c("ENTREZID", "ENSEMBL"),
                OrgDb = org.Hs.eg.db)

IRF3_BP <- enrichGO(gene     = IRF3_mRNA_df$ENTREZID,
                   OrgDb         = org.Hs.eg.db,
                   ont           = "BP",
                   pAdjustMethod = "BH",
                   pvalueCutoff  = 0.01,
                   qvalueCutoff  = 0.05,
                   readable = TRUE)

IRF3_RA<-enrichPathway(IRF3_mRNA_df$ENTREZID, organism = "human", readable = TRUE, pvalueCutoff  = 0.1)

naive_RA_df<- as.data.frame(naive_RA)
write_csv(naive_RA_df, "SupFigure2_naive_specific_ReactomeAnalysis")

barplot(IRF3_RA) + theme_minimal() +
  coord_flip() + 
  theme(axis.text.x = element_text(size=12), 
        axis.text.y = element_text(),
        panel.grid.major = element_blank(),
        strip.text.y = element_text(size=12)) +
  labs(y= "Number of transcripts", title = "GO analysis IRF3 specific")


shared_mRNA <- shared_enriched$gene_short_name
shared_mRNA_df <- bitr(shared_mRNA, fromType = "SYMBOL",
                     toType = c("ENTREZID", "ENSEMBL"),
                     OrgDb = org.Hs.eg.db)
shared_BP <- enrichGO(gene     = shared_mRNA_df$ENTREZID,
                    OrgDb         = org.Hs.eg.db,
                    ont           = "BP",
                    pAdjustMethod = "BH",
                    pvalueCutoff  = 0.01,
                    qvalueCutoff  = 0.05,
                    readable = TRUE)

shared_RA<-enrichPathway(shared_mRNA_df$ENTREZID, organism = "human", readable = TRUE, pvalueCutoff  = 0.1)

barplot(shared_RA) + theme_minimal() +
  coord_flip() + 
  theme(axis.text.x = element_text(size=12), 
        axis.text.y = element_text(),
        panel.grid.major = element_blank(),
        strip.text.y = element_text(size=12)) +
  labs(y= "Number of transcripts", title = "GO analysis IRF3 specific")

shared_RA_df<- as.data.frame(shared_RA)
#write_csv(shared_RA_df, "SupFigure2_shared_enriched_ReactomeAnalysis")
#############################################################################
functional_targets<- read_csv("functional_targets")
func_mRNA <- spec_IRF3_enrich_static$gene_short_name

func_mRNA_df <- bitr(functional_targets$gene, fromType = "SYMBOL",
                       toType = c("ENTREZID", "ENSEMBL"),
                       OrgDb = org.Hs.eg.db)

func_mRNA_BP <- enrichGO(gene     = func_mRNA_df$ENTREZID,
                      OrgDb         = org.Hs.eg.db,
                      ont           = "BP",
                      pAdjustMethod = "BH",
                      pvalueCutoff  = 0.1,
                      qvalueCutoff  = 0.5,
                      readable = TRUE)

shared_RA<-enrichPathway(func_mRNA_df$ENTREZID, organism = "human", readable = TRUE)

barplot(shared_RA) + theme_minimal() +
  coord_flip() + 
  theme(axis.text.x = element_text(size=12), 
        axis.text.y = element_text(),
        panel.grid.major = element_blank(),
        strip.text.y = element_text(size=12)) +
  labs(y= "Number of transcripts", title = "GO analysis IRF3 specific")


func_mRNA_df_BP <- enrichKEGG(gene     = func_mRNA_df$ENTREZID,
                             #OrgDb         = org.Hs.eg.db,
                             #ont           = "BP",
                             pAdjustMethod = "BH",
                             pvalueCutoff  = 0.1,
                             qvalueCutoff  = 0.5)
#readable = TRUE)












#Reactome PA enrichment 
IRF3_RA<-enrichPathway(IRF3_mRNA_df$ENTREZID, organism = "human", 
                       pAdjustMethod = "BH", readable = TRUE)
naive_RA<- enrichPathway(naive_mRNA_df$ENTREZID, organism = "human", readable = TRUE)

emapplot(naive_RA, showCategory = 10, layout = "kk")
cnetplot(naive_RA, showCategory = 3, colorEdge=TRUE)
dotplot(naive_RA)
dotplot(IRF3_RA)
emapplot(IRF3_RA, pie_scale= 1.5, layout = "kk")
cnetplot(IRF3_RA, showCategory = 5, circular=TRUE, colorEdge=TRUE) 

#visualizations of molecularfunction/biological process
naive_BP_simple<-clusterProfiler::simplify(naive_mRNA_BP)

p1 <- emapplot(IRF3_BP_simple, pie_scale = 1.5, layout = "kk")
p2 <- emapplot(IRF3_RA, pie_scale = 1.5, layout = "kk")
p3 <- emapplot(naive_mRNA_BP, pie_scale = 1.5, layout = "kk")
p4 <- emapplot(naive_RA, pie_scale = 1.5, layout = "kk")

cowplot::plot_grid(p1, p4, ncol=2, labels= LETTERS[1:4])


p1 <- barplot(naive_RA, showCategory = 10)
p2 <- dotplot(IRF3_RA, showCategory = 10)
p3 <- dotplot(ego_RA, showCategory = 10)

dotplot(naive_RA, IRF3_RA, ego_RA, showCategory = 20)

library(Cairo)
# Force R-studio to write output file type using Cairo, typically results in better formatting of image
require(Cairo)
# Initiate the file write session
CairoSVG(file = "KR_go.pdf", width = 7, height = 7)

cowplot::plot_grid(p1, p2, p3, ncols= 3, label_size = .1)

cairo_pdf(28, 8, file="KR_go.pdf", bg="transparent")
KR_go
dev.off()


#################################################################
#Go ontology for RIP and PARCLIP overlap 

#write_csv(IRFdata, "RIP_PAR_overlap_IRF3.csv")
#write_csv(naivedata, "RIP_PAR_overlap_naive.csv")
library(tidyverse)

IRF3_data<- read_csv("RIP_PAR_overlap_IRF3.csv")
naive_data<-read_csv("RIP_PAR_overlap_naive.csv")

potential_functional<- read_csv("functional_targets")
glimpse(potential_functional)

IRF3_RIP_df <- bitr(potential_functional$gene, fromType = "SYMBOL",
                    toType = c("ENTREZID", "ENSEMBL"),
                    OrgDb = org.Hs.eg.db)


IRF3_RIP_df_BP <- enrichKEGG(gene     = IRF3_RIP_df$ENTREZID,
                           rgDb         = org.Hs.eg.db,
                           ont           = "BP",
                           pAdjustMethod = "BH",
                           pvalueCutoff  = 0.1,
                           qvalueCutoff  = 0.5)
                           readable = TRUE, )

glimpse(IRF3_RIP_df_BP)
?enrichKEGG
IRF3_RIP_df_RA<-enrichPathway(IRF3_RIP_df$ENTREZID, organism = "human", 
                              pAdjustMethod = "BH", readable = TRUE,pvalueCutoff = .5, qvalueCutoff  = 0.5)

glimpse(IRF3_RIP_df_BP)

write.csv(IRF3_RIP_df_BP,"halflife_KEGG_functional_targets")

emapplot(IRF3_RIP_df_BP, showCategory = 10, layout = "kk")
cnetplot(IRF3_RIP_df_RA,colorEdge=TRUE, showCategory = 10, node_label="category", layout= "kk")
clusterProfiler::plot(IRF3_RIP_df_RA)

barplot(IRF3_RIP_df_BP) + 
  geom_bar(stat = "identity") +
  theme_minimal() +
  coord_flip() + 
  theme(axis.text.x = element_text(size=12), 
        axis.text.y = element_text(),
        panel.grid.major = element_blank(),
        strip.text.y = element_text(size=12)) +
  labs(y= "Number of transcripts", title = "GO analysis IRF3 specific")

IRF3_enriched<-data.frame(IRF3_RIP_df_RA)

### Trying to get enrichment values based on GO ontology. Binning by Go ontolgy 
library(tidyr)
library(tidyverse)
glimpse(IRF3_enriched)

IRF3_enriched <-IRF3_enriched %>% 
  mutate(SYMBOL =strsplit(geneID, "/")) %>% 
  unnest(SYMBOL) 

IRF3_enriched<-IRF3_enriched %>%
  dplyr::select(SYMBOL, ID, Description, p.adjust, qvalue, Count)

IRF3_data <- IRF3_data %>% 
  dplyr::rename(SYMBOL = gene_short_name) %>%
  dplyr::select(SYMBOL, normalized_Flag_naive,normalized_Flag_IRF3,naive_3UTR,IRF3_3UTR)

RIP_IRF3_enriched<-full_join(IRF3_data, IRF3_enriched, by = "SYMBOL")
RIP_IRF3_enriched<-na.omit(RIP_IRF3_enriched)

RIP_IRF3_enriched<-RIP_IRF3_enriched %>% 
  group_by(Description) %>%
  mutate(average = mean(normalized_Flag_naive))

RIP_IRF3_enriched_top<- RIP_IRF3_enriched %>%
  dplyr::select(Description, p.adjust, qvalue, Count, average) %>%
  unique() %>%
  arrange(p.adjust) %>%
  ungroup() %>%
  top_n(n = -20, wt= p.adjust)

ggplot(RIP_IRF3_enriched_top, aes(x= Description , y= average)) +
  geom_col() +
  coord_flip()
  theme (
    axis.text.x = element_text(size= 5, angle = 90))

glimpse(RIP_IRF3_enriched)

test2<-RIP_IRF3_enriched %>% 
  arrange(p.adjust) 
  ungroup() %>%
  top_n(n =  , wt= p.adjust)

length(unique(test2$Description))
  
ggplot(aes(x= normalized_Flag_IRF3, color = Description))+
  stat_ecdf()
  
###################################### 
#GO of high conficence 3'UTR targets with decreased half-life in KO  
halflife<-read_csv("functional_targets") %>% as.data.frame()
halflife<-halflife %>% filter(normalized_Flag_IRF3 > 0)
halflife_mRNA<- halflife$gene

halflife_df <- bitr(halflife_mRNA, fromType = "SYMBOL",
                     toType = c("ENTREZID", "ENSEMBL"),
                     OrgDb = org.Hs.eg.db)

halflife_df_BP <- enrichGO(gene     = halflife_df$ENTREZID,
                    OrgDb         = org.Hs.eg.db,
                    ont           = "BP",
                    pAdjustMethod = "BH",
                    readable = TRUE)

halflife_df_KEGG <- enrichKEGG(gene     = halflife_df$ENTREZID,
                           pvalueCutoff  = 0.1,
                           qvalueCutoff  = 0.5)

halflife_df_KEGG <-enrichPathway(halflife_df$ENTREZID, organism = "human",  readable = TRUE)

emapplot(halflife_df_BP, showCategory = 10, layout = "kk")
cnetplot(halflife_df_KEGG, showCategory = 10, colorEdge=TRUE, node_label="category" )
cnetplot(halflife_df_BP, showCategory = 10, colorEdge=TRUE, node_label="category")
barplot(halflife_df_RA)

barplot(halflife_df_KEGG) + 
  geom_bar(stat = "identity") +
  theme_minimal() +
  coord_flip() + 
  theme(axis.text.x = element_text(size=12), 
        axis.text.y = element_text(size=10),
        panel.grid.major = element_blank(),
        strip.text.y = element_text(size=12)) +
  labs(y= "Number of transcripts", title = "REACTOME ELAVL1 Targets affected half-life")


halflife_shiny_KEGG<-read_csv("KEGG_analaysis_3UTR_halflifetargets.csv")

halflife_shiny_KEGG<- halflife_shiny_KEGG %>%
  top_n(n = 20, wt= halflife_shiny_KEGG$`Enrichment FDR`)

glimpse(halflife_shiny_KEGG)

ggplot(halflife_shiny_KEGG, aes(y= `Genes in list`, x=`Functional Category`)) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  coord_flip() +
  theme(axis.text.x = element_text(size=12), 
        axis.text.y = element_text(size=10),
        panel.grid.major = element_blank(),
        strip.text.y = element_text(size=12)) +
  labs(y= "Number of transcripts", title = "REACTOME ELAVL1 Targets affected half-life")
  

halflife_shiny_KEGG<-halflife_shiny_KEGG

