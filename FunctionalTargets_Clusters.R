library(tidyverse)
library(devtools)
library(dplyr)
library(GenomicRanges)
library(genomation)

###Import PAr-CLIP cluster data from IRF3 condition
IRF3_clusters <-read_csv("full_cluster_IRF3")
data<-read.csv("RIP_PAR_overlap_IRF3.csv")

### Import cluster data based on critieria. "I want all clusters that have CPM > 5" 
top_enrich <- data %>% filter(normalized_Flag_IRF3 > 2)
bottom_enrich<- data %>% filter(normalized_Flag_IRF3 < .5)
non_func<- read.csv("non_functional") #non functional targets 
func_targets<-read_tsv("test") #import list of function data

## From cluster data grabbing info needed to make a bed file (Chr, Start, End, Strand)
func_targets <- func_targets %>%
  rename(gene_short_name = SYMBOL) %>% 
  filter(IRF3_mean_counts/ naive_mean_counts > 1)

func_targets_cluster<- right_join(IRF3_clusters, func_targets, by = "gene_short_name" )
func_targets_cluster<-na.omit(func_targets_cluster)

non_func <- non_func %>%
  rename(gene_short_name = SYMBOL) 

non_func_targets <- right_join(IRF3_clusters, non_func, by = "gene_short_name")
non_func_targets_cluster<- na.omit(non_func_targets)

topenrich_cluster<-right_join(IRF3_clusters, top_enrich, by ="gene_short_name")
bottom_enrich <- right_join(IRF3_clusters, bottom_enrich, by ="gene_short_name")

topenrich_cluster<- na.omit(topenrich_cluster) %>% 
  filter(Aligned.to == "utr3") %>%
  dplyr::select(Chr, Start, End, Strand)

bottomenrich_cluster<- na.omit(bottom_enrich) %>% 
  filter(Aligned.to == "utr3") %>%
  dplyr::select(Chr, Start, End, Strand)

write_tsv(topenrich_cluster, "topenrich_cluster.bed")
write_tsv(bottomenrich_cluster, "bottomenrich_cluster.bed")

# Making a bed file of the clusters that are in the 3'UTR of our functional targets
func_target_bed<- func_targets_cluster %>%
  filter(Aligned.to == "utr3") %>%
  dplyr::select(Chr, Start, End, Strand)

non_func_targets_cluster<- non_func_targets_cluster %>%
  filter(Aligned.to == "utr3") %>%
  dplyr::select(Chr, Start, End, Strand)

write_tsv(non_func_targets_cluster, "non_func_targets_cluster.bed")

library(RCAS)
library(genomation)

#Import bed files using RCAS importBed 
queryRegions_func<-importBed("func_target_bed.bed")
queryRegions_topenrich<-importBed("topenrich_cluster.bed")
queryRegions_bottomenrich<- importBed("bottomenrich_cluster.bed")
queryRegions_allclusters<- importBed("IRF3_all_cluster.bed")
queryRegions_nonfunc<- importBed("non_func_targets_cluster.bed")

#Importing GTF that annotates the regions of interest (3'UTR) we want to overlap with our bed file 
gff<-import.gff("gencode.v19.chr_patch_hapl_scaff.annotation.gtf")
txdb <- getTxdbFeaturesFromGRanges(gff)
annotations<-unlist(txdb)
targetRegions_ThreeUTR <- txdb$threeUTRs
three_anno<-as.data.frame(targetRegions_ThreeUTR@elementMetadata@listData)
windows <- targetRegions_ThreeUTR[GenomicRanges::width(targetRegions_ThreeUTR) >= 10]

#ScoreMatrixBin
sm_func <-ScoreMatrixBin(target = queryRegions_func,
                            windows = windows,
                            bin.num = 100,
                            bin.op = "mean",
                            strand.aware = TRUE)


sm_topenrich <-ScoreMatrixBin(target = queryRegions_topenrich,
                         windows = windows,
                         bin.num = 100,
                         bin.op = "mean",
                         strand.aware = TRUE)

sm_bottomenrich <-ScoreMatrixBin(target = queryRegions_bottomenrich,
                         windows = windows,
                         bin.num = 100,
                         bin.op = "mean",
                         strand.aware = TRUE)

sm_all <-ScoreMatrixBin(target = queryRegions_allclusters,
                                 windows = windows,
                                 bin.num = 100,
                                 bin.op = "mean",
                                 strand.aware = TRUE)

sm_nonfunctional <-ScoreMatrixBin(target = queryRegions_nonfunc,
                        windows = windows,
                        bin.num = 100,
                        bin.op = "mean",
                        strand.aware = TRUE)

#S4 to dataframe
mdata_func <- data.frame(sm_func@.Data); dim(mdata_func)
mdata_topenrich<- data.frame(sm_topenrich@.Data)
mdata_bottomenrich <- data.frame(sm_bottomenrich@.Data)
mdata_all<- data.frame(sm_all@.Data)
mdata_non_fun<- data.frame(sm_nonfunctional)

#taking out rows with no information and annotating genes that overlapped with our windows 
mdata_func<- mdata_func %>% 
  filter(rowSums(mdata_func[,1:100]) !=0)
mfunc<-as.numeric(rownames(mdata_func))
annotated<-(as.data.frame(windows[mfunc,]))
anno_func_boundThreeUTRs<-cbind(annotated, mdata_func)
anno_func_boundThreeUTRs<- anno_func_boundThreeUTRs %>%
  mutate(source = "func_targets")

mdata_topenrich<- mdata_topenrich %>%
  filter(rowSums(mdata_topenrich[,1:100]) !=0)
mtopenrich<-as.numeric(rownames(mdata_topenrich))
annotated<-(as.data.frame(windows[mtopenrich,]))
anno_topenrich_boundThreeUTRs<-cbind(annotated, mdata_topenrich)
anno_topenrich_boundThreeUTRs <- anno_topenrich_boundThreeUTRs %>%
  mutate(source = "top_enriched")

mdata_bottomenrich<- mdata_bottomenrich %>% filter(rowSums(mdata_bottomenrich[,1:100]) !=0)
mbotenrich<-as.numeric(rownames(mdata_bottomenrich))
annotated<-(as.data.frame(windows[mbotenrich,]))
anno_bottomenrich_boundThreeUTRs<-cbind(annotated, mdata_bottomenrich)
anno_bottomenrich_boundThreeUTRs <- anno_bottomenrich_boundThreeUTRs %>%
  mutate(source = "bottom_enriched")

mdata_all<- mdata_all %>% filter(rowSums(mdata_all[,1:100]) !=0)
mall<-as.numeric(rownames(mdata_all))
annotated<-(as.data.frame(windows[mall,]))
anno_all_boundThreeUTRs<-cbind(annotated, mdata_all)
anno_all_boundThreeUTRs <- anno_all_boundThreeUTRs %>%
  mutate(source = "all")

mdata_non_fun<- mdata_non_fun %>% filter(rowSums(mdata_non_fun[,1:100]) !=0)
mnonfun<-as.numeric(rownames(mdata_non_fun))
annotated<-(as.data.frame(windows[mnonfun,]))
anno_nonfunctional_boundThreeUTRs<-cbind(annotated, mdata_non_fun)
anno_nonfunctional_boundThreeUTRs <- anno_nonfunctional_boundThreeUTRs %>%
  mutate(source = "non_functional")

#### wrangling data for ggplot 
all_anno_bins<- rbind(anno_bottomenrich_boundThreeUTRs, anno_topenrich_boundThreeUTRs, anno_func_boundThreeUTRs,
                      anno_all_boundThreeUTRs, anno_nonfunctional_boundThreeUTRs)

all_anno_bins <- all_anno_bins %>% 
  gather(bin, coverage, X1:X100) %>%
  group_by(bin, source) %>% 
  mutate(bin_mean= mean(coverage), 
         bin_sd = sd(coverage))

unique_all_anno_bins<- all_anno_bins[!duplicated(all_anno_bins[c('gene_name', 'bin', 'source')]),]

unique_to_plot<- unique_all_anno_bins %>% 
  as.data.frame() %>%
  mutate(bins = str_remove(bin, "X")) %>%
  mutate(bins= as.integer(bins)) 

#Plotting meta gene analysis with geom_smooth function. Similar to RCAS metagene analysis 
unique_to_plot %>% filter(source == c("func_targets" , "non_functional", "all")) %>%
ggplot(aes(x= bins, y = coverage, colour= source)) +
  geom_smooth(aes(fill= source, colour= source)) +
  theme_minimal() 

unique_to_plot %>% 
  filter(gene_name == "YTHDF1") %>%
ggplot( aes(x= bins, y = coverage)) +
  geom_line(color = "black") +
  geom_smooth(aes(x= bins, y = bin_mean), colour = "blue") +
  theme_minimal() +
  labs(x= "Distribution of Sites in the 3'UTR") +
  labs(y= "coverage") +
  labs(title = "Gene Specific 3'UTR binding")

test1<- anno_bottomenrich_boundThreeUTRs[ 8: 107] %>% transpose()
test2<- anno_func_boundThreeUTRs[8:107]

to_calculate<- group_by(unique_to_plot, source, bins) %>%
  summarise(
    count = n(),
    mean = mean(coverage, na.rm = TRUE),
    sd = sd(coverage, na.rm = TRUE))

to_heatmap<- to_calculate %>% 
  dplyr::select(source, mean, bins) %>%
  spread(source, mean)


#######################
#getting rid of duplicate gene values
unique_to_heatmap<- unique_to_plot %>% 
  filter(source == c("func_targets", " non_functional")) %>% 
  dplyr::select(source, coverage, bins) %>%
  spread(source, coverage)

str(to_heatmap)
datares_naive <- pheatmap(to_heatmap[,c(4:5)], cluster_rows = F,
                          color = colorRampPalette(rev(brewer.pal(n = 11, name = "RdYlBu")))(100), 
                          cellwidth = 50,cellheight = 3, cluster_cols= F)

datares_naive <- pheatmap(to_heatmap[,4:5], cluster_rows = F,
                          color = colorRampPalette(rev(brewer.pal(n = 11, name = "RdYlBu")))(1000), 
                          cellwidth = 50,cellheight = 3, cluster_cols= F)


func<- unique_to_plot  %>%
  dplyr::filter(source == "func_targets")

forheatmap<- func %>% 
  as.data.frame() %>%
  mutate(bins = str_remove(bin, "X")) %>%
  mutate(bins= as.integer(bins)) %>% 
  dplyr::select(bins, gene_name, coverage ) %>% 
  spread(bins, coverage)

library(pheatmap)
library(RColorBrewer)
library(Cairo)

query_regions_func <- pheatmap(forheatmap[,2:101], 
                          color = colorRampPalette(rev(brewer.pal(n = 11, name = "RdYlBu")))(1000), cellwidth = 3, cellheight = 30, cluster_cols= F, kmeans_k = 3, cutree_rows=3)

datares_naive <- pheatmap(forheatmap[,2:101], color = colorRampPalette(rev(brewer.pal(n = 11, name = "RdYlBu")))(1000), cellwidth = 5, cellheight = .1, cluster_cols= F)

summary(func_boundThreeUTRs)

#########################################################################
txdbFeatures <- getTxdbFeaturesFromGRanges(gff)
transcriptCoords <- GenomicFeatures::transcripts(txdb)

### coverage profile of query regions for all transcript features from RCAS. sampleN = 10000 for now to speed things up 
cvgList_func <- calculateCoverageProfileList(queryRegions = queryRegions_func,  targetRegionsList = unlist(windows), sampleN  = 10000)
write.csv(cvgList_naive, "cvgList_Naive")

### making seperate dataframe for particular transcript feature
cvg_3UTR_func<-cvgList_func %>% 
  filter(cvgList_func$feature == "threeUTRs") %>%
  mutate (source = "function targets")

#### binding data from naive and IRF3 so that I can used faucet in ggplot

cvg_all<-rbind(cvg_3UTR_naive, cvg_3UTR_IRF3)

ggplot2::ggplot(cvg_3UTR_func, aes(x = bins , y = meanCoverage)) + 
  geom_ribbon(fill = 'lightblue', 
              aes(ymin = meanCoverage - standardError * 1.96, 
                  ymax = meanCoverage + standardError * 1.96)) + 
  geom_line(color = 'black') + theme_bw(base_size = 14)


ggplot2::ggplot(to_plot, aes(x = as.integer(bin) , y = bin_mean))  + 
  geom_line(color = 'black') + theme_bw(base_size = 14)


ggplot(cvg_3UTR_naive, aes(x = bins , y = meanCoverage)) +
  geom_tile(aes(fill=meanCoverage, group=bins)) #, #color = "white") +
#scale_fill_gradient(low = "white", high = "steelblue")


head(cvg_3UTR_naive)


p2 <- ggplot( cvgT, aes(x = binds, y = ~meanCoverage), 
              type = 'scatter', mode = 'lines', 
              name = paste(featureName, "3' end coverage")) %>% 
  add_ribbons(x = ~bases, 
              ymin = cvgT$meanCoverage - cvgT$standardError*1.96, 
              ymax = cvgT$meanCoverage + cvgT$standardError*1.96, 
              line = list(color = 'rgba(7, 164, 181, 0.05)'),
              fillcolor = 'rgba(7, 164, 181, 0.2)',
              name = paste(featureName, "3' standard error (95% conf. int.)"))





#'
#' This function overlaps the input query regions with a target list of
#' annotation features and calculates the coverage profile along the target
#' regions.
#'
#' @param queryRegions GRanges object imported from a BED file using
#'   \code{importBed} function
#' @param txdb A txdb object obtained by using \code{GenomicFeatures::makeTxDb}
#'   family of functions
#' @param sampleN If set to a positive integer, the targetRegions will be
#'   downsampled to \code{sampleN} regions
#' @param type A character string defining the type of gene feature for which a
#'   profile should be calculated. The options are: transcripts, exons, introns,
#'   promoters, fiveUTRs, threeUTRs, and cds.
#' @return A data.frame object consisting of two columns: 1. coverage level 2.
#'   bins. The target regions are divided into 100 equal sized bins and coverage
#'   level is summarized in a strand-specific manner using the
#'   \code{genomation::ScoreMatrixBin} function.
#' @examples
#' data(gff)
#' data(queryRegions)
#' txdb <- GenomicFeatures::makeTxDbFromGRanges(gff)
#' df <- calculateCoverageProfileFromTxdb(queryRegions = queryRegions,
#'                                                type = 'exons',
#'                                                txdb = txdb,
#'                                             sampleN = 1000)
#' @import GenomicFeatures
#' @import GenomicRanges
#' @importFrom BiocGenerics unlist
#' @export
calculateCoverageProfileFromTxdb <- function (queryRegions,
                                              txdb,
                                              type,
                                              sampleN = 0) {
  
  if (type == 'transcripts') {
    targetRegions <- GenomicFeatures::transcripts(txdb)
  } else if (type == 'exons') {
    tmp <- GenomicFeatures::exonsBy(x = txdb, by = "tx", use.names = TRUE)
    targetRegions <- BiocGenerics::unlist(tmp)
  } else if (type == 'introns') {
    tmp <- GenomicFeatures::intronsByTranscript(x = txdb, use.names = TRUE)
    targetRegions <- BiocGenerics::unlist(tmp)
  } else if (type == 'promoters') {
    targetRegions <- GenomicFeatures::promoters(txdb)
  } else if (type == 'fiveUTRs') {
    tmp <- GenomicFeatures::fiveUTRsByTranscript(x = txdb, use.names = TRUE)
    targetRegions <- BiocGenerics::unlist(tmp)
  } else if (type == 'threeUTRs') {
    tmp <- GenomicFeatures::threeUTRsByTranscript(x = txdb, use.names = TRUE)
    targetRegions <- BiocGenerics::unlist()
  } else if (type == 'cds') {
    tmp <- GenomicFeatures::cdsBy(x = txdb, by = "tx", use.names = TRUE)
    targetRegions <- BiocGenerics::unlist(tmp)
  } else {
    stop ("Can calculate coverage profiles for only:
          transcripts, exons, introns,
          promoters, fiveUTRs, threeUTRs, and cds")
  }
  result = calculateCoverageProfile(queryRegions = queryRegions,
                                    targetRegions = targetRegions,
                                    sampleN = sampleN)
  return (result)
}

findLongLines <- function (myfile, lineLimit = 80) {
  counts <- lapply(X = readLines(myfile), FUN = function (x) {l=nchar(x)})
  names(counts) <- c(1:length(counts))
  return(names(counts)[counts > lineLimit])
}






