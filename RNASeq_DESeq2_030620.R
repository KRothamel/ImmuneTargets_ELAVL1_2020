#R code: Katie Rothamel, Kelly Barnett, and Tim Scott 
### Using Rsubread to calculate counts and DESEQ2 to calculate differentially expressed (DE) genes 
### Thank you Kelly Barnett and Tim Scott for script of on DE and lfcshrink using DeSeq2

library(tidyverse)
library(Rsubread)
library(dplyr)
library(pheatmap)
library(RColorBrewer)
library(DESeq2)
library(Cairo)

###Import STAR aligned BAM files

counts_0H <- featureCounts("THP1_0H_1_PE.bam", isPairedEnd=TRUE, annot.ext = "genes.gtf", isGTFAnnotationFile = TRUE)
counts_0H_2 <- featureCounts("THP1_0H_2_PE.bam",isPairedEnd=TRUE, annot.ext = "genes.gtf", isGTFAnnotationFile = TRUE)
counts_0H_3 <- featureCounts("THP1_0H_3_PE.bam",isPairedEnd=TRUE, annot.ext = "genes.gtf", isGTFAnnotationFile = TRUE)

counts_16H <- featureCounts("THP1_16H_1_PE.bam", isPairedEnd=TRUE, annot.ext = "genes.gtf", isGTFAnnotationFile = TRUE)
counts_16H_2 <- featureCounts("THP1_16H_2_PE.bam", isPairedEnd=TRUE, annot.ext = "genes.gtf", isGTFAnnotationFile = TRUE)
counts_16H_3 <- featureCounts("THP1_16H_3_PE.bam",isPairedEnd=TRUE, annot.ext = "genes.gtf", isGTFAnnotationFile = TRUE)


# Setup of experimental design and files
data<-cbind(counts_0H$counts, counts_0H_2$counts,counts_0H_3$counts,
            counts_16H$counts, counts_16H_2$counts, counts_16H_3$counts)

cts <- as.matrix(data,row.names=counts_0H$annotation$GeneID) 

#saveRDS(cts, "cts_030620") ##saves "cts" in R environment 

#load in counts data saved on Rproject folder
cts<-readRDS("cts")

colnames(cts)<-c("THP1_0H","THP1_0H_2", "THP1_0H_3", 
                 "THP1_16H", "THP1_16H_2", "THP1_16H_3")

RNames=c("THP1_0H","THP1_0H_2", "THP1_0H_3", 
         "THP1_16H", "THP1_16H_2", "THP1_16H_3")

condition <- c( "0hr", "0hr", "0hr",
                "16hr", "16hr", "16hr")

coldata <- data.frame(row.names = RNames, condition)

# Load DESeq2 library
library(DESeq2)

# Run standard differential expression analysis to measure the effect of "treatment"
dds <- DESeqDataSetFromMatrix(countData = cts,
                              colData = coldata,
                              design = ~ condition)
dds<-DESeq(dds)

#Removing lowly expressed genes. Filter genes with low counts (require >1 cpm in 1/2 of the samples) & normalize
keep <- rowSums(counts(dds)>1) >= 5 # identify the rows to keep that have cpm() > 1
dds_1 <- dds[keep, ] # remove the rows that do not have sufficient cpm

#List results
res<-results(dds_1)
resultsNames(dds_1)

# Log transforms the DESeq2 count table so that it is more appropriate for use in clustering or PCA
rld <- rlog(dds_1, blind=F) 

head(assay(rld))
head(assay(dds_1)) # this is just to look at the log2 data vs. count data

#Build new dataframe of shrunken log fold changes
resLFC<-lfcShrink(dds_1, coef = "condition_16hr_vs_0hr")

resLFC_df<-data.frame(rownames(resLFC), resLFC)

#Arrange dataframe by descending value of the log2 foldchange
resLFC_df<-arrange(resLFC_df, desc(resLFC_df$log2FoldChange))
head(resLFC_df, n = 20)
length(which(resLFC_df$log2FoldChang < -2 & resLFC_df$padj <=.005))

# Build similar dataframe but take the absolute value of log fold change
lfcShrink_df_absvalue <- data.frame(rownames(resLFC), abs(resLFC$log2FoldChange), resLFC$padj)
colnames(lfcShrink_df_absvalue) <- c("GeneID", "log2FoldChange", "padj")

# Filter lfcShrink_df_absvalue for padj value to remove low confidence genes 
lfcShrink_df_absvalue <- filter(lfcShrink_df_absvalue, padj <= .005)

#Filter for LogFC >= 1 and >= 3
lfcShrink_df_absvalue_filt <- filter(lfcShrink_df_absvalue, log2FoldChange >= 1)
lfcShrink_df_absvalue_filt2 <- filter(lfcShrink_df_absvalue, log2FoldChange >= 2)

# Sort data by abs(logFC) maxima
lfcShrink_df_absvalue_filt2 <- arrange(lfcShrink_df_absvalue_filt2, desc(log2FoldChange))

head(lfcShrink_df_absvalue_filt2)

# Filters the log transformed data against the p-adj values in the resTC table
rld.sig <- rld[ which(res$padj < 5*10^-3), ]
rld.sigtest2 <- rld[ which(lfcShrink_df_absvalue$log2FoldChange >= 1), ] %>% glimpse()
rld.sigtest3 <- rld[ which(lfcShrink_df_absvalue$log2FoldChange >= 3), ]

# Generate matrix that we can feed to pheatmap for clustering
log2_count_matrix <- as.data.frame(assay(rld), row.names=rownames(assay(rld)))
log2_count_matrix <- log2_count_matrix %>% mutate(gene_id = rownames(log2_count_matrix))
write.csv(log2_count_matrix, "DeSeq2_log2Counts")

#mat2 <- assay(rld.sigtest2)
#mat2.1 <- mat2[ which(rowSums2(mat2) != 0), ]

#mat3 <- assay(rld.sigtest3)
#mat3.1 <- mat3[ which(rowSums2(mat3) != 0), ]

# Re-order rows of  matrix for input into pheatmap
list <- c( "0hr", "0hr", "0hr", "16hr", "16hr", "16hr")

#datares <- pheatmap(log2_count_matrix, color = colorRampPalette(rev(brewer.pal(n = 11, name = "RdYlBu")))(100), cellwidth = 30, cellheight = 30, cluster_cols=F, kmeans_k = 7, scale = "row", cutree_rows=7)
#matrixwithclusters <- cbind(log2_count_matrix, datares$kmeans$cluster)
#cat(matrixwithclusters, sep="\t", file="matrixwithclusters.kmeans.txt")

#datares_HC <- pheatmap(mat3, color = colorRampPalette(rev(brewer.pal(n = 11, name = "RdYlBu")))(100), cluster_cols=F, scale = "row", annotation_names_row=F)

# Get gene lists that belong to each cluster 

#matrixwithclusters <- cbind(mat, datares$kmeans$cluster)
#write.table(matrixwithclusters, sep="\t", file="matrixwithclusters.kmeans.txt")

# Use any DESeqDataSet that contains the gene IDs as rownames and transform to a matrix
mat_all<-as.data.frame(log2_count_matrix)

# Optionally you can write this resulting table as a text file
write.table(mat_all, sep="\t", file="CompleteGeneList_DESeq2.txt")

############################################### Modified from ahkats Github (https://gist.github.com/akhats/73fdbbf45ebef560fc8c380cc11e6c00) 
# Paper: Mol Cell Biol 2006 Nov;26(21):7913-28.
# PMID: 16940178 (http://www.ncbi.nlm.nih.gov/pubmed/16940178)
# R code: Ahmed Moustafa

#Ploting expression expression plot and Volcano plot below
mat_all<- mat_all %>% mutate_at(vars(starts_with("T")), funs(as.numeric))
str(mat_all)

Zero = mat_all[,2:4]
Sixteen = mat_all[,5:7]

Zero.mean = apply(Zero, 1, mean)
Sixteen.mean = apply(Sixteen, 1, mean)
# Just get the maximum of all the means
limit = max(Zero.mean, Sixteen.mean)
plot(Zero.mean, Sixteen.mean, xlab = "baseline", ylab = "IRF3_stimulation",
     main = "GSE5583 - Scatter", xlim = c(-2, limit), ylim = c(-2, limit), pch=20)
fold = Sixteen.mean - Zero.mean
hist(fold, col = "gray", xlim=c(-5,5))

# Compute statistical significance (using t-test)
pvalue = NULL # Empty list for the p-values
tstat = NULL # Empty list of the t test statistics

for(i in 1:nrow(mat_all)) { # For each gene : 
  a = Zero[i,] # Zero of gene number i
  b = Sixteen[i,] # 16hrs of gene number i
  
  # Compute t-test between the two conditions
  t = t.test(unlist(a), unlist(b))
  
  # Put the current p-value in the pvalues list
  pvalue[i] = t$p.value
  # Put the current t-statistic in the tstats list
  tstat[i] = t$statistic
}

# Histogram of p-values (-log10)
hist(-log10(pvalue), col = "gray")
plot(fold, -log10(pvalue), main = "GSE5583 - Volcano", pch=20)
fold_cutoff = 2
pvalue_cutoff = 0.05
abline(v = fold_cutoff, col = "red", lwd = 3)
abline(v = -fold_cutoff, col = "blue", lwd = 3)
abline(h = -log10(pvalue_cutoff), col = "green", lwd = 3)
identify(fold, -log10(pvalue), rownames(mat_all))

#Fold-change filter for "biological" significance
filter_by_fold = abs(fold) >= fold_cutoff
x<-mat_all[filter_by_fold, ]

# P-value filter for "statistical" significance
filter_by_pvalue = pvalue <= pvalue_cutoff
filt_mat<-mat_all[filter_by_pvalue, ]

# Combined filter (both biological and statistical)
filter_combined = filter_by_fold & filter_by_pvalue
filtered = mat_all[filter_combined,]

# Let's generate the volcano plot again, highlighting the significantly differential expressed genes
plot(fold, -log10(pvalue), main = "GSE5583 - Volcano #2", pch=20)
points (fold[filter_combined], -log10(pvalue[filter_combined]),
        pch = 16, col = "red")

# Highlighting up-regulated in red and down-regulated in blue
plot(fold, -log10(pvalue), main = "GSE5583 - Volcano #3", pch=20)
points (fold[filter_combined & fold < 0],
        -log10(pvalue[filter_combined & fold < 0]),
        pch = 20, col = "red")
points (fold[filter_combined & fold > 0],
        -log10(pvalue[filter_combined & fold > 0]),
        pch = 20, col = "blue")


##################################
