### bam files were aligned using STAR (see READme file in STAR_aligned folder)
###Using Rsubread to calculate counts and DESEQ2 to calculate differentially expressed (DE) genes 

library(tidyverse)
library(Rsubread)
library(dplyr)

###Import BAM files

counts_FLAG_0H_1 <- featureCounts("RIPseq_data/RIPseq_Sample3_PE.bam", isPairedEnd=TRUE,  annot.ext = "genes.gtf", isGTFAnnotationFile = TRUE)
counts_FLAG_0H_2 <- featureCounts("RIPseq_data/RIPseq_Sample6_PE.bam",isPairedEnd=TRUE, annot.ext = "genes.gtf", isGTFAnnotationFile = TRUE)

counts_FLAG_16H_1 <- featureCounts("RIPseq_data/RIPseq_Sample4_PE.bam",isPairedEnd=TRUE, annot.ext = "genes.gtf", isGTFAnnotationFile = TRUE)
counts_FLAG_16H_2 <- featureCounts("RIPseq_data/RIPseq_Sample5_PE.bam", isPairedEnd=TRUE, annot.ext = "genes.gtf", isGTFAnnotationFile = TRUE)

counts_IgG_0H_1 <- featureCounts("RIPseq_data/RIPseq_Sample7_PE.bam", isPairedEnd=TRUE, annot.ext = "genes.gtf", isGTFAnnotationFile = TRUE)
counts_IgG_0H_2 <- featureCounts("RIPseq_data/RIPseq_Sample9_PE.bam", isPairedEnd=TRUE, annot.ext = "genes.gtf", isGTFAnnotationFile = TRUE)

counts_IgG_16H_1 <- featureCounts("RIPseq_data/RIPseq_Sample11_PE.bam", isPairedEnd=TRUE,annot.ext = "genes.gtf", isGTFAnnotationFile = TRUE)
counts_IgG_16H_2 <- featureCounts("RIPseq_data/RIPseq_Sample12_PE.bam",isPairedEnd=TRUE, annot.ext = "genes.gtf", isGTFAnnotationFile = TRUE)


# Setup of experimental design and files
data<-cbind(counts_FLAG_0H_1$counts, counts_FLAG_0H_2$counts, 
            counts_FLAG_16H_1$counts,counts_FLAG_16H_2$counts,
            counts_IgG_0H_1$counts,counts_IgG_0H_2$counts,
            counts_IgG_16H_1$counts,counts_IgG_16H_2$counts)


rip_cts <- as.matrix(data,row.names=rownames(data))

colnames(rip_cts)<-c("FLAG_0H_1","FLAG_0H_2", 
                 "FLAG_16H_1","FLAG_16H_2", 
                 "IgG_0H_1","IgG_0H_2",
                 "IgG_16H_1", "IgG_16H_2")

RNames=c("FLAG_0H_1","FLAG_0H_2",  
           "FLAG_16H_1","FLAG_16H_2", 
           "IgG_0H_1","IgG_0H_2",
           "IgG_16H_1", "IgG_16H_2")
         
condition <- c("FLAG_0H", "FLAG_0H", 
                "FLAG_16H", "FLAG_16H", 
                "IgG_0H", "IgG_0H",  
                "IgG_16H", "IgG_16H")
coldata <- data.frame(row.names = RNames, condition)

head(rip_cts)


library(DESeq2)

# Run standard differential expression analysis
dds_rip <- DESeqDataSetFromMatrix(countData = rip_cts,
                              colData = coldata,
                              design = ~ condition)
dds_rip<-DESeq(dds_rip)
saveRDS(dds_rip)
#Removing lowly expressed genes. Filter genes with low counts (require >1 cpm in all of the samples) & normalize
keep <- rowSums(counts(dds_rip)>1) >= 3 # identify the rows to keep that have cpm() > 1
dds_rip_1 <- dds_rip[keep, ] # remove the rows that don't have sufficient cpm

#List results
res<-results(dds_rip_1)
resultsNames(dds_rip_1)

# Log transforms the DESeq2 count table so that it is more appropriate for use in clustering or PCA
rld_rip <- rlog(dds_rip_1, blind=F) 

RIP_data<-as.data.frame(assay(rld_rip))
RIP_data<- RIP_data %>% mutate(gene_id = rownames(RIP_data))

head(assay(dds_rip_1)) # this is just to look at the log2 data vs. count data

Zero_Flag = RIP_data[,1:2]
Sixteen_Flag = RIP_data[,3:4]
Zero_IgG = RIP_data[,5:6]
Sixteen_IgG = RIP_data[,7:8]

Flag.Zero.mean = apply(Zero_Flag, 1, mean)
Flag.Sixteen.mean = apply(Sixteen_Flag, 1, mean)
IgG.Zero.mean = apply(Zero_IgG, 1, mean)
IgG.Sixteen.mean = apply(Sixteen_IgG, 1, mean)

average_RIPdata<-cbind(Flag.Zero.mean, Flag.Sixteen.mean, IgG.Zero.mean, IgG.Sixteen.mean, RIP_data$gene_id) 
average_RIPdata<-as_data_frame(average_RIPdata)
write.csv(average_RIPdata, "averageRIPdata")

# Just get the maximum of all the means
limit = max(Zero.mean, Sixteen.mean)

plot(Flag.Sixteen.mean - IgG.Sixteen.mean, Flag.Zero.mean- IgG.Zero.mean)
abline(v = 1,col= "red")
abline(v = -1,col= "red")
abline(h = 1,col= "red")
abline(h = -1,col= "red")
identify(Flag.Sixteen.mean - IgG.Sixteen.mean, Flag.Zero.mean- IgG.Zero.mean, RIP_data$gene_id)


plot(average_RIPdata$Flag.Sixteen.mean, average_RIPdata$IgG.Sixteen.mean)
abline(a = 0, b= 1)

fold = Sixteen.mean - Zero.mean
hist(fold, col = "gray")

write.csv(average_RIPdata, "average_RIPdata")

