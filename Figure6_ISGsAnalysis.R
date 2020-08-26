library(tidyverse)

#####Performing CDF plot to see if ELAVL1 targets are ISGs

alldata<- read_csv("Halflife_alldata")
nontarget_ISGs <- read_tsv("nontargets_ISGs.txt", col_names = FALSE)
ThreeUTR_ISGs <- read_tsv("ISG_3UTR_targets.txt", col_names = FALSE)
nonThree_IsGs <- read_tsv("non3UTR_ISG_targets.txt", col_names = FALSE)

ISGs<-rbind(nontarget_ISGs, ThreeUTR_ISGs, nonThree_IsGs) %>% as.list()
glimpse(alldata)

alldata <- alldata %>% mutate(ISG= SYMBOL %in% ISGs$X1)

alldata2<-alldata  %>% 
  mutate(Category = case_when(
    IRF3_3UTR == 0 & IRF3_Intron == 0 & IRF3_Exon ==0 & IRF3_5UTR== 0 ~ "non-targets",
    IRF3_3UTR > 0 & normalized_Flag_IRF3 > 0 ~ "> 0 3'UTR sites_enriched"
  )) %>% na.omit()

alldata3<- alldata2 %>%  unite("z", ISG:Category, remove = FALSE)

#CDF plot looking at if the half-lives of ISG mRNA are more affected by the KO
alldata3 %>% filter(log2(halflife.y/halflife.x) != 0) %>% 
ggplot(aes(x = log2(halflife.y/halflife.x), color= z)) +
  stat_ecdf() +
  stat_ecdf( aes(y= log2(halflife.y/halflife.x), color = "black")) +
  theme_minimal() +
  xlim(-1.5,1.5)

#As a boxplot
alldata3 %>% 
  ggplot() +
  geom_boxplot(aes(y = log2(halflife.y/halflife.x), color= z, width= .5)) +
  #geom_boxplot( aes(y= log2(halflife.y/halflife.x), color = "black"))+
  theme_minimal() +
  ylim(-1,1)

#Filtering data
targets_ISGs <- alldata3 %>% filter(z == "TRUE_ > 0 3'UTR sites_enriched")
targets_notISGS<- alldata3 %>% filter(z == "FALSE_ > 0 3'UTR sites_enriched")
nontargets_ISGs<- alldata3 %>% filter(z == "TRUE_non-targets")

library(RColorBrewer)
#Looking at enrichment of ISGs vs their expression level (mRNA level)
ggplot(alldata3, aes(x= rank(normalized_Flag_IRF3)- rank(normalized_Flag_naive),
                    y = IRF3_mean_counts-naive_mean_counts, color = z )) +
  geom_point() +
  scale_color_manual(values=c("#a6611a", "#dfc27d", "#80cdc1", "#018571")) +
  theme_minimal()

### testing whether the change in half-life from KO/wt is significant from non-targets
x<- wilcox.test((nontargets_ISGs$halflife.y/nontargets_ISGs$halflife.x), (nontargets$halflife.y/nontargets$halflife.x))
p.adjust(x$p.value, method = "bonferroni", n = 650) 
x

alldata2  %>% filter(log2(halflife.y/halflife.x) != 0) %>% 
  mutate(Category = case_when(
    IRF3_3UTR == 0 & IRF3_Intron == 0 & IRF3_Exon ==0 & IRF3_5UTR== 0 ~ "non-targets",
    IRF3_3UTR > 0  ~ " > 0 3'UTR sites"
  )) %>%
  ggplot(aes(x = log2(halflife.y/halflife.x),  color= Category)) +
  stat_ecdf()+
  theme_minimal()

ISG_bound<-alldata3 %>% dplyr::filter(z == "TRUE_> 0 3'UTR sites_enriched" )
ISG_nonbound<-alldata3 %>% filter(z == "TRUE_non-targets") 
nonISG_bound <- alldata3 %>% filter(z == "FALSE_> 0 3'UTR sites_enriched")

ISG_bound_topenrich<- ISG_bound %>% filter(halflife.x - halflife.y > 1.5)

write_csv(ISG_bound_topenrich, "ISG_bound_topenrich")

#############################




