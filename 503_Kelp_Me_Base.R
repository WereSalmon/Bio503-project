#### BIOL503 Kelp Core Microbiome ####
#### Packages ####
library(tidyverse)
library(phyloseq)
library(vegan)
library(plyr)
library(qualpalr)
library(ggpubr)
library(dplyr)
library(indicspecies)
library(data.table)
library(readr)


####Data Load In####
lemay2018 <- readRDS("lemay2018_unfiltered_phyloseq.RDS")
View(lemay2018@otu_table)
View(lemay2018@sam_data)
View(lemay2018@tax_table)

newMeta = separate(as.data.frame(as.matrix(lemay2018@sam_data)), 
                              col="sample_id",
                              into=c("kelp_location", "number"),
                              sep="-")


newMeta$location = str_sub(newMeta$kelp_location, -2,-1)
newMeta$kelp= str_sub(newMeta$kelp_location, 1,-3)

new_phyloseq = phyloseq(sample_data(newMeta),
                        otu_table(lemay2018@otu_table, taxa_are_rows = F),
                        tax_table(lemay2018@tax_table))



king2022 <- readRDS("king2022_unfiltered_phyloseq-2.RDS")
View(king2022@otu_table)
View(king2022@sam_data)
View(king2022@tax_table)

#### Binding dataframes together ####
kelp_me = merge_phyloseq(new_phyloseq, king2022)
write_rds(kelp_me, "kelp_me_unfiltered.RDS")

#### Cleaning ####
#Removing off-target taxa
kelp_me <- subset_taxa(kelp_me,
                       domain != "Unassigned" &
                         domain != "Eukaryota" &
                         order != "Chloroplast"&
                         order != "Mitochondria")
#Add unifiltered sample read numbers in 
kelp_me@sam_data$sample_sums_unfiltered <- as.numeric(sample_sums(kelp_me))
#sort(as.numeric(sample_sums(kelp_me))) : Most samples > 1000,
kelp_me_high <- prune_samples(sample_sums(kelp_me) >= 1000, kelp_me)
#extracting OTUs and filtering low frequency ASV's 
otutab <- as.data.frame((t(as.matrix(otu_table(kelp_me_high@otu_table)))))
otutab$asv_abundance <- rowSums(otutab)
#removing ASV's that occur fewer than 100 times in the dataset
otu.pruned <- subset(otutab, otutab$asv_abundance>100)
widthotu = ncol(otu.pruned)
otu.pruned <- otu.pruned[,-c(widthotu)]
#removing low freq samples 
ASVoccur = function(x){return(sum(x>0))}
otu.pruned$asv_occur_count = apply(otu.pruned,1, ASVoccur)
otu.highfreq = subset(otu.pruned, otu.pruned$asv_occur_count > 2)
otu.highfreq = otu.highfreq[,-c(widthotu)]

#### De Noising the data #### 
#converting reads that are likely barcode switching to 0
otu.clean <- mutate_all(otu.highfreq, funs(ifelse(. <5, 0, .)))

kelp_me_clean = phyloseq(sample_data(kelp_me_high),
                          tax_table(kelp_me_high),
                          otu_table(as.matrix(otu.clean), taxa_are_rows = T))

kelp_me_clean@sam_data$sample_sums_filtered <- sample_sums(kelp_me_clean)

