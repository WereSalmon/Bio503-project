library(dplyr)
library(tidyr)
library(tidyverse)
####Data Load In####
lemay2018 <- readRDS("lemay2018_unfiltered_phyloseq.RDS")
#View(lemay2018@otu_table)
#View(lemay2018@sam_data)
#View(lemay2018@tax_table)
#Cleaning the sample data from lemay 2018
newMeta = separate(as.data.frame(as.matrix(lemay2018@sam_data)), 
                   col="sample_id",
                   into=c("sample_id", "Sample"),
                   sep="-")
newMeta$Region = str_sub(newMeta$sample_id, -2,-1)
newMeta$genus= str_sub(newMeta$sample_id, 1,-3)
newMeta$Species = newMeta[,"Species"] = NA
# Subsetting the king data to sample info that we need 
king2022 <- readRDS("king2022_unfiltered_phyloseq-2.RDS")
#View(king2022@otu_table)
#View(king2022@sam_data)
#View(king2022@tax_table)
kingsubsam <- subset(king2022@sam_data, select = 
                       c("sample_id","Sample", "Region","Species"))
kingsubsam = separate(as.data.frame(as.matrix(kingsubsam)),
                      col= "Species",
                      into=c("genus", "Species"),
                      sep=" ")
#combining sample data into 1 table
SamData <- rbind(newMeta,kingsubsam)
# new phyloseq object with all the data 
taxa1 <-  read.table("KingLemay_tax_SILVAv138_16s.txt", header = TRUE, row.names = 1)

OTUs1 <- t(readRDS("seqtab_nochim_KingLemay.RDS"))

combine_phylo = phyloseq(otu_table(OTUs1, taxa_are_rows = T),
                         tax_table(as.matrix(taxa1)),
                         sample_data((SamData)))

#### Cleaning ####

#### Cleaning ####
#Removing off-target taxa
combine_phylo <- subset_taxa(combine_phylo,
                       Kingdom != "Unassigned" &
                         Kingdom != "Eukaryota" &
                         Order != "Chloroplast"&
                         Order != "Mitochondria")
#Add unifiltered sample read numbers in 
combine_phylo@sam_data$sample_sums_unfiltered <- as.numeric(sample_sums(combine_phylo))
#sort(as.numeric(sample_sums(combine_phylo))) : Most samples > 1000,
combine_phylo_high <- prune_samples(sample_sums(combine_phylo) >= 1000, combine_phylo)
#extracting OTUs and filtering low frequency ASV's 
otutab <- as.data.frame((t(as.matrix(otu_table(combine_phylo_high@otu_table)))))
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

combine_phylo_clean = phyloseq(sample_data(combine_phylo_high),
                         tax_table(combine_phylo_high),
                         otu_table((otu.clean), taxa_are_rows = F))


combine_phylo_clean@sam_data$sample_sums_filtered <- sample_sums(combine_phylo_clean)

write_rds(combine_phylo_clean, "combine_phylo_clean.RDS")

#### Data now cleaned, ready for analysis #### 

