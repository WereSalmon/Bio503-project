##### set up #####
setwd("/parfreylab/BIOL403_503_datasets/king_lemay")

library(phyloseq)
library(tidyverse)
library(plyr)
library(dada2)


king = as.data.frame(t(readRDS("king2022_unfiltered_phyloseq-2.RDS")))
lemay = as.data.frame(t(readRDS("lemay2018_unfiltered_phyloseq.RDS")))

##### MATCH LENGTHS ####
king = king |> rownames_to_column(var="full_seq")
lemay = lemay |>  rownames_to_column(var="full_seq")

king$short_seq = str_sub(king$full_seq, start=1, end=212)
lemay$short_seq = str_sub(lemay$full_seq, start=1, end=212)

king = king[,-c(1)]
lemay = lemay[,-c(1)]

#https://stackoverflow.com/questions/10180132/consolidate-duplicate-rows
king.sum = ddply(king,"short_seq", numcolwise(sum))
lemay.sum = ddply(lemay,"short_seq", numcolwise(sum))


king.sum = king.sum |> 
  column_to_rownames(var="short_seq") |> as.matrix() |> t()
lemay.sum = lemay.sum |> 
  column_to_rownames(var="short_seq") |> as.matrix() |> t()


##### combine and run taxonomy #####

st.all <- mergeSequenceTables(king.sum, lemay.sum, tryRC=TRUE)
## remove bimeras
seqtab.nochim <- removeBimeraDenovo(st.all, method="consensus", multithread=TRUE,verbose = FALSE)
dim(seqtab.nochim)


## save files as RDS
write_rds(seqtab.nochim, "seqtab_nochim_KingLemay.RDS")


##### ASSIGN TAXONOMY ######

## SILVA needs to be downloaded in the directory to run
taxa <- assignTaxonomy(as.matrix(seqtab.nochim), "/parfreylab/schenk/taxonomy_databases/Silva138/silva_nr_v138_train_set.fa", 
                       multithread=TRUE, tryRC=TRUE)

# Save taxonomy table
write.table(data.frame("row_names"=rownames(taxa),taxa),"KingLemay_tax_SILVAv138_16s.txt", 
            row.names=FALSE, quote=F, sep="\t")

## add species information to taxonomy data
taxa_sp <- addSpecies(taxa, "/parfreylab/schenk/taxonomy_databases/Silva138/silva_species_assignment_v138.fa")

#NA taxa are hard to separate later if they have no label. apply "Unassigned" label here now.
unique(taxa_sp[,1]) #possible labels here: eukaryotic, archaeal, bacterial, and "NA" taxa. 
NAs <- is.na(taxa_sp[,1]) #test for NA
NAs <- which(NAs == TRUE) #get indices of NA values
taxa_sp[NAs,1] <- "Unassigned" #apply new label to identified indices
colnames(taxa_sp) <- c("domain", "phylum", "class", "order", "family", "genus", "species") 


# propagates taxonomy from left
tax_propagated = taxa_sp %>%
  t() %>% #transpose (moves taxonomy from column names to row names)
  na.locf() %>% #fill the NAs with the values from the cell to the left (higher taxonomic rank)
  t() |> # transpose back to have column names be taxonomy and row names be ASV
  as.data.frame() |> 
  rownames_to_column(var="asv_sequence") |> 
  rownames_to_column(var="asv_id") |> 
  column_to_rownames(var="asv_sequence")

# Save taxonomy table
write.table(data.frame("row_names"=rownames(tax_propagated),tax_propagated),
            "KingLemay_tax_noNAs_SILVAv138_16s.txt", 
            row.names=FALSE, quote=F, sep="\t")

#