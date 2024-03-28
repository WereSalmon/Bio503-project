##BIOL503 Kelp Core Microbiome Project

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

##Load in cleaned combined phyloseq
kelp_me_clean <- readRDS("combine_phylo_clean.RDS")
View(kelp_me_clean@otu_table)
View(kelp_me_clean@sam_data)
View(kelp_me_clean@tax_table)

kelp_me_clean@sam_data$Group = ifelse(is.na(kelp_me_clean@sam_data$Species) != T, "UK", "BC")


####IndVal####

dephyloseq = function(phylo_obj){ 
  ## get the metadata
  meta = as.data.frame(as.matrix(phylo_obj@sam_data)) 
  ## how many metadata columns you have
  metacols = ncol(meta)+1
  ## get out the otu table
  #otu = as.data.frame(t(as.matrix(phylo_obj@otu_table)))
  otu = as.data.frame(as.matrix(phylo_obj@otu_table))
  ## merge the metadata and otu table by the rownames (sample ids from the Illumina sequencing data)
  mo = merge(meta, otu, by=0) 
  ## get out the taxonomy file
  tax = as.data.frame(phylo_obj@tax_table)
  ## get the ASV ID out. This the matches the placeholder ASV ID in the OTU table
  tax = tax %>% rownames_to_column(var="asv_id")
  tax = tax %>% rownames_to_column(var="asv_id2") # This is what I added to create the second column with asv_id as a just a number!
  ## pivot longer to be able to match the ASVs in the OTU table to the taxonomy table
  mo = mo %>% pivot_longer(cols = -c(1:metacols), names_to = "asv_id", values_to="asv_abundance") 
  ## Join the metadata and otu table with the taoxnomy table
  mot = full_join(mo, tax)
  ## Specify the output for the dephyloseq function
  output = mot 
}

#calculate sample sums of kelp me clean dataset
kelp_me_clean@sam_data$read_depth = sample_sums(kelp_me_clean)

##Formatting data for IndVal analysis
#sample IDs should be rows, metadata and ASVs should be columns

#remove samples with "NA" in genus in kelpme dataset - keep as phyloseq object
kelp_me_clean = subset_samples(kelp_me_clean, !(is.na(kelp_me_clean@sam_data$genus)))

#get the data you need
metadata = as.data.frame(kelp_me_clean@sam_data)
otu = as.data.frame(kelp_me_clean@otu_table)

#create a value with the comparison stored
comparison = metadata$genus

#Run indval to compare the host species
indval <- multipatt(otu, comparison, duleg = TRUE, control = how(nperm=999))

#Extracting the indval information - get output in a more convenient form

# Get indval statistic
indval.str <- as.data.frame(indval$str)
indval.str$rn <- rownames(indval.str)
# get p-value
indval.stat <- as.data.frame(indval$sign) #get dataframe of indval statistic
indval.stat$rn <- rownames(indval.stat) # make column of ASVs
# Prevalence as dataframe
indval.prev <- as.data.frame(indval$A)
# extract rownames into column
setDT(indval.prev, keep.rownames = TRUE)[]

## rename columns
colnames(indval.prev) <- paste0("prev.", colnames(indval.prev))
names(indval.prev)[names(indval.prev) == 'prev.rn'] <- 'rn'
# Specificity as dataframe
indval.sp <- as.data.frame(indval$B)
# extract rownames into column
setDT(indval.sp, keep.rownames = TRUE)[]

## rename columns
colnames(indval.sp) <- paste0("sp.", colnames(indval.sp)) 
names(indval.sp)[names(indval.sp) == 'sp.rn'] <- 'rn'
# Join statistics together
str.and.stat = full_join(indval.str, indval.stat, by="rn")
prev.and.fid = full_join(indval.prev, indval.sp, by="rn")
indval_table = full_join(str.and.stat, prev.and.fid, by="rn")

## get taxonomy table from phyloseq object
tax = as.data.frame(kelp_me_clean@tax_table)
## get the ASV ID into a new column
tax = tax %>% rownames_to_column(var="asv_id")
## rename columns fromt he IndVal output to join with taxonomy
names(indval_table)[names(indval_table) == 'rn'] <- 'asv_id'
## merge with taxonomy data
indval_table= inner_join(indval_table, tax)

view(indval_table)

##Cleaning up the IndVal data

#keep only ASVs that meet the filtering criteria >50% stat value
coreLaminaria = subset(indval_table,
                       #subset Laminaria
                       indval_table$index==5 &
                         #subset prevalence > 50%
                         indval_table$stat >= 0.5) 

coreSaccharina = subset(indval_table,
                        #subset Saccharina
                        indval_table$index==9 &
                          #subset prevalence > 50%
                          indval_table$stat >= 0.5) 

#save a .csv with the indval output
write.csv(coreLaminaria, "indval_output_core_laminaria.csv")
write.csv(coreSaccharina, "indval_output_core_saccharina.csv")

##Plotting IndVal data

##Plotting LAMINARIA core

#get a dataframe of the core microbiome datasets
## get a dataframe of the dataset
kelpme_df = dephyloseq(kelp_me_clean)

#use the Laminaria dataframe to only keep the taxa in the kelpmedf that are associated with Laminaria

laminariataxa = inner_join(kelpme_df, coreLaminaria)

## calculate the relative abundance of each ASV in each sample
laminariataxa$relative_abundance = as.numeric(laminariataxa$asv_abundance)/as.numeric(laminariataxa$read_depth)
## make a presence/abscence variable for each asv in each sample
laminariataxa$presabs <- if_else(laminariataxa$relative_abundance == 0, "asb.", "pres.")
## use this to make the 0s white later
laminariataxa$relative_abundance <- ifelse(laminariataxa$relative_abundance==0,NA,laminariataxa$relative_abundance)

#only plot one Group to keep the graph a reasonable size

UK_laminaria = subset(laminariataxa, laminariataxa$Group == "UK")

#save as .RDS
write_rds(UK_laminaria, "UK_laminaria_subset_indval.RDS")

#make bubble plot
ggplot(UK_laminaria, aes(x=as.character(Row.names),
                         y=paste0(genus, " (", asv_id2,")"),
                         size=relative_abundance))+
  geom_point()+
  theme(panel.grid = element_blank(),
        axis.text.y = element_text(size = 8, colour = "black", face="bold"),
        axis.text.x = element_blank(),
        axis.title = element_text(size=15),
        strip.text = element_text(color="black", size=12),
        legend.text=element_text(size=12),
        axis.line = element_line(colour = "black"),
        strip.text.y = element_text(angle = 0, size=8))+
  facet_grid(Family~Genus, space="free", scales="free")+
  labs(x="Samples", y = "Genus (asv_id)", size="relative
abundance")

##??????
#Convert filtered laminariataxa dataframe to phyloseq object

#remove NAs from genus
laminariataxa = na.omit(laminariataxa, Genus)



