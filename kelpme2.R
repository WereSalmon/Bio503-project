## load libraries
library(phyloseq) 
library(tidyverse)

## set your working directory
setwd("/Users/alice/Desktop/UBC/503/Lab")
# read in merged file
## read in non-rarefied data
kelp_me_clean = readRDS("combine_phylo_clean.RDS")

#view(kelp_me_clean@sam_data)


# remove NA from england region, replace with UK and replace NA from BC with BC
# Filter rows
kelp_me_clean_filtered2 <- subset_samples(kelp_me_clean, !(sample_names(kelp_me_clean) %in% c("ERR8442991", "ERR8442995", "ERR8442996", "ERR8443009","ERR8442960","ERR8442941")))
#view(kelp_me_clean_filtered2@sam_data)
kelp_me_clean_filtered2@sam_data[, "Region"] <- lapply(kelp_me_clean_filtered2@sam_data[, "Region"], function(x) {
  ifelse(x %in% c("W Scotland", "N Scotland", "S England", "Wales"), "UK", 
         ifelse(x %in% c("WB", "SI", "TB","GI"), "BC", x))
})
#view(kelp_me_clean_filtered2@sam_data)

#asv_id
#view(kelp_me_clean_filtered@tax_table)
#view(kelp_me_clean_filtered2@otu_table)


####indval####
## load libraries
library(phyloseq) 
library(tidyverse)
#install.packages("indicspecies")
library(indicspecies)
#install.packages("data.table")
library(data.table)
library(dplyr)

#### Dephyloseq function #### 
dephyloseq <- function(phylo_obj){
  meta = as.data.frame(as.matrix(phylo_obj@sam_data))
  metacols = ncol(meta)+1
  otu = as.data.frame(as.matrix(phylo_obj@otu_table))
  mo = merge(meta,otu, by=0)
  tax = as.data.frame(phylo_obj@tax_table)
  tax = tax %>% rownames_to_column(var = "asv_id2")
  mo = mo %>% pivot_longer(cols = -c(1:metacols), names_to = "asv_id2", 
                           values_to =  "asv_abundance")
  mot = full_join(mo,tax)
  output = mot
}


## calculate sample sums of the dataset
kelp_me_clean_filtered2@sam_data$read_depth = sample_sums(kelp_me_clean_filtered2)  
#view(kelp_me_clean_filtered2@sam_data)
## get data you need, converted with sample ID's as rows, remove NA's from region in phyloseq
kelp_me_clean_filtered2@sam_data$read_depth = sample_sums(kelp_me_clean_filtered2)
metadata = as.data.frame(as.matrix(kelp_me_clean_filtered2@sam_data))
otu = as.data.frame(kelp_me_clean_filtered2@otu_table)

## create a value with the comparison stored
comparison = metadata$Region


## run indval to compare the region
indval <- multipatt(otu, comparison, duleg = TRUE, control = how(nperm=999))

# Get indval statistic
indval.str <- as.data.frame(indval$str) 
indval.str$rn <- rownames(indval.str)
# get p-value
indval.stat <- as.data.frame(indval$sign) 
#get dataframe of indval statistic
indval.stat$rn <- rownames(indval.stat) # make column of ASVs
# Prevalence as dataframe
indval.prev <- as.data.frame(indval$A)
# extract rownames into column
setDT(indval.prev, keep.rownames = TRUE)[]


#rename columns 
colnames(indval.prev) <- paste0("spec.", colnames(indval.prev))
names(indval.prev)[names(indval.prev) == 'spec.rn'] <- 'rn' 

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


#get taxonomy table from phyloseq object 
#kelp_me_clean_tax <- as.data.frame(kelp_me_clean_filtered2@tax_table)
##get ASV ID into new column 
#kelp_me_clean_tax = kelp_me_clean_tax %>% rownames_to_column(var="asv_name")
## rename columns fromt he IndVal output to join with taxonomy
#names(indval_table)[names(indval_table) == 'rn'] <- 'asv_name'
## merge with taxonomy
#indval_table= inner_join(indval_table, kelp_me_clean_tax)

## get taxonomy table form phyloseq object
tax = as.data.frame(kelp_me_clean_filtered2@tax_table)
## get the ASV ID into a new column
tax = tax %>% rownames_to_column(var="asv_id2")
## rename columns fromt he IndVal output to join with taxonomy
names(indval_table)[names(indval_table) == 'rn'] <- 'asv_id2'
## merge with taxonomy
indval_table= inner_join(indval_table, tax)




## keep only ASVs that meet our filtering criteria
kelpmenow = subset(indval_table, 
                   indval_table$stat >= 0.8)
## save a .csv with the indval output
write.csv(kelpmenow, "indval_output_kelpme.csv")

## get a dataframe of the dataset
kelpmenowdf = dephyloseq(kelp_me_clean_filtered2)
## calcuate the relative abundance of each ASV in each sample
kelpmenowdf$relative_abundance = as.numeric(kelpmenowdf$asv_abundance)/as.numeric(kelpmenowdf$read_depth)

## use the dataframe to only keep the taxa in the df that meet filtering criteria
taxa = inner_join(kelpmenow, kelpmenowdf)


## get the ASV ID  in
taxa2 = taxa %>% rownames_to_column(var="asv_idd")


## calcuate the relative abundance of each ASV in each sample
taxa2$relative_abundance = as.numeric(taxa2$asv_abundance)/as.numeric(taxa2$read_depth)
## make a presence/abscence variable for each asv in each sample
# you don't need this today but it is helpful
taxa2$presabs <- if_else(taxa2$relative_abundance == 0, "asb.", "pres.")
## use this to make the 0s white later
taxa2$relative_abundance <- ifelse(taxa2$relative_abundance==0,NA,taxa2$relative_abundance)

## lets only plot one kelp to keep the graph a reasonable size
Saccharina = subset(taxa2,taxa2$genus =="Saccharina")
## lets only plot one kelp to keep the graph a reasonable size
Laminaria = subset(taxa2,taxa2$genus =="Laminaria")
#adapted from this to take out asv id - ggplot(Saccharina, aes(x=as.character(Row.names), y=paste0(species.x",(",asv_idd,")"),size=relative_abundance)) # y = "Genus (asv_id)"

#plot to  look at the class by order distribution of Saccharina
ggplot(Saccharina, aes(x=as.character(Row.names), y=paste0(Order),
                       size=relative_abundance))+
  geom_point()+
  theme(panel.grid = element_blank(),
        axis.text.y = element_text(size = 8, colour = "black", face="bold"), axis.text.x = element_blank(),
        axis.title = element_text(size=15),
        strip.text = element_text(color="black", size=12), legend.text=element_text(size=12),
        axis.line = element_line(colour = "black"),
        strip.text.y = element_text(angle = 0, size=8))+ facet_grid(Class~Region, space="free", scales="free")+ labs(x="Samples", y = "Order", size="relative_abundance")

#plot to just look at the class distribution of Saccharina
ggplot(Saccharina, aes(x=as.character(Row.names),y=paste0(Class),
                       size=relative_abundance))+
  geom_point()+
  theme(panel.grid = element_blank(),
        axis.text.y = element_text(size = 8, colour = "black", face="bold"), axis.text.x = element_blank(),
        axis.title = element_text(size=15),
        strip.text = element_text(color="black", size=12), legend.text=element_text(size=12),
        axis.line = element_line(colour = "black"),
        strip.text.y = element_text(angle = 0, size=8))+ facet_grid(Class~Region, space="free", scales="free")+ labs(x="Samples", y = "Class", size="relative_abundance")

#plot to  look at the class by order distribution of Laminaria
ggplot(Laminaria, aes(x=as.character(Row.names), y=paste0(Order),
                       size=relative_abundance))+
  geom_point()+
  theme(panel.grid = element_blank(),
        axis.text.y = element_text(size = 8, colour = "black", face="bold"), axis.text.x = element_blank(),
        axis.title = element_text(size=15),
        strip.text = element_text(color="black", size=12), legend.text=element_text(size=12),
        axis.line = element_line(colour = "black"),
        strip.text.y = element_text(angle = 0, size=8))+ facet_grid(Class~Region, space="free", scales="free")+ labs(x="Samples", y = "Order", size="relative_abundance")

#plot to just look at the class distribution of laminaria
ggplot(Laminaria, aes(x=as.character(Row.names),y=paste0(Class),
                       size=relative_abundance))+
  geom_point()+
  theme(panel.grid = element_blank(),
        axis.text.y = element_text(size = 8, colour = "black", face="bold"), axis.text.x = element_blank(),
        axis.title = element_text(size=15),
        strip.text = element_text(color="black", size=12), legend.text=element_text(size=12),
        axis.line = element_line(colour = "black"),
        strip.text.y = element_text(angle = 0, size=8))+ facet_grid(Class~Region, space="free", scales="free")+ labs(x="Samples", y = "Class", size="relative_abundance")



