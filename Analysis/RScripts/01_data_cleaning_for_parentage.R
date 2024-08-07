#####################
#     Libraries     #
#####################

library(adegenet)
library(poppr)
library(PopGenReport)
library(tidyverse)

###########################
#     Load Data Files     #
###########################

#set working directoy
setwd("/Users/mikaelyevans/Documents/GitHub/USBGHybridAcornsREU2023")

#load in genepop file as a genind object
UHA_genind <- read.genepop("Data_Files/Genotype_Files/2024_07_genepop.gen", ncode = 2)

#load score df
UHA_scores_df <- read.csv("Data_Files/CSV_Files/2024_07_UHA_database.csv")

###############################
#     Data Cleaning steps     #
###############################

#reduce genind file for individuals with greater than 25% missing data 
UHA_genind_nomd <- missingno(UHA_genind, type = "geno", 
                                cutoff = 0.25, quiet = FALSE, freq = FALSE)

#write out genind object as a genalex file
genind2genalex(UHA_genind_nomd,
               "Data_Files/CSV_Files/UHA_genalex_clean.csv", overwrite = TRUE)

#limit by the cleaned individuals
UHA_scores_clean_df <- UHA_scores_df[UHA_scores_df[,1] %in% 
                                       rownames(UHA_genind_nomd@tab),]

#write out 
write.csv(UHA_scores_clean_df, "Data_Files/CSV_Files/UHA_score_clean_df.csv")

###########################
#     Score Analysis      #
###########################

####Null alleles
#run null allele calculations over all genind objects
UHA_null_all <- null.all(UHA_genind_nomd)

#store in a data frame 
null_all_df <- signif(data.frame(UHA_null_all$null.allele.freq$summary2),3)*100

#write out null allele data frame
write.csv(null_all_df, "Results/Preliminary_Genotyping_Analysis/null_all_df.csv")

####Linkage Disequilibrium
#calculate linkage disequilbrium
UHA_ld <- pair.ia(UHA_genind_nomd, sample = 1000)

#convert to a data frame
UHA_ld_df <- data.frame(round(UHA_ld, digits = 2))

#write out ld df 
write.csv(UHA_ld_df, "Results/Preliminary_Genotyping_Analysis/ld_df.csv")

#########################################
#     Prep Data Files for Parentage     #
#########################################

##Create data frames with and without loci with 
#large number of null alleles

#first identify which alleles have > 10% null alleles
na_reorg_df <- as.data.frame(t(null_all_df))

#list of loci with > 15% null alleles 
na_loci <- rownames(na_reorg_df[which(na_reorg_df$`Observed frequency` > 15),])
#remove underscores 
na_loci <- gsub('_[^_]*$', '', na_loci)

#generate input data frame for parentage - all loci
write.csv(UHA_scores_clean_df, "Analysis/Parentage_Analysis/All_Loci/Input_Files/UHA_all_loci_genotype_df.csv",
          row.names = FALSE)

#generate a data frame removing loci with high null allele %
UHA_red_loci_df <- UHA_scores_clean_df %>%
                    dplyr::select(-contains(na_loci))

#generate parentage input data without loci with high null allele %
write.csv(UHA_red_loci_df, "Analysis/Parentage_Analysis/Red_Loci/Input_Files/UHA_red_loci_genotype_df.csv",
          row.names = FALSE)

