# NB: this code needs the function defined in b_load_clean.R

# packages needed
library(dplyr)
library(tidyr)
library(xlsx)
library(ComplexHeatmap)

### load data into R 
clavatus_ptn_all_raw <- read.csv("data/clavatus_tblastn.csv", header = F)
flavus_ptn_all_raw <- read.csv("data/flavus_tblastn.csv", header = F)
fumigatus_ptn_all_raw <- read.csv("data/fumigatus_tblastn.csv", header = F)
nidulans_ptn_all_raw <- read.csv("data/nidulans_tblastn.csv", header = F)
niger_ptn_all_raw <- read.csv("data/niger_tblastn.csv", header = F)
oryzae_ptn_all_raw <- read.csv("data/oryzae_tblastn.csv", header = F)
terreus_ptn_all_raw <- read.csv("data/terreus_tblastn.csv", header = F)

### clean data fro each species
# CLAVATUS
df_asp <- clavatus_ptn_all_raw
df_asp$query_subject_frame <- paste(df_asp$V14, df_asp$V15, sep = "/")
df_asp[, 14:15] <- NULL

names(df_asp) <- c("name", "subject_acc_ver", "identity", "align_length", "mismatches", 
                   "gap_open", "query_start", "query_end", "sub_start", "sub_end", "e_value", "bit_score", 
                   "%positives", "query_subject_frame")
clavatus_ptn_all_clean <- df_asp

# FLAVUS
df_asp <- flavus_ptn_all_raw
df_asp$query_subject_frame <- paste(df_asp$V14, df_asp$V15, sep = "/")
df_asp[, 14:15] <- NULL
names(df_asp) <- c("name", "subject_acc_ver", "identity", "align_length", "mismatches", 
                   "gap_open", "query_start", "query_end", "sub_start", "sub_end", "e_value", "bit_score", 
                   "%positives", "query_subject_frame")
flavus_ptn_all_clean <- df_asp

# FUMIGATUS
df_asp <- fumigatus_ptn_all_raw
df_asp$query_subject_frame <- paste(df_asp$V14, df_asp$V15, sep = "/")
df_asp[, 14:15] <- NULL
names(df_asp) <- c("name", "subject_acc_ver", "identity", "align_length", "mismatches", 
                   "gap_open", "query_start", "query_end", "sub_start", "sub_end", "e_value", "bit_score", 
                   "%positives", "query_subject_frame")
fumigatus_ptn_all_clean <- df_asp

# NIDULANS
df_asp <- nidulans_ptn_all_raw
df_asp$query_subject_frame <- paste(df_asp$V14, df_asp$V15, sep = "/")
df_asp[, 14:15] <- NULL
names(df_asp) <- c("name", "subject_acc_ver", "identity", "align_length", "mismatches", 
                   "gap_open", "query_start", "query_end", "sub_start", "sub_end", "e_value", "bit_score", 
                   "%positives", "query_subject_frame")
nidulans_ptn_all_clean <- df_asp

# NIGER
df_asp <- niger_ptn_all_raw
df_asp$query_subject_frame <- paste(df_asp$V14, df_asp$V15, sep = "/")
df_asp[, 14:15] <- NULL
names(df_asp) <- c("name", "subject_acc_ver", "identity", "align_length", "mismatches", 
                   "gap_open", "query_start", "query_end", "sub_start", "sub_end", "e_value", "bit_score", 
                   "%positives", "query_subject_frame")
niger_ptn_all_clean <- df_asp

# ORYZAE
df_asp <- oryzae_ptn_all_raw
df_asp$query_subject_frame <- paste(df_asp$V14, df_asp$V15, sep = "/")
df_asp[, 14:15] <- NULL
names(df_asp) <- c("name", "subject_acc_ver", "identity", "align_length", "mismatches", 
                   "gap_open", "query_start", "query_end", "sub_start", "sub_end", "e_value", "bit_score", 
                   "%positives", "query_subject_frame")
oryzae_ptn_all_clean <- df_asp

# TERREUS
df_asp <- terreus_ptn_all_raw
df_asp$query_subject_frame <- paste(df_asp$V14, df_asp$V15, sep = "/")
df_asp[, 14:15] <- NULL
names(df_asp) <- c("name", "subject_acc_ver", "identity", "align_length", "mismatches", 
                   "gap_open", "query_start", "query_end", "sub_start", "sub_end", "e_value", "bit_score", 
                   "%positives", "query_subject_frame")
terreus_ptn_all_clean <- df_asp

# Add species col
clavatus_ptn_all_clean$species <- "clavatus"
flavus_ptn_all_clean$species <- "flavus"
fumigatus_ptn_all_clean$species <- "fumigatus"
nidulans_ptn_all_clean$species <- "nidulans"
niger_ptn_all_clean$species <- "niger"
oryzae_ptn_all_clean$species <- "oryzae"
terreus_ptn_all_clean$species <- "terreus"

# R BIND into a single object
ptn_complete <- rbind(clavatus_ptn_all_clean, flavus_ptn_all_clean, fumigatus_ptn_all_clean, nidulans_ptn_all_clean, 
                      niger_ptn_all_clean, oryzae_ptn_all_clean, terreus_ptn_all_clean)
# write object to a file for reference
write.csv(ptn_complete, "ptn_complete.csv")

### Protein name change
banana <- ptn_complete
banana$name <- gsub("_.*","", banana$name)
ptn_complete <- banana
nrow(ptn_complete)

### EXCEL WITH GENE/PROT INFO - Load and merge with dataframe
gp_info <- read.xlsx("table1.xlsx", sheetIndex = 1)
# nb the gene pigA is only pig cause it was wrong in fasta header

### MERGE INFORMATION FROM EXCEL AND DF FROM BLAST
info_blast_merge <- merge(ptn_complete, gp_info, by = "name")
nrow(info_blast_merge)
View(info_blast_merge)

### ADD QUERY COVER
# % of query lenght in alignment
names(info_blast_merge)

info_blast_merge$query_cover <- (info_blast_merge$query_end - info_blast_merge$query_start)/info_blast_merge$prot_length_uni*100
names(info_blast_merge)

summary(info_blast_merge$query_cover)
summary(info_blast_merge$identity)

### select coloumns for heatmap
hmap_df <- info_blast_merge %>%
        select("name", "species", "identity", "query_cover", "e_value") %>%
        group_by(name, species) %>%
        slice_max(identity, with_ties = F)
nrow(hmap_df)

thesis_heat <- hmap_df[, c("name", "species", "identity", "query_cover")]
thesis_heat <- unique(thesis_heat)
checkNA <- lapply(hmap_df, is.na)
lapply(checkNA, sum) # NA check

### cast the dataset into proper shape for heatmap, remove col with names, add rownames
cast_thesis_heat <- dcast(thesis_heat, species ~ name, value.var = "identity")
row.names(cast_thesis_heat) <- cast_thesis_heat$species
cast_thesis_heat <- subset(cast_thesis_heat, select = -c(species))

### heatmap generation (default is hierarchical clustering)
Heatmap(cast_thesis_heat)
