library(dplyr)
library(xlsx)

b_load_clean <- function(file, species) {
        df_asp <- read.csv(paste0("data/", file), header = F)
        df_asp[,"query_subject_frame"] <- paste(df_asp[ , 14], df_asp[, 15], sep = "/")
        df_asp[, 14:15] <- NULL
        names(df_asp) <- c("name", "subject_acc_ver", "identity", "align_length", "mismatches", 
                           "gap_open", "query_start", "query_end", "sub_start", "sub_end", "e_value", "bit_score", 
                           "positives", "query_subject_frame")
        df_asp[, "species"] <- species
        return(df_asp)
}