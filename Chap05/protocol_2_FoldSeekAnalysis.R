library(dplyr)
library(tidyr)
library(ggplot2)
library("ggbreak")


### Load FoldSeek Data for Afu1g17270
# NB: unknow_value may be TM-score of alignment, irrelevant for analysis
fseek_table_afu1g17270_proteome <- read.table(
        "FoldSeek/Foldseek_afu1g17270/alis_afdb-proteome.m8",
        sep="\t", header = F, 
        col.names = c("query","target","fident","alnlen","mismatch","gapopen","qstart","qend","tstart","tend","evalue","bits", 
                      "querylength", "targetlength", "queryalign", "targetalign", "unkown_values", "fulltarget", "unkown_score", "speciesname"))

fseek_table_afu1g17270_swisspros <- read.table(
        "FoldSeek/Foldseek_afu1g17270/alis_afdb-swissprot.m8",
        sep="\t", header = F, 
        col.names = c("query","target","fident","alnlen","mismatch","gapopen","qstart","qend","tstart","tend","evalue","bits", 
                      "querylength", "targetlength", "queryalign", "targetalign", "unkown_values", "fulltarget", "unkown_score", "speciesname"))

fseek_table_afu1g17270_afdb50 <- read.table(
        "FoldSeek/Foldseek_afu1g17270/alis_afdb50.m8",
        sep="\t", header = F,
        col.names = c("query","target","fident","alnlen","mismatch","gapopen","qstart","qend","tstart","tend","evalue","bits", 
                      "querylength", "targetlength", "queryalign", "targetalign", "unkown_values", "fulltarget", "unkown_score", "speciesname"))

fseek_table_afu1g17270_pdb <- read.table(
        "FoldSeek/Foldseek_afu1g17270/alis_pdb100.m8",
        sep="\t", header = F, 
        col.names = c("query","target","fident","alnlen","mismatch","gapopen","qstart","qend","tstart","tend","evalue","bits", 
                      "querylength", "targetlength", "queryalign", "targetalign", "unkown_values", "fulltarget", "unkown_score", "speciesname"))

fseek_table_afu1g17270_gmcl_noLast2Cols <- read.table(
        "FoldSeek/Foldseek_afu1g17270/alis_gmgcl_id.m8",
        sep="\t", header = F, 
        col.names = c("query","target","fident","alnlen","mismatch","gapopen","qstart","qend","tstart","tend","evalue","bits", 
                      "querylength", "targetlength", "queryalign", "targetalign", "unkown_values", "fulltarget"))
fseek_table_afu1g17270_gmcl_noLast2Cols$unkown_score <- 0
fseek_table_afu1g17270_gmcl_noLast2Cols$speciesname <- "unknown"
head(fseek_table_afu1g17270_gmcl_noLast2Cols)

fseek_table_afu1g17270_mgnify_noLast2Cols <- read.table( # some col missing
        "FoldSeek/Foldseek_afu1g17270/alis_mgnify_esm30.m8",
        sep="\t", header = F, 
        col.names = c("query","target","fident","alnlen","mismatch","gapopen","qstart","qend","tstart","tend","evalue","bits", 
                      "querylength", "targetlength", "queryalign", "targetalign", "unkown_values", "fulltarget"))
fseek_table_afu1g17270_mgnify_noLast2Cols$unkown_score <- 0
fseek_table_afu1g17270_mgnify_noLast2Cols$speciesname <- "unknown"
head(fseek_table_afu1g17270_mgnify_noLast2Cols)


fseek_table_afu1g17270 <- rbind(fseek_table_afu1g17270_afdb50, fseek_table_afu1g17270_pdb,
                                fseek_table_afu1g17270_proteome, fseek_table_afu1g17270_swisspros,
                                fseek_table_afu1g17270_gmcl_noLast2Cols, fseek_table_afu1g17270_mgnify_noLast2Cols)

fseek_table_afu1g17270$query <- "afu1g17270"

# exploratory barplot for afu1g17270 data
ggplot(fseek_counts_table_afu1g17270, aes(x = qstart)) +
        geom_bar()

### Load FoldSeek Data for Afu8g01310
fseek_table_afu8g01310_proteome <- read.table(
        "FoldSeek/Foldseek_afu8g01310/alis_afdb-proteome.m8",
        sep="\t", header = F, 
        col.names = c("query","target","fident","alnlen","mismatch","gapopen","qstart","qend","tstart","tend","evalue","bits", 
                      "querylength", "targetlength", "queryalign", "targetalign", "unkown_values", "fulltarget", "unkown_score", "speciesname"))

fseek_table_afu8g01310_swisspros <- read.table(
        "FoldSeek/Foldseek_afu8g01310/alis_afdb-swissprot.m8",
        sep="\t", header = F, 
        col.names = c("query","target","fident","alnlen","mismatch","gapopen","qstart","qend","tstart","tend","evalue","bits", 
                      "querylength", "targetlength", "queryalign", "targetalign", "unkown_values", "fulltarget", "unkown_score", "speciesname"))

fseek_table_afu8g01310_afdb50 <- read.table(
        "FoldSeek/Foldseek_afu8g01310/alis_afdb50.m8",
        sep="\t", header = F,
        col.names = c("query","target","fident","alnlen","mismatch","gapopen","qstart","qend","tstart","tend","evalue","bits", 
                      "querylength", "targetlength", "queryalign", "targetalign", "unkown_values", "fulltarget", "unkown_score", "speciesname"))

fseek_table_afu8g01310_pdb <- read.table(
        "FoldSeek/Foldseek_afu8g01310/alis_pdb100.m8",
        sep="\t", header = F, 
        col.names = c("query","target","fident","alnlen","mismatch","gapopen","qstart","qend","tstart","tend","evalue","bits", 
                      "querylength", "targetlength", "queryalign", "targetalign", "unkown_values", "fulltarget", "unkown_score", "speciesname"))

fseek_table_afu8g01310_gmcl_noLast2Cols <- read.table(
        "FoldSeek/Foldseek_afu8g01310/alis_gmgcl_id.m8",
        sep="\t", header = F, 
        col.names = c("query","target","fident","alnlen","mismatch","gapopen","qstart","qend","tstart","tend","evalue","bits", 
                      "querylength", "targetlength", "queryalign", "targetalign", "unkown_values", "fulltarget"))
fseek_table_afu8g01310_gmcl_noLast2Cols$unkown_score <- 0
fseek_table_afu8g01310_gmcl_noLast2Cols$speciesname <- "unknown"
head(fseek_table_afu8g01310_gmcl_noLast2Cols)

fseek_table_afu8g01310_mgnify_noLast2Cols <- read.table( # some col missing
        "FoldSeek/Foldseek_afu8g01310/alis_mgnify_esm30.m8",
        sep="\t", header = F, 
        col.names = c("query","target","fident","alnlen","mismatch","gapopen","qstart","qend","tstart","tend","evalue","bits", 
                      "querylength", "targetlength", "queryalign", "targetalign", "unkown_values", "fulltarget"))
fseek_table_afu8g01310_mgnify_noLast2Cols$unkown_score <- 0
fseek_table_afu8g01310_mgnify_noLast2Cols$speciesname <- "unknown"
head(fseek_table_afu8g01310_gmcl_noLast2Cols)

fseek_table_afu8g01310 <- rbind(fseek_table_afu8g01310_afdb50, fseek_table_afu8g01310_pdb,
                                fseek_table_afu8g01310_proteome, fseek_table_afu8g01310_swisspros, 
                                fseek_table_afu8g01310_gmcl_noLast2Cols, fseek_table_afu8g01310_mgnify_noLast2Cols)
fseek_table_afu8g01310$query <- "afu8g01310"

## row bind 2 datasets
fseek_table_complete <- rbind(fseek_table_afu1g17270, fseek_table_afu8g01310)
fseek_table_complete <- fseek_table_complete[!duplicated(fseek_table_complete[,"target"]),]
summary(fseek_table_complete)

# Quantification of % of alignments starting in different positions
fsqstart <- fseek_table_complete$qstart
length(fsqstart[fsqstart < 21])
length(fsqstart[fsqstart < 21]) / length(fsqstart) * 100
length(fsqstart[fsqstart > 20 & fsqstart < 131]) / length(fsqstart) * 100
length(fsqstart[fsqstart > 130 & fsqstart < 431]) / length(fsqstart) * 100

length(fsqstart[fsqstart < 21]) / length(fsqstart) / 20
length(fsqstart[fsqstart > 20 & fsqstart < 131]) / length(fsqstart) * 100 /110
length(fsqstart[fsqstart > 130 & fsqstart < 431]) / length(fsqstart) * 100 / (430-130)


### GRAPHS FOR FOLDSEEK
test <- fseek_table_complete
test <-test[, c("qstart", "qend")]
colnames(test) <- c("Hit start", "Hit end")
test2 <- test %>%
        gather(1:2, key = "Alignment", value = "val")
head(test2)

test2$breaks <- cut(test2$val, breaks = c(0, 20, 130, 430, 600, 744))
test2$Domains <- cut(test2$val, breaks = c(0, 20, 130, 430, 600, 744))
levels(test2$Domains) <- c("Signal Peptide",
                           "DUF",
                           "TM Domain",
                           "FAD-binding",
                           "NAD-binding")

### Histogram generation
ggplot(test2, aes(x = val, colour = Alignment, fill = Domains)) +
        geom_histogram(binwidth = 20, alpha = 0.5, size = 1.2) +
        scale_colour_manual(values = c("sienna4", "grey10")) +
        labs(title = "Foldseek indication of domains preservation in the Fre protein family",
             x = "Aminoacidic Position",
             y = "Number of Hits") +
        scale_fill_manual(values = c("lightblue4", "springgreen4", "royalblue4", "gold", "khaki4")) +
        theme_bw(18) +
        scale_x_continuous(breaks = scales::pretty_breaks(n = 20)) +
        theme(plot.title = element_text(face = "bold", size = 15, hjust = 0.4))

### Alternative histogram with breaks (unused)
ggplot(test2, aes(x = val, colour = Alignment, fill = Domains)) +
        geom_histogram(binwidth = 10, alpha = 0.5, size = 1) +
        scale_colour_manual(values = c("brown4", "grey10")) +
        scale_y_break(c(520, 1000)) +
        labs(title = "Foldseek indication of domains preservation in the Fre protein family",
             x = "Aminoacidic Position",
             y = "Number of Hits") +
        scale_fill_manual(values = c("springgreen4", "lightblue2", "royalblue4", "gold", "seagreen4")) +
        theme_bw(18) +
        scale_x_continuous(breaks = scales::pretty_breaks(n = 20)) +
        theme(plot.title = element_text(face = "bold", size = 15, hjust = 0.4))