# Packages needed
library(edgeR)
library(DESeq2)
library(pcaExplorer)

### SET WD
setwd('/localRaid/DATA/CompBiolGroup/lucainvernizzi/project_inv_afum_mrna/data')
basedir <- getwd()

############# DIFFERENTIAL EXPRESSED GENES - DEGS ##########
# NB: this uses the non trimmed version of the data

counts.star.allrna <- read.table("counts_star_allrna_fungidb") # this was generated at line 78 of part 3_edgeR_DESeq_graphs.R

sampledescription <- read.table("GEOcode_name.csv", header = T, sep = ",")
row.names(sampledescription) <- paste(sampledescription$description, sampledescription$rep, sep = "_")

colnames(counts.star.allrna)
colnames(counts.star.allrna) <- paste(sampledescription$description, sampledescription$rep, sep = "_") # a bit crude


sum(rownames(sampledescription) != colnames(counts.star.allrna))


### DEGS analysis using DESeq2
asp_fum_star_notrim <- DESeqDataSetFromMatrix(countData = counts.star.allrna,
                                              colData = sampledescription,
                                              design = ~ description)
# converts some variable from character to factors

asp_fum_star_notrim <- DESeq2::DESeq(asp_fum_star_notrim)
vsdasp_fum_star_notrim <- DESeq2::vst(asp_fum_star_notrim)


# Exploration of the data with pcaExplorer
pcaExplorer(asp_fum_star_notrim)

### comparison
tmpresults_asp_fum_star <- results(asp_fum_star_notrim, contrast = c("description", "fe_starved_no_ox", "control"))

# some command to extract informations, can be adapted for other parameters
data_starved_vs_control <- subset(tmpresults_asp_fum_star, log2FoldChange > 1 )
data_starved_vs_control_sorted <- data_starved_vs_control[order(data_starved_vs_control$log2FoldChange), ]
head(data_starved_vs_control_sorted, 50)
tail(data_starved_vs_control_sorted, 50)

# creates a plot with the counts for the DEGS, to check a specific gene
plotCounts(gene = "Afu5g09730", asp_fum_star_notrim, intgroup = "description")

