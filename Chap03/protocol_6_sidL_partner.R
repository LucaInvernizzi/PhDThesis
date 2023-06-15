### --- PLAN --- ###
# strategy is different with fungiDB
# 1 - Stats of everything with annotation (from step 5 and 7)
# 2 - DEG control vs -Fe-OX (lfc > 1)
# 3 - DEG control vs +Fe+OX (lfc > 2)
# 4 - extract granges for genes
# 5 - Search SreA sequence
# 6 - 
# 4 - Search annotation different from transacetylase

### --- 1 stats from all dataset --- ###

library(edgeR)
library(DESeq2)
library("biomaRt")

setwd('/localRaid/DATA/CompBiolGroup/lucainvernizzi/project_inv_afum_mrna/data')
basedir <- getwd()


# load the data and do some basic cleaning
counts.star.allrna <- read.table("counts_star_allrna_fungidb")
sampledescription <- read.table("GEOcode_name.csv", header = T, sep = ",")

row.names(sampledescription) <- paste(sampledescription$description, 
                                      sampledescription$rep, sep = "_")
colnames(counts.star.allrna)
colnames(counts.star.allrna) <- paste(sampledescription$description, sampledescription$rep, sep = "_") 
# ATTENTION THIS ASSUME THE COLS ARE IN THE RIGHT ORDER
# BETTER REWRITE IN A WAY THAT IT EXTRACT THE CORRECT ONE FROM SAMPLE DESCRIPTION


sum(sort(rownames(sampledescription)) != sort(colnames(counts.star.allrna)))

sampledescription$sampleID <- paste(sampledescription$description, 
                                    sampledescription$rep, sep = "_")

# actual analysis
asp_fum_star_notrim <- DESeqDataSetFromMatrix(countData = counts.star.allrna,
                                              colData = sampledescription,
                                              design = ~ description)

asp_fum_star_notrim <- DESeq2::DESeq(asp_fum_star_notrim)
vstasp_fum_star_notrim <- DESeq2::vst(asp_fum_star_notrim)

### 2 ask values Silvano
result_ctr_vs_nofe <- results(asp_fum_star_notrim, contrast = c("description", "fe_starved_no_ox", "control"))
topgenes_ctr_vs_nofe_index <- which((result_ctr_vs_nofe[,'padj']<0.05 & result_ctr_vs_nofe[ ,'log2FoldChange'] > 1))
print(paste("Ndiff=", length(topgenes_ctr_vs_nofe_index)))
topgenes_ctr_vs_nofe <- rownames(result_ctr_vs_nofe[topgenes_ctr_vs_nofe_index, ])

### 3
result_ctr_vs_ox <- results(asp_fum_star_notrim, contrast = c("description", "ox_stress_with_fe", "control"))
topgenes_ctr_vs_ox_index <- which((result_ctr_vs_ox[,'padj']<0.05 & result_ctr_vs_ox[ ,'log2FoldChange'] > 1.5))
print(paste("Ndiff=", length(topgenes_ctr_vs_ox_index)))
topgenes_ctr_vs_ox <- rownames(result_ctr_vs_ox[topgenes_ctr_vs_ox_index, ])

### ALTRO
result_ox_vs_nofeox <- results(asp_fum_star_notrim, contrast = c("description", "fe_starved_ox_stress", "ox_stress_with_fe"))
topgenes_ox_vs_nofeox_index <- which((result_ox_vs_nofeox[,'padj']<0.05 & result_ox_vs_nofeox[ ,'log2FoldChange'] > 1))
print(paste("Ndiff=", length(topgenes_ox_vs_nofeox_index)))
topgenes_ox_vs_nofeox <- rownames(result_ox_vs_nofeox[topgenes_ox_vs_nofeox_index, ])

result_nofe_vs_nofeox <- results(asp_fum_star_notrim, contrast = c("description", "fe_starved_ox_stress", "fe_starved_no_ox"))
topgenes_nofe_vs_nofeox_index <- which((result_nofe_vs_nofeox[,'padj']<0.05 & result_nofe_vs_nofeox[ ,'log2FoldChange'] > 1.5))
print(paste("Ndiff=", length(topgenes_nofe_vs_nofeox_index)))
topgenes_nofe_vs_nofeox <- rownames(result_nofe_vs_nofeox[topgenes_nofe_vs_nofeox_index, ])

# topgenes_ctr_vs_nofe - topgenes_ctr_vs_ox = 1018 with lfc 1 - 1.5
topgenes <- topgenes_ctr_vs_nofe[! topgenes_ctr_vs_nofe %in% topgenes_ctr_vs_ox]
length(topgenes)

### get ranges from gene name using the gff3 file ###
# read gff into R and make into df
gff <- rtracklayer::import("/localRaid/DATA/database/Aspergillusgenome/sequenceFungiDB/FungiDB-52_AfumigatusAf293.gff3")
gff_df <- as.data.frame(gff)
# extract only rows and cols needed
gff_df_genes <- gff_df[which(gff_df$type == "protein_coding_gene"),]
gff_df_genes <- gff_df_genes[, c("ID", "start", "end", "strand", "seqnames")]
head(gff_df_genes)

topgenes_gff_df <- gff_df_genes[gff_df_genes$ID %in% topgenes,]

library(GenomicRanges)
gr_topgenes <- GRanges(Rle(topgenes_gff_df$seqnames),
                       IRanges(topgenes_gff_df$start,
                               topgenes_gff_df$end,
                               names = topgenes_gff_df$ID),
                       Rle(strand(topgenes_gff_df$strand)))


table(seqnames(gr_topgenes))
table(strand(gr_topgenes))

### Ranges for target sreA
library(Biostrings)
library(BSgenome.Afumigatus.fungiDB.Af293)
Afumigatus

srea_A <- vmatchPattern("ATCAGATAA", Afumigatus) 
srea_T <- vmatchPattern("ATCTGATAA", Afumigatus)

srea_combined <- c(srea_A, srea_T)
table(seqnames(srea_combined))

# tst to see if there is some overlapping in the sreA sequences
srea_A_nostrand <- srea_A
strand(srea_A_nostrand) <- Rle(rep("*")) # access the strand and insert repetition of *
srea_T_nostrand <- srea_T
strand(srea_T_nostrand) <- Rle(rep("*"))

srea_overlaps <- findOverlaps(srea_A_nostrand, srea_T_nostrand)
# 14 are hit
(527-14)/16

countOverlaps(srea_A_test, srea_T_test)

srea_combined_nostrand <- c(srea_A_nostrand, srea_T_nostrand)
# if needed

### overlaps with and without strand
promoter_topgenes <- promoters(gr_topgenes) # promotors from topgenes

overlaps <- findOverlaps(srea_combined, promoter_topgenes)
# promoter WITH strand -> 40 hits
subjectHits(overlaps)
gr_best_genes <- gr_topgenes[subjectHits(overlaps),]
gr_best_genes
best_genes <- names(gr_best_genes)

overlaps_nostrand <- findOverlaps(srea_combined_nostrand, promoter_topgenes)
# promoter WITHOUT strand -> 74
subjectHits(overlaps_nostrand)
gr_best_genes_nostrand <- gr_topgenes[subjectHits(overlaps_nostrand),]
gr_best_genes_nostrand
best_genes_nostrand <- names(gr_best_genes)

### Add annotation, file location needs to be modified accordingly
gff <- rtracklayer::import("/localRaid/DATA/database/Aspergillusgenome/sequenceFungiDB/FungiDB-52_AfumigatusAf293.gff3")
gff_df <- as.data.frame(gff)
# generated at line 74 first time, just to be sure the object exists

gff_df_genes_ann <- gff_df[which(gff_df$type == "protein_coding_gene"),]
gff_df_genes_ann <- gff_df_genes_ann[, c("ID", "description")]
head(gff_df_genes_ann)

### Extract gene with specific annotation
acyl <- gff_df_genes_ann[grep(pattern = "acyl", gff_df_genes_ann$description),]
acetyl <- gff_df_genes_ann[grep(pattern = "acetyl", gff_df_genes_ann$description),]
hypo <- gff_df_genes_ann[grep(pattern = "hypo", gff_df_genes_ann$description),]
putative <- gff_df_genes_ann[grep(pattern = "putative", gff_df_genes_ann$description),]

topgenes_ann_regex <- unique(rbind(hypo, acyl, acetyl, putative))
# unique dataframe with all unique proteins

### extract only the hits from strand and nostrand
final_genes <- topgenes_ann_regex[topgenes_ann_regex$ID %in% best_genes,]
final_genes_nostrand <- topgenes_ann_regex[topgenes_ann_regex$ID %in% best_genes_nostrand,]
# this gives the same result... 7 genes...
# lot of them are already known/can't be correct

df_result_ctr_vs_nofe <- as.data.frame(result_ctr_vs_nofe)
df_result_ctr_vs_nofe$ID <- rownames(df_result_ctr_vs_nofe) # adding an ID col

final_genes_pvalue <- df_result_ctr_vs_nofe[df_result_ctr_vs_nofe$ID %in% final_genes$ID,]
final_gene_sidL_table <- cbind(final_genes[, c("ID", "description")], final_genes_pvalue)
write.table(final_gene_sidL_table, "final_gene_sidL_table.csv", sep = ";", quote = F)