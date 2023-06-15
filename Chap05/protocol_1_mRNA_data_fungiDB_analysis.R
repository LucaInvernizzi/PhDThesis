library(ggplot2)
library(ggfortify)
library(ComplexHeatmap)
library(dplyr)
library(hrbrthemes)

### HEATMAP
# load csv file and calculated log2 values
getwd()
mRNA_fungiDB_data <- read.csv("Heatmap_mRNAdata/mRNA_data.csv",
                              sep = ";",
                              row.names = 1)

mRNA_fungiDB_data_log <- log(mRNA_fungiDB_data + 1)

# DEGs calculation on natural condition (no mutants or antifungal treatments)
mRNA_fungiDB_DEG_sel <- data.frame(row.names = row.names(mRNA_fungiDB_data_log))
mRNA_fungiDB_DEG_sel$hillmann_30_0 <- mRNA_fungiDB_data_log$hillmann_30_5 - mRNA_fungiDB_data_log$hillmann_0_100_CTR
mRNA_fungiDB_DEG_sel$irmer_b180_m180 <- mRNA_fungiDB_data_log$irmer_blood_b180 - mRNA_fungiDB_data_log$irmer_blood_m180
mRNA_fungiDB_DEG_sel$hagiwara_hyp_con <- mRNA_fungiDB_data_log$hagiwara_hyphae - mRNA_fungiDB_data_log$hagiwara_conidia
mRNA_fungiDB_DEG_sel$hagiwara_germ_con <- mRNA_fungiDB_data_log$hagiwara_germ - mRNA_fungiDB_data_log$hagiwara_conidia
mRNA_fungiDB_DEG_sel$kowalski_noox_CTR <- mRNA_fungiDB_data_log$kowalski_hyp_WT - mRNA_fungiDB_data_log$kowalski_nor_WT_CTR
mRNA_fungiDB_DEG_sel$kuruc_festarv_CTR <- mRNA_fungiDB_data_log$kurucz_festarv - mRNA_fungiDB_data_log$kurucz_CTR
mRNA_fungiDB_DEG_sel$kuruc_oxstress_CTR <- mRNA_fungiDB_data_log$kurucz_oxstress - mRNA_fungiDB_data_log$kurucz_CTR
mRNA_fungiDB_DEG_sel$losada_36_CTR <- mRNA_fungiDB_data_log$losada_t36 - mRNA_fungiDB_data_log$losada_t0_CTR

colnames(mRNA_fungiDB_DEG_sel) <- c("Oxigen limitation 1", #hillman
                                    "Grown in medium vs blood", # Irmer
                                    "Conidia vs Hyphae", # Hagiwara
                                    "Germinated conidia vs conidia", # Hagiwara
                                    "Oxigen limitation 2", # kowalski
                                    "Iron starvation", # Kuruc
                                    "Oxidative stress", # Kuruc
                                    "Oxigen limitation 3")

### Heatmap generation from DEGs values
heatmap_DEG_sel_fungidb <- Heatmap(as.matrix(mRNA_fungiDB_DEG_sel))

png("Heatmap_mRNAdata/complexheatmap_DEG_SEL_fungiDB.jpg", width = 960)
heatmap_DEG_sel_fungidb
dev.off()

### PCA generation from DEGs values
pca_DEG_SEL_fungiDB <- prcomp(mRNA_fungiDB_DEG_sel, 
                              center = TRUE,
                              scale. = TRUE)

summary(pca_DEG_SEL_fungiDB)
autoplot(pca_DEG_SEL_fungiDB, label = T)


### localisation data from deeploc 2.0
localization_signal <- read.csv("Localization/Summary_localization_15Genes_short.csv") # different file name
localization_signal$Protein_ID <- gsub("*_1", "", localization_signal$Protein_ID)
row.names(localization_signal) <- tolower(localization_signal$Protein_ID)
localization_signal <- subset(localization_signal, select = -Protein_ID)
View(localization_signal)

### Complete dataframe generation
# N term presence
n_term_factor <- read.csv("Heatmap_mRNAdata/mRNA_nterm_factor.csv",
                                  sep = ";",
                                  row.names = 1)

## presence of sreA/hapX target sequence
# starting from 
unique(names(gr_FRE_genes_sreA_mis))
unique(names(gr_FRE_genes_hapX_mis))
# this has been copied into a csv file manually
sreA_hapX_genes <- read.csv("Heatmap_mRNAdata/mRNA_sreA_hapX_target_sequence.csv",
                            sep = ";",
                            row.names = 1)


## Complete dataframe generation
merged_mRNA_factors <- cbind(mRNA_fungiDB_DEG_sel, localization_signal[,1:2])
merged_mRNA_factors <- cbind(merged_mRNA_factors, n_term_factor)
merged_mRNA_factors <- cbind(merged_mRNA_factors, sreA_hapX_genes)

merged_mRNA_factors$Localizations <- as.factor(merged_mRNA_factors$Localizations)

write.csv(merged_mRNA_factors, "Complete_DF_15proteins_info.csv")

### Lollipop plot
# Df preparation
Lol_2Genes_log <- mRNA_fungiDB_DEG_sel[c("afu1g17270", "afu8g01310"),]

Lol_2Genes_log2 <- data.frame(t(Lol_2Genes_log))
colnames(Lol_2Genes_log2) <- row.names(Lol_2Genes_log)
Lol_2Genes_log2$Experiments <- rownames(Lol_2Genes_log2)
head(Lol_2Genes_log2, Lol_2Genes_log2$Diff)

Lol_2Genes_log2$Diff <- abs(Lol_2Genes_log2$afu1g17270 - Lol_2Genes_log2$afu8g01310)
Lol_2Genes_log3 <- arrange(Lol_2Genes_log2, -Lol_2Genes_log2$Diff)

Lol_2Genes_log3 %>%
mutate(Diff=factor(Diff, levels=Diff))

# Lollipop plot generation
ggplot(Lol_2Genes_log3) +
        geom_segment( aes(x=reorder(Experiments, Diff), xend=Experiments, y=afu1g17270, yend=afu8g01310)) +
        geom_point( aes(x=Experiments, y=afu1g17270), color= "deepskyblue4", size=4 ) +
        geom_point( aes(x=Experiments, y=afu8g01310), color= "orange", size=4 ) +
        coord_flip()+
        theme_bw() +
        xlab("") +
        ylab("Difference in Log2FC")+
        labs(title = "FREB vs Afu8g01310 log2FC in Different Trascriptomic Experiments") +
        theme(plot.title = element_text(face = "bold", size = 10, hjust = 0.6), legend.position = "bottom")

# Lollipop plot saved into a file
ggsave("Lollipop.png", device = "png", dpi = 300, units = "cm", height = 11.39, width = 16.39)
