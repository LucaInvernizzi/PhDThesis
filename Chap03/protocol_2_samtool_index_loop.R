### Index the bam files

setwd('/localRaid/DATA/CompBiolGroup/lucainvernizzi/project_inv_afum_mrna/data/STAROUT_notrim/')

folder_vec <- dir(pattern = "SRR")
command <- paste("samtools", "index")

for (i in seq(from = 1, to = length(folder_vec))) {
        
        file_bam <- paste0(folder_vec[i], "/Aligned.sortedByCoord.out.bam")
        command_bam_to_index <- paste("samtools", "index", file_bam)
        system(command_bam_to_index)
        
}

