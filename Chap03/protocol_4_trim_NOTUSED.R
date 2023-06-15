### SET WD
setwd('/localRaid/DATA/CompBiolGroup/lucainvernizzi/project_inv_afum_mrna/data')
basedir <- getwd()

# as usual create some constant
files <- dir(path = "originalfastq/", pattern = ".fastq.gz$", full.names=T)
outnames <- dir(path="originalfastq/", pattern = "fastq.gz$")

# remove extension
outnames <- gsub(".fastq.gz", "", outnames)
# outnames <- unique(outnames)   (this makes outnames unique if double)

system("rm -r TRIMGALOREOUT/*") # remove folder if it exists already
dir.create("TRIMGALOREOUT/", recursive = T) # recursive here not really needed

# check if trimgalore is installed!

########## REMOVE 3 BASES AT 5' ########
for(i in seq(from = length(outnames),to = 1)){
        
        # again, if pair ended this needs to be modified
        FILE1 <- files[grepl(outnames[i], files)][1]
        
        outname <- outnames[i]
        if(file.exists(paste("TRIMGALOREOUT/", outname, "/", sep = ""))){
                print(paste(outname, "already done"))   
        }else{  dir.create(paste("TRIMGALOREOUT/" , outname, "/", sep=""))
                system(paste("trim_galore --gzip --cores 4  --illumina --clip_R1 3 --output_dir TRIMGALOREOUT/" , 
                             outname,  " ", FILE1, " ", 
                             sep=""), # --clip_R1 removes at 5', check trimgalore documentation for more specifications
                       wait=T)
                }
        print(date())
} # this takes some time

# 
system("rm -r trimmedfastq/*")
system("mkdir trimmedfastq") # sarebbe meglio: dir.create("trimmedfastq/", recursive = T)
system("find /localRaid/DATA/CompBiolGroup/lucainvernizzi/project_inv_afum_mrna/data/TRIMGALOREOUT/ -iname *.gz  -exec ln -s '{}' //localRaid/DATA/CompBiolGroup/lucainvernizzi/project_inv_afum_mrna/data/trimmedfastq/ ';'")

########## REMOVE  POLY-A ########
# everything could be done in one step, but this way you can check the output of the first one before the second step

setwd('/localRaid/DATA/CompBiolGroup/lucainvernizzi/project_inv_afum_mrna/data')
basedir <- getwd()
system("rm -r TRIMGALOREOUT_step2/*") # remove files if they exists


files <- dir(path = "trimmedfastq/", pattern = ".fq.gz$", full.names=T)
outnames <- dir(path="trimmedfastq/", pattern="fq.gz$")
outnames <- gsub(".fq.gz", "", outnames)
#outnames<-gsub("_L002_R1_001_trimmed.fq.gz","",outnames) # this is there are the paired end

dir.create("TRIMGALOREOUT_step2/",recursive=T)


# removes poly-A
for(i in seq(from = length(outnames), to=1)){
        
        #for(i in seq(from=1,to=2)){
        FILE1 <- files[grepl(outnames[i], files)][1]
        #FILE2 <-  files[grepl(outnames[i],files)][2]
        
        outname <- outnames[i]
        if(file.exists(paste("TRIMGALOREOUT_step2/" , outname, "/", sep = ""))){
                print(paste(outname, "already done"))   
        }else{  dir.create(paste("TRIMGALOREOUT_step2/" , outname, "/", sep=""))
                #system(paste("trim_galore --gzip --cores 2 -a 'A{10}' --hardtrim3 3 --output_dir TRIMGALOREOUT/" , outname, " ", FILE1, " ", sep=""),wait=F)
                system(paste("trim_galore --gzip --cores 4 -a 'A{10}'  --output_dir TRIMGALOREOUT_step2/" , outname, " ", FILE1, " ", sep=""), wait=T)
                
                
        }
        print(date())
        
}


system("rm -r trimmedfastq2/*")
system("mkdir trimmedfastq2")
system("find //localRaid/DATA/CompBiolGroup/lucainvernizzi/project_inv_afum_mrna/data/TRIMGALOREOUT_step2/ -iname *.gz  -exec ln -s '{}' /localRaid/DATA/CompBiolGroup/lucainvernizzi/project_inv_afum_mrna/data/trimmedfastq2/ ';'")


system("rm -r STAROUT_trim/*")
dir.create("STAROUT_trim/", recursive=T)

files <- dir(path = "trimmedfastq2/", pattern = ".fq.gz$", full.names=T)
outnames <- dir(path = "trimmedfastq2/", pattern = ".fq.gz$")
outnames <- gsub(".fq.gz", "", outnames)



STARbin="/localRaid/APPS/star/STAR"
STARindex <- paste(OUTdir,"STARindex_star",STARversion,"_", GTFversion, "_ReadLength", ReadLength,sep="")

for (i in seq(from=1,to=length(outnames))) {
        
        FILE1 <- files[grepl(outnames[i],files)][1]
        
        
        
        outname<-outnames[i]
        if(file.exists(paste("STAROUT_trim/" , outname, "/", sep=""))){
                print(paste(outname, "STAROUT_trim done"))   
        }else{  dir.create(paste("STAROUT_trim/" , outname, "/", sep=""))
                
                command <- paste(STARbin,  " --runThreadN 30  --genomeDir ", STARindex,  " --outFilterType BySJout --outSAMtype BAM SortedByCoordinate  --quantMode TranscriptomeSAM GeneCounts --outSAMstrandField intronMotif    --outFilterIntronMotifs RemoveNoncanonical   --alignEndsType Local  --outSAMattributes XS  --readFilesCommand zcat --outFileNamePrefix STAROUT_trim/" , outname, "/ --readFilesIn " ,  FILE1, " ", sep="")
                system(wait=T,command)
                #system(wait=T,paste(STARbin,  " --runThreadN 70  --genomeDir ", STARindex,  " --outFilterType BySJout --outSAMtype BAM SortedByCoordinate  --quantMode TranscriptomeSAM GeneCounts --outSAMstrandField intronMotif    --outFilterIntronMotifs RemoveNoncanonical   --alignEndsType Local  --outSAMattributes XS  --readFilesCommand zcat --outFileNamePrefix STAROUT_notrim/" , outname, "/ --readFilesIn " ,  FILE2,  sep=""))
                
        }
        print(date())
        
}

dir.create("QC")  # ?????? WHat does this does?
system('find STAROUT_trim -name "Aligned.sortedByCoord.out.bam" | parallel fastqc {} ',wait=F)