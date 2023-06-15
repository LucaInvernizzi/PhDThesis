# DATA for this version of the code is from:
# https://fungidb.org/fungidb/app/downloads/Current_Release/AfumigatusAf293/
# Locations are to be set accordingly by the user

### Create Genome Folder for STAR

STARbin <- "/localRaid/APPS/star/STAR"
SEQUENCEfile <- "/localRaid/DATA/database/Aspergillusgenome/sequenceFungiDB/FungiDB-52_AfumigatusAf293_Genome.fasta"
GFF3file <- "/localRaid/DATA/database/Aspergillusgenome/sequenceFungiDB/FungiDB-52_AfumigatusAf293.gff3"
OUTdir <- "/localRaid/DATA/database/Aspergillusgenome/STARindexFungiDB/"
ReadLength = 50
NumberOfThreads = 70


STARversion <- system(paste(STARbin, "--version"), intern = TRUE)
GTFversion = "VEuPathDB"

dir.create(OUTdir, recursive=T)
STARindex <- paste(OUTdir,"STARindex_star",STARversion,"_", GTFversion, "_ReadLength", ReadLength,sep="")

command <- paste(STARbin,
	"--runThreadN ", NumberOfThreads ,
	"--runMode ", "genomeGenerate", 
	"--genomeDir ", STARindex, 
	"--genomeFastaFiles ", SEQUENCEfile,
	"--sjdbGTFfile ", GFF3file,
	"--sjdbGTFtagExonParentTranscript Parent ",
	# in gtf there is a "transcript_id" col (the default of this command), the reference col here is "Parent"
	"--sjdbOverhang ", ReadLength-1,
	"--genomeSAindexNbases 11"
)
system(command)



### STAR allignment
setwd('/localRaid/DATA/CompBiolGroup/lucainvernizzi/project_inv_afum_mrna/data')

dir.create("STAROUT_notrim/", recursive=T)

files <- dir(path="originalfastq/", pattern=".fastq.gz$", full.names=T)
outnames <- dir(path="originalfastq/", pattern="fastq.gz$")
outnames <- gsub(".fastq.gz", "", outnames)


STARbin <- "/localRaid/APPS/star/STAR"
STARindex <- paste(OUTdir,"STARindex_star",STARversion,"_", GTFversion, "_ReadLength", ReadLength,sep="")

for (i in seq(from = 1, to = length(outnames))) {
  FILE1 <- files[grepl(outnames[i],files)][1] 
  # this is FILE 1 cause this way you can do FILE2 if it's pair ended
  
     outname <- outnames[i]
     
  if(file.exists(paste("STAROUT_notrim/" , outname, "/", sep=""))){
    print(paste(outname, "already done"))   
  }else{  dir.create(paste("STAROUT_notrim/" , outname, "/", sep=""))
               system(wait=T,paste(STARbin,  
                                   " --runThreadN 70  --genomeDir ", 
                                   STARindex,  " --outFilterType BySJout --outSAMtype BAM SortedByCoordinate  --quantMode TranscriptomeSAM GeneCounts --outSAMstrandField intronMotif    --outFilterIntronMotifs RemoveNoncanonical   --alignEndsType Local  --outSAMattributes XS  --readFilesCommand zcat --outFileNamePrefix STAROUT_notrim/" , 
                                   outname, "/ --readFilesIn " ,  
                                   FILE1, " ", 
                                   sep="")) # NB: Sorted by Coordinate is needed for IGV
                #system(wait=T,paste(STARbin,  " --runThreadN 70  --genomeDir ", STARindex,  " --outFilterType BySJout --outSAMtype BAM SortedByCoordinate  --quantMode TranscriptomeSAM GeneCounts --outSAMstrandField intronMotif    --outFilterIntronMotifs RemoveNoncanonical   --alignEndsType Local  --outSAMattributes XS  --readFilesCommand zcat --outFileNamePrefix STAROUT_notrim/" , outname, "/ --readFilesIn " ,  FILE2,  sep=""))

  }
print(date())
}