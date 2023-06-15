### SET WD
basewd <- c('/localRaid/DATA/CompBiolGroup/lucainvernizzi/project_inv_afum_mrna/data/')
setwd(basewd)
samples.star <- dir("STAROUT_notrim/") # this is a vector of all files/folder in this directory

# here some line to clean the names if needed (to be modified)
# samples.star.simple <- samples.star
# samples.star.simple<-gsub(pattern="_S.+","",samples.star.simple)

### MAKE A LIST OF TABLES FROM THE READSxGENE
counts.star.list<-list()

for (i in seq(from=1,to=length(samples.star))){
        
        tmpfile <- read.table(file=paste(basewd,"/STAROUT_notrim/",
                                         samples.star[i],
                                         "/ReadsPerGene.out.tab",
                                         sep=""),
                              skip=4,header=F,
                              sep="\t",stringsAsFactors=F)
        #tmpfile<-tmpfile[seq(from=1,to=(nrow(tmpfile)-4)),]
        counts.star.list[[i]] <- tmpfile
}

### HERE IT EXTRACT THE COL (2-, 3 and 4) AND MAKES A MATRIX
# col 1 is the geneIDs (not needed)

# COL 4
# extract
counts4.star.allrna.list <- lapply(counts.star.list, function(x)return(x[,4])) 
# cbind do.call is to passa a list to cbind
counts4.star.allrna <- do.call(cbind, counts4.star.allrna.list)
# use names of sample as colnames
colnames(counts4.star.allrna) <- samples.star
# use geneID as rownames
rownames(counts4.star.allrna) <- counts.star.list[[1]][,1]
# reorder (not needed but to keep everything more "clean")
counts4.star.allrna <- counts4.star.allrna[order(rownames(counts4.star.allrna)),]

# COL 3
counts3.star.allrna.list <- lapply(counts.star.list,function(x)return(x[,3]))
counts3.star.allrna <- do.call(cbind,counts3.star.allrna.list)
colnames(counts3.star.allrna) <- samples.star
rownames(counts3.star.allrna) <- counts.star.list[[1]][,1]
counts3.star.allrna <- counts3.star.allrna[order(rownames(counts3.star.allrna)),]

# COL 2
counts2.star.allrna.list <- lapply(counts.star.list,function(x)return(x[,2]))
counts2.star.allrna <- do.call(cbind,counts2.star.allrna.list)
colnames(counts2.star.allrna) <- samples.star
rownames(counts2.star.allrna) <- counts.star.list[[1]][,1]
counts2.star.allrna <- counts2.star.allrna[order(rownames(counts2.star.allrna)),]

####### Sum the col of the 3 matrixes, result is 3 vectors
# combine vectors in a matrix
# each col in the matrix is data from same col (fw, rv or unstranded)
# but from different experiments
counts.comparison <- cbind(
        count2_unstrand = apply(counts2.star.allrna, 2, sum),
        count3_FW = apply(counts3.star.allrna, 2, sum),
        count4_RV = apply(counts4.star.allrna, 2, sum)
)
counts.comparison

######### sum the different experiments to see where we have more data
counts.comparison_global <- apply(counts.comparison,2,sum)
counts.comparison_global
# find the highest one
max(counts.comparison_global)
# here is the unstranded

#### The best one is used as data
counts.star.allrna <- counts2.star.allrna
# saved as table with:
write.table(counts.star.allrna, file = "counts_star_allrna_fungidb.csv")

# to re read:
counts.star.allrna <- read.table("counts_star_allrna_fungidb.csv")
# but check wd