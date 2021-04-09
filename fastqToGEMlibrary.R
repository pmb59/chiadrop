#
# create GEMcode library from FASTQ file (R1) for 1 sample
#
# Input: FASTQ file with R1 reads
#

args = commandArgs(trailingOnly=TRUE)

library(ShortRead)

setwd('/fastq')

fastqfile <- args[1]
outFile   <- args[2]   

x <- readFastq(fastqfile, full=TRUE, withIds=TRUE)  
class(x)
head(sread(x ))
N <- length( sread(x ))
head(id(x ))
N


A <- substr(x= as.character(sread(x) )  , start=1, stop=16   ) 
#print(nchar(A) )
head(A)
B <- lapply( as.character(id(x)) , strsplit, split=" ", fixed=TRUE )
head(B)

B <- unlist(B)
SEL <- seq(from=1, to=length(B)-1, by=2)
B <- B[SEL]
length(B)
#B <- strsplit(x= as.character(id(x) )  , split=" ", fixed=TRUE )[[1]][1] 
#B


# The output returns 2 columns: GEMcode and ReadID
write.table(x=cbind(A,B), file = outFile, append = TRUE, quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
# Example:
#GGTATTGAGGTAGCTG	NB501241:284:HJNWLBGXB:1:11101:21567:1115
#CGATGGCGTTATCACG	NB501241:284:HJNWLBGXB:1:11101:18967:1127
#CGATCAACAGACGTAG	NB501241:284:HJNWLBGXB:1:11101:14287:1128
#AGCAGTTAGCGTTGTT	NB501241:284:HJNWLBGXB:1:11101:18807:1150
#....

#check at the end that IDs are unique
print( length(unique(  A  )) )
print( length(unique(  B  )) )

print( N )

rm(A)
rm(B)

