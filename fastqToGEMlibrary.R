##########################################################################################
# This script creates a GEMcode library file from a FASTQ file (R1) for 1 sample
#
# Input: FASTQ file with R1 reads
##########################################################################################
args = commandArgs(trailingOnly=TRUE)

library(ShortRead)

setwd('/fastq')

fastqfile <- args[1]
outFile   <- args[2]   

x <- readFastq(fastqfile, full=TRUE, withIds=TRUE)  

N <- length( sread(x ))

A <- substr(x= as.character(sread(x) ) , start=1, stop=16  ) 

B <- lapply( as.character(id(x)) , strsplit, split=" ", fixed=TRUE )

B <- unlist(B)
SEL <- seq(from=1, to=length(B)-1, by=2)
B <- B[SEL]

# The output returns 2 columns: GEMcode and ReadID
write.table(x=cbind(A,B), file = outFile, append = TRUE, quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
# Example:
#GGTATTGAGGTAGCTG	NB501241:284:HJNWLBGXB:1:11101:21567:1115
#CGATGGCGTTATCACG	NB501241:284:HJNWLBGXB:1:11101:18967:1127
#CGATCAACAGACGTAG	NB501241:284:HJNWLBGXB:1:11101:14287:1128
#AGCAGTTAGCGTTGTT	NB501241:284:HJNWLBGXB:1:11101:18807:1150
#....

