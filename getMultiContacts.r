#----------------------------------------------------------------------------------------------------
# Chromatin Interaction Detection at Single Molecule Resolution
#----------------------------------------------------------------------------------------------------

args = commandArgs(trailingOnly=TRUE)

#-------------------------------------------------------------------
# Configured parameters: 
#-------------------------------------------------------------------
GapWidth <- as.numeric( args[1] )  #bp
NCORES=2
setwd('/chiadrop/sampleX')
outdir <- '/chiadrop/sampleX/outs/'
header <- '/chiadrop/sampleX/HEADER_miaSig.txt'
SampleID <- 'sampleX'
GEMraw <- 'sampleX.gemLibrary'  # file with GEM barcodes and readIDs 
BED    <- 'sampleX.bed'   # BED file with alignments (do not have GEM barcode)
BAM <- "/chiadrop/sampleX/q40_sampleX_sorted.bam"
genome <- 'Mus_musculus.GRCm38.dna.primary_assembly.fa.fai'  # genome chr sizes
QC_barcodes <- FALSE  # {FALSE,TRUE}   ( computationally very slow, might rescue only up to 13% ) 
ReadLengthMin <- 10    #35 #50 # bp
MinReadsInEachFragment <- 2
FragmentNumberCutoff <- 2    # min MUST be 2
#-------------------------------------------------------------------
#-------------------------------------------------------------------


library(parallel)  
options(cores = NCORES) 
getOption('cores')


start_time <- format(Sys.time(), "%a %b %d %X %Y")

write.table("--------------------------" , file = paste0(outdir,SampleID,'_GapWidth_',GapWidth,".log"), append = TRUE, quote = FALSE, sep = "\t", row.names = FALSE, col.names =FALSE)
write.table(paste0('Starting analysis of sample: ',SampleID) , file = paste0(outdir,SampleID,'_GapWidth_',GapWidth,".log"), append = TRUE, quote = FALSE, sep = "\t", row.names = FALSE, col.names =FALSE)
write.table(paste0('GapWidth: ',GapWidth) , file = paste0(outdir,SampleID,'_GapWidth_',GapWidth,".log"), append = TRUE, quote = FALSE, sep = "\t", row.names = FALSE, col.names =FALSE)
write.table(paste0('ReadLengthMin: ',ReadLengthMin) , file = paste0(outdir,SampleID,'_GapWidth_',GapWidth,".log"), append = TRUE, quote = FALSE, sep = "\t", row.names = FALSE, col.names =FALSE)
write.table(paste0('FragmentNumberCutoff: ',FragmentNumberCutoff) , file = paste0(outdir,SampleID,'_GapWidth_',GapWidth,".log"), append = TRUE, quote = FALSE, sep = "\t", row.names = FALSE, col.names =FALSE)
write.table(paste0('Rescue Barcodes: ',QC_barcodes) , file = paste0(outdir,SampleID,'_GapWidth_',GapWidth,".log"), append = TRUE, quote = FALSE, sep = "\t", row.names = FALSE, col.names =FALSE)
write.table("--" , file = paste0(outdir,SampleID,'_GapWidth_',GapWidth,".log"), append = TRUE, quote = FALSE, sep = "\t", row.names = FALSE, col.names =FALSE)
write.table(start_time , file = paste0(outdir,SampleID,'_GapWidth_',GapWidth,".log"), append = TRUE, quote = FALSE, sep = "\t", row.names = FALSE, col.names =FALSE)

#-------------------------------------------------------------------
# aux functions
#-------------------------------------------------------------------
library(seqinr)
library(GenomicRanges)

GEMnumber <- function(INPUTSEQ){
  return(  paste0( s2n(seq=strsplit(tolower(INPUTSEQ), split='')[[1]], levels = s2c("acgtn"), base4 = TRUE, forceToLower = FALSE) , collapse='') )
}

setMethod("$<-", "GRanges", function(x, name, value) { 
  elementMetadata(x)[[ name ]] <- value
  return(x)
}) 


#-------------------------------------------------------------------
# raw data stats, and reading GEM barcodes of raw reads 
#-------------------------------------------------------------------
library(data.table)
mydata <- fread(GEMraw, header=FALSE)
head(mydata)
dim(mydata)

TA <- table(mydata$V1)
head(TA)


write.table(paste0('Number of total raw reads in FASTQ is....', nrow(mydata)) , file = paste0(outdir,SampleID,'_GapWidth_',GapWidth,".log"), append = TRUE, quote = FALSE, sep = "\t", row.names = FALSE, col.names =FALSE)
write.table(paste0('Number of unique raw GEM barcodes is....', length(   unique(names(TA))   )) , file = paste0(outdir,SampleID,'_GapWidth_',GapWidth,".log"), append = TRUE, quote = FALSE, sep = "\t", row.names = FALSE, col.names =FALSE)
write.table("--" , file = paste0(outdir,SampleID,'_GapWidth_',GapWidth,".log"), append = TRUE, quote = FALSE, sep = "\t", row.names = FALSE, col.names =FALSE)

png(paste0(outdir ,GEMraw, '_GapWidth_',GapWidth , '.png'), height=400, width=705)
  par(mfrow=c(1,2))
  let <- length(sort(TA, decreasing = TRUE) )
  plot( as.numeric(sort(TA, decreasing = TRUE) )[1:let ] , lwd=2, type='h', log='x', col='magenta', ylab='Number of reads', xlab='',
      main='Frequency of GEM barcodes \n (Log x)', sub=paste0('\n Number of unique GEM barcodes: ', length(   unique(names(TA))   ), '\n Number of total reads in FastQ: ', nrow(mydata))) 
  plot( as.numeric(sort(TA, decreasing = TRUE) )[1:let ] , lwd=2, type='h', log='xy', col='magenta', ylab='Number of reads', xlab='',
      main='Frequency of GEM barcodes \n (Log xy)', sub=paste0('\n Number of unique GEM barcodes: ', length(   unique(names(TA))   ), '\n Number of total reads in FastQ: ', nrow(mydata))) 
  rm(let)
dev.off()

#Number of unique GEM barcodes
length(   unique(names(TA))   )

write.table(x=as.data.frame(sort(TA, decreasing = TRUE) ), file =  paste0(outdir ,GEMraw, '_GapWidth_',GapWidth, '.txt'), 
            append = FALSE, quote = FALSE, sep = "\t",
            row.names = FALSE, col.names = TRUE)

write.table("Top 10 GEMcodes in raw data:" , file = paste0(outdir,SampleID,'_GapWidth_',GapWidth,".log"), append = TRUE, quote = FALSE, sep = "\t", row.names = FALSE, col.names =FALSE)
write.table( as.data.frame(sort(TA, decreasing = TRUE) )[1:10,] , file = paste0(outdir,SampleID,'_GapWidth_',GapWidth,".log"), append = TRUE, quote = FALSE, sep = "\t", row.names = FALSE, col.names =FALSE)


rm(TA)



#read file with length of each chromosome
genomeSizes <- read.table(genome, head=FALSE)[,c(1,2)]
head(genomeSizes)


#read bed file
myBED <- fread(BED, header=FALSE)
myBED <- myBED 
myBED <- myBED[, c(1:4,6,11) ]
nrow(myBED)

write.table("--" , file = paste0(outdir,SampleID,'_GapWidth_',GapWidth,".log"), append = TRUE, quote = FALSE, sep = "\t", row.names = FALSE, col.names =FALSE)
write.table(paste0("Number of R1/R2 reads to process: ", nrow(myBED) ), file = paste0(outdir,SampleID,'_GapWidth_',GapWidth,".log"),append = TRUE, quote = FALSE, sep = "\t", row.names = FALSE, col.names =FALSE)
write.table(paste0('read length MIN = ',min(myBED$V11), '; read length MAX = ', max(myBED$V11)) , file = paste0(outdir,SampleID,'_GapWidth_',GapWidth,".log"),append = TRUE, quote = FALSE, sep = "\t", row.names = FALSE, col.names =FALSE)
write.table(paste0("Chromosome included: ", unique(myBED$V1) ), file = paste0(outdir,SampleID,'_GapWidth_',GapWidth,".log"),append = TRUE, quote = FALSE, sep = "\t", row.names = FALSE, col.names =FALSE)
write.table("--" , file = paste0(outdir,SampleID,'_GapWidth_',GapWidth,".log"), append = TRUE, quote = FALSE, sep = "\t", row.names = FALSE, col.names =FALSE)



#retain reads with length >= ReadLengthMin
myBED <- myBED[ which(myBED$V11 >= ReadLengthMin ) , ]
nrow(myBED)  #all reads passed
head(myBED)

#Get rid of read-length column
myBED <- myBED[, c(1:5) ]
head(myBED)

#-------------------------------------------------------------------
# Extend BED 500bp from its 3' end -> NO
#-------------------------------------------------------------------
forward <- which(myBED$V6 == '+')
reverse <- which(myBED$V6 == '-')
length(forward)
length(reverse)
length(forward) + length(reverse) == nrow(myBED) # check

#myBED$V3[forward] <- myBED$V3[forward] + 500
#myBED$V2[reverse] <- myBED$V2[reverse] - 500

#sanity check
length( which( myBED$V2 <= 0) )
myBED$V2 [  which( myBED$V2 <= 0)  ] <- 1


#ITERATE IN EACH CHROMOSOME:
head(genomeSizes)
UC <- unique(as.character(myBED$V1) )

for (j in length(UC)  ){
  temp1 <- which(as.character(myBED$V1) == UC[j] )
  
  temp2 <- which(   myBED$V3[temp1] > genomeSizes$V2 [ which(as.character( genomeSizes$V1) == UC[j]  )]     )
  if (length(temp2) > 0){
    myBED$V3[temp1][temp2] <- genomeSizes$V2 [ which(as.character( genomeSizes$V1) == UC[j]  )]  
  }
  
  rm(temp1, temp2)
  
}
rm(UC)

rm(reverse,forward)


#-------------------------------------------------------------------
# Assign GEMcode
#-------------------------------------------------------------------
colnames(myBED)  <- c('chr','start','end', 'readID', 'strand')
head(myBED)
dim(myBED)
colnames(mydata) <- c( 'GEM', 'readID' )

X <- merge(x=myBED, y=mydata, by='readID', sort = FALSE)
head(X)
dim(X)

#Get rid of 'mydata' and 'myBED'
rm(myBED)   #mydata was necessary for the Chr
rm(mydata)


#-------------------------------------------------------------------
# QC the barcodes. 
# RULE. Any barcode not matching the white list WILL BE removed as it will introduce mistakes
#-------------------------------------------------------------------
# Download 10X barcode white list
#https://github.com/OpenGene/fastp/issues/91
#wget https://raw.githubusercontent.com/10XGenomics/supernova/master/tenkit/lib/python/tenkit/barcodes/4M-with-alts-february-2016.txt
barcodeWhiteList <- fread('10X_barcodes_white_list_4M-with-alts-february-2016.txt', header=FALSE)
head(barcodeWhiteList )
dim(barcodeWhiteList )
length(unique(barcodeWhiteList$V1 )) #all are unique barcodes

write.table(paste0("Number of unique GEMcodes in 10X barcode white list: ", length(unique(barcodeWhiteList$V1 )) ), file = paste0(outdir,SampleID,'_GapWidth_',GapWidth,".log"),append = TRUE, quote = FALSE, sep = "\t", row.names = FALSE, col.names =FALSE)



if(QC_barcodes ==TRUE) {   # very slow. NOT TESTED
  #------------------------------------------------------------------------------------------------
  # Modify X$GEM. If hamming distance d=1 with only one 10X reference barcode, assign it
  #-------------------------------------------------------------------------------------------------
  library(dartR)
  hammingFun <- function(SEQ, b){  
  
  #d <- stringdist(a=SEQ, b=b, method ="hamming" ) ;
  #d <- as.numeric( adist(x=SEQ, y=b ) )
  d <- unlist(lapply(b, utils.hamming, str2=SEQ, r=0))

  minD <- min( d[which(d!=0)] ) 
  
  #if (minD != 1) { return(  SEQ   )   }
  if (minD == 1/16) { 
    LminD <- which(d == minD  )
    if (length(LminD) == 1 ) {  return( U[which(d == 1/16)] )  }
    #if (length(minD) > 1  )   return( SEQ )
    else { return(  SEQ   )   }
  }
  else { return(  SEQ   )   }
  
  rm(d, minD,LminD)
  }

  U <- unique(barcodeWhiteList$V1 )
  newBarcodes <- lapply(X$GEM, hammingFun , b=U)
  newBarcodes <- unlist (newBarcodes )
  #utils.hamming(str1="CGACACGGTTTGGGCA", str2="CGACACGGTTTGGGCC", r=0)
  #SEQ="CGACACGGTTTGGGCA"
  X$GEM <- newBarcodes 
  rm(U, newBarcodes)
#------------------------------------------------------------------------------------------------
}

 
library(fastmatch) 
MatchedBarcodes <- fmatch(x=X$GEM, table=barcodeWhiteList$V1, nomatch = NA_integer_, incomparables = NULL)

write.table(paste0("Percentage of GEMs with perfect match in the 10X white list: ", 100*(length( which(is.na(MatchedBarcodes)==FALSE) )/ nrow(X) ) ), file = paste0(outdir,SampleID,'_GapWidth_',GapWidth,".log"),append = TRUE, quote = FALSE, sep = "\t", row.names = FALSE, col.names =FALSE)
write.table("--" , file = paste0(outdir,SampleID,'_GapWidth_',GapWidth,".log"), append = TRUE, quote = FALSE, sep = "\t", row.names = FALSE, col.names =FALSE)

X <- X[which(is.na(MatchedBarcodes)==FALSE) ,]

rm(barcodeWhiteList, MatchedBarcodes)




#--------------------- not necessary IF RUNNING THE CODE ABOVE (it should be 0)----
### here we remove barcodes which will have an N
haveN <- unlist(lapply('N' , grepl, X$GEM ))  # grepl returns logical value
#proportion of GEMs with Ns in the mapped reads
100* ( length(which(haveN == TRUE))/nrow(X)  )

if( length(which(haveN == TRUE) ) > 0 ){
  X <- X[which(haveN == FALSE)  ,]
  head(X)
  dim(X)
}

rm(haveN)
#-----------------------------------------------------------------------------------


# Assign GEMcode
X$GEMcode <- unlist(lapply(X$GEM, GEMnumber))

#Get rid of unnecesary columns
X <- X[, c(2:5,7)]

#plot stats for the QC'ed bam file with mapped reads
TA <- table(X$GEMcode)


png(paste0(outdir ,GEMraw,'_GapWidth_',GapWidth,'_PROCESSED.png'), height=400, width=705)

  par(mfrow=c(1,2))
  let <- length(sort(TA, decreasing = TRUE) )
  plot( as.numeric(sort(TA, decreasing = TRUE) )[1:let ] , lwd=2, type='h', log='x', col='orange', ylab='Number of reads', xlab='',
      main='Frequency of GEM barcodes \n (Log x)', sub=paste0('\n Number of unique GEM barcodes: ', length(   unique(names(TA))   ), '\n Number of total reads after Processing: ', nrow(X))) 
  plot( as.numeric(sort(TA, decreasing = TRUE) )[1:let ] , lwd=2, type='h', log='xy', col='orange', ylab='Number of reads', xlab='',
      main='Frequency of GEM barcodes \n (Log xy)', sub=paste0('\n Number of unique GEM barcodes: ', length(   unique(names(TA))   ), '\n Number of total reads after Processing: ', nrow(X))) 
  rm(let)
dev.off()


write.table(paste0('Number of reads after Processing is....', nrow(X) ) , file = paste0(outdir,SampleID,'_GapWidth_',GapWidth,".log"), append = TRUE, quote = FALSE, sep = "\t", row.names = FALSE, col.names =FALSE)
write.table(paste0('Number of unique GEM barcodes is....', length(   unique(names(TA))  )) , file = paste0(outdir,SampleID,'_GapWidth_',GapWidth,".log"), append = TRUE, quote = FALSE, sep = "\t", row.names = FALSE, col.names =FALSE)
write.table("--" , file = paste0(outdir,SampleID,'_GapWidth_',GapWidth,".log"), append = TRUE, quote = FALSE, sep = "\t", row.names = FALSE, col.names =FALSE)


rm(TA)


write.table('Merging Reads if same GEMCode and closer than GapWidth to Create PUCs....', file = paste0(outdir,SampleID,'_GapWidth_',GapWidth,".log"), append = TRUE, quote = FALSE, sep = "\t", row.names = FALSE, col.names =FALSE)


#-----------------------------------------------------------------------------------
# Merge Read if same GEMCode and closer to GapWidth bp to Create CHROMATIN FRAGMENTS
#-----------------------------------------------------------------------------------

library(GenomicRanges)
library(Rsamtools)
library(bamsignals)

uniGEM <- unique(X$GEMcode)
head(uniGEM)
length(uniGEM)

# tapply function is the most convenient here to avoid for loop

Xtest <- X #[1:10000,]
Xtest


#genomic ranges
gr <-  GRanges(
  seqnames = Rle(Xtest$chr),
  ranges = IRanges(Xtest$start, end = Xtest$end ),
  strand = Rle(strand(c(Xtest$strand)) ) #,
  #GEM= rep('0', length(sel))
)
gr
mcols(gr)$GEM <-  Xtest$GEMcode #[1:10000] #[k]
gr

rm(Xtest)




MergeGenomicRanges <- function(bedData){
  
  #bedData is Genomic ranges with metadata colum GEM (barcode)
  
  grMerged <- reduce(x=bedData, drop.empty.ranges=FALSE, min.gapwidth= GapWidth, with.revmap=FALSE,  # merge if closer than 3kb
                     with.inframe.attrib=FALSE, ignore.strand=TRUE)
  
  # add GEMcode to GRanges
  U <- unique( mcols(bedData)$GEM )
  mcols(grMerged)$GEMcode <-  U
  
  #add read counts and keep only in ReadCounts pass a CutOff
  mcols(grMerged)$ReadCount <- bamCount(bampath=BAM, gr=grMerged , mapqual=0, verbose=FALSE, shift = 0) #bamsignals
  grMerged <- grMerged [  which(mcols(grMerged)$ReadCount >= MinReadsInEachFragment )  ]

  Le <- length(grMerged)
  
  if (Le >= 1) {
  
    #compute GEM purity index
    mcols(grMerged)$GEM_purity <- 0.0
    temp_chr   <- as.character(seqnames(grMerged ))
    temp_table <- table(temp_chr )
    #J          <- length( unique( temp_chr ) )
    for (h in 1: length( temp_table)){
      mcols(grMerged)$GEM_purity <-  mcols(grMerged)$GEM_purity +  ( as.numeric(temp_table[h])  / Le )^2
    }
    rm( temp_chr )
    #rm( temp_table )   WILL USE IT BELOW
    temp_sum <- sum (temp_table)
  

    #categorisation based on the number of chromatin fragments
    if( Le ==1  ) { 
      mcols(grMerged)$chromatin_fragment  <-  'singleton' 
      #mcols(grMerged)$GEM_purity <- -1    #    GEM purity for singletons = -1
      #mcols(grMerged)$multiplexType       <-  'NA' 
      #mcols(grMerged)$PUC       <-  U
      mcols(grMerged)$fragments  <-  1
      #
      mcols(grMerged)$multiplexType   <- 'none'  
    }
    if( Le >=2  ) { 
      mcols(grMerged)$chromatin_fragment <- 'multiplex' # 'multiplex-intra-chromosomal' 
      mcols(grMerged)$fragments  <-  length(grMerged)
      #now each multiplex fragment can be categorised into intra-chromosomal, inter-chromosomal  , or MIXED
      #CHRS <- unique(seqnames(grMerged)) # distinct chrs
      
      if( length(temp_table) == 1 ){  
          mcols(grMerged)$multiplexType   <- 'intra-chromosomal'
      #  mcols(grMerged)$PUC        <-  U
      #  mcols(grMerged)$fragments  <-  length(grMerged)
      }
      if( (length(temp_table) > 1) &  (length(temp_table) == temp_sum ) ){  
        mcols(grMerged)$multiplexType   <- 'inter-chromosomal'   
      #  for (k in 1:length(CHRS) ){
      #    mcols(grMerged)$PUC        <-  paste(U, CHRS[k] , sep='_' )
      #    mcols(grMerged)$fragments  <-  length( which( as.character(seqnames(grMerged) ) == as.character(CHRS) ) )
      #    
      #  }
      }
      if( (length(temp_table) > 1) &  (length(temp_table) != temp_sum ) ){  
        mcols(grMerged)$multiplexType   <- 'mixed'   
      }
        
    }
    
    rm( temp_table ) 
    rm( temp_sum  ) 
    return( as.data.table(grMerged , row.names = FALSE)[,c(1:4,6:11) ]   )  # DT object
  
    
  } 
}



mergedX <- tapply(gr, mcols(gr)$GEM,  MergeGenomicRanges  )


#There 2: should be equal
length(mergedX)
length( unique(mcols(gr)$GEM) ) # number of unique GEMs




mergedX2 <- rbindlist(mergedX, use.names=TRUE, fill=FALSE) #, idcol=TRUE)


unique(mergedX2$chromatin_fragment)
unique(mergedX2$multiplexType)


write.table('PUC Type:' , file = paste0(outdir,SampleID,'_GapWidth_',GapWidth,".log"), append = TRUE, quote = FALSE, sep = "\t", row.names = FALSE, col.names =FALSE)
write.table(table(mergedX2$chromatin_fragment) , file = paste0(outdir,SampleID,'_GapWidth_',GapWidth,".log"), append = TRUE, quote = FALSE, sep = "\t", row.names = FALSE, col.names =FALSE)
write.table('PUC multiplexType:' , file = paste0(outdir,SampleID,'_GapWidth_',GapWidth,".log"), append = TRUE, quote = FALSE, sep = "\t", row.names = FALSE, col.names =FALSE)
write.table(table(mergedX2$multiplexType) , file = paste0(outdir,SampleID,'_GapWidth_',GapWidth,".log"), append = TRUE, quote = FALSE, sep = "\t", row.names = FALSE, col.names =FALSE)

write.table("--" , file = paste0(outdir,SampleID,'_GapWidth_',GapWidth,".log"), append = TRUE, quote = FALSE, sep = "\t", row.names = FALSE, col.names =FALSE)
write.table('Distribution of Fragments in GEMs:' , file = paste0(outdir,SampleID,'_GapWidth_',GapWidth,".log"), append = TRUE, quote = FALSE, sep = "\t", row.names = FALSE, col.names =FALSE)
write.table(table(mergedX2$fragments)   /  as.numeric(names(table(mergedX2$fragments))) , file = paste0(outdir,SampleID,'_GapWidth_',GapWidth,".log"), append = TRUE, quote = FALSE, sep = "\t", row.names = FALSE, col.names =FALSE)
write.table("--" , file = paste0(outdir,SampleID,'_GapWidth_',GapWidth,".log"), append = TRUE, quote = FALSE, sep = "\t", row.names = FALSE, col.names =FALSE)



#Save PUTATIVE CHROMATIN COMPLEXES
#write.csv( mergedX2, file = paste0( 'chr',Chr, '_PUCs_', nrow(mergedX2),'_GEMs_',length(mergedX), ".csv" )  , row.names = FALSE)
#mergedX2$GEMcode <- as.integer( mergedX2$GEMcode )
#
#library(WriteXLS)
#WriteXLS(x=mergedX2, ExcelFileName = paste0( 'chr',Chr, '_PUCs_', nrow(mergedX2),'_GEMs_',length(mergedX), ".xls" ), SheetNames = NULL)
write.csv(x=mergedX2, file = paste0(outdir , SampleID, '_PUCs' ,"_GapWidth_",GapWidth, ".csv" ) , row.names = FALSE )


write.table(paste0("Number of distinct Fragments: ", nrow(mergedX2) ), file = paste0(outdir,SampleID,'_GapWidth_',GapWidth,".log"),append = TRUE, quote = FALSE, sep = "\t", row.names = FALSE, col.names =FALSE)
write.table(paste0("Number of distinct GEMs: ", length(mergedX) ), file = paste0(outdir,SampleID,'_GapWidth_',GapWidth,".log"),append = TRUE, quote = FALSE, sep = "\t", row.names = FALSE, col.names =FALSE)


rm(mergedX)


mergedX3 <- mergedX2 [  which( mergedX2$fragments >= FragmentNumberCutoff ) ,] 
head(mergedX3 )
dim( mergedX3 )


write.table(paste0("Number of Multiplex-Fragments: ", nrow(mergedX3) ), file = paste0(outdir,SampleID,'_GapWidth_',GapWidth,".log"),append = TRUE, quote = FALSE, sep = "\t", row.names = FALSE, col.names =FALSE)


rm(mergedX2)
#save for chia-view
cv <- data.frame(Chrom=paste0("chr",mergedX3$seqnames),	Start	=mergedX3$start,End=mergedX3$end,	Frag_num=mergedX3$fragments	,GEM_ID=paste0(SampleID,"-", mergedX3$GEMcode)  )
head(cv)
write.table(x=cv , file = paste0(outdir,  SampleID, "_GapWidth_",GapWidth, "_ChIA-view.region"), append = FALSE, quote = FALSE, sep = "\t",
            eol = "\n", na = "NA", dec = ".", row.names = FALSE,
            col.names = FALSE)
rm(cv)
##########################################

# finishing...

end_time <- format(Sys.time(), "%a %b %d %X %Y")

#print(paste0( ' ...and finished on ', end_time  ) )  
write.table(paste0("Analysis started at.... ", start_time  ), file = paste0(outdir,SampleID,'_GapWidth_',GapWidth,".log"),append = TRUE, quote = FALSE, sep = "\t", row.names = FALSE, col.names =FALSE)
write.table(paste0("Analysis has finished at ... ", end_time  ), file = paste0(outdir,SampleID,'_GapWidth_',GapWidth,".log"),append = TRUE, quote = FALSE, sep = "\t", row.names = FALSE, col.names =FALSE)


rm(start_time, end_time)

rm(mergedX3,X)
rm(uniGEM, gr)

write.table("--------------------------" , file = paste0(outdir,SampleID,'_GapWidth_',GapWidth,".log"), append = TRUE, quote = FALSE, sep = "\t", row.names = FALSE, col.names =FALSE)


