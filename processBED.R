setwd('/chiadrop/sampleX')

library(data.table)

mydata <- fread('q40_sampleX_sorted_filtered.bed', header=FALSE)

#ONLY APPLY THIS TO PE read alignments
temp <- tstrsplit(x=mydata$V4, split="/", fixed=TRUE)[[1]]
mydata$V4 <- temp

mydata2 <- mydata [ which(mydata$V1 =='1' | mydata$V1 =='2' | mydata$V1 =='3' |  mydata$V1 =='4' |
                    mydata$V1 =='5' | mydata$V1 =='6' | mydata$V1 =='7' |  mydata$V1 =='8' |
                    mydata$V1 =='9' | mydata$V1 =='10' | mydata$V1 =='11' |  mydata$V1 =='12' |
                    mydata$V1 =='13' | mydata$V1 =='14' | mydata$V1 =='15' |  mydata$V1 =='16' |
                    mydata$V1 =='17' | mydata$V1 =='18' | mydata$V1 =='19' |  
                    mydata$V1 =='X' | mydata$V1 =='Y' )  , ]

# keep reads with R1 and R2
i1 <- which( duplicated(mydata2$V4) ==FALSE )
i2 <- which( duplicated(mydata2$V4) ==TRUE )

y = merge(x=mydata2[i1,], y=mydata2[i2,], by="V4", all=FALSE)

#Combine R1 and R2

mydata3 <- data.table(
    V1 <- y$V1.x ,
    V2 <-  apply(y[,c(3,4,14,15)], 1, FUN=min) ,
    V3 <-  apply(y[,c(3,4,14,15)], 1, FUN=max) ,
    V4 <- y$V4 ,
    V5 <- 0 ,
    V6 <- '+' ,
    V7 <- 68 ,
    V8 <- 188,
    V9 <- NA ,
    V10 <- 1 ,
    V11 <- 120,
    V12 <- 0
    ) 

mydata3$V5 <- mydata3$V3 - mydata3$V2

fwrite(x=mydata3, file = "sampleX.bed", row.names=FALSE, col.names=FALSE, sep="\t")

