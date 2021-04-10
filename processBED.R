setwd('/chiadrop/SAMPLEX')

library(data.table)
mydata <- fread('q40_SAMPLEX_sorted_filtered.bed', header=FALSE)
head(mydata)
dim(mydata)

#ONLY APPLY THIS TO PE read alignments
temp <- tstrsplit(x=mydata$V4, split="/", fixed=TRUE)[[1]]
mydata$V4 <- temp
rm(temp)



sort(table(mydata$V1))

# chrs for test
#mydata2 <- mydata [ which(mydata$V1 =='MT' | mydata$V1 =='GL456221.1' |   mydata$V1 =='JH584304.1' | 
#                           mydata$V1 == 'GL456216.1' | mydata$V1 == 'GL456392.1')  , ]

#mydata2 <- mydata [  which(mydata$V1 =='16' | mydata$V1 =='17' | mydata$V1 =='18' |   mydata$V1 =='19')  , ]


mydata2 <- mydata [ which(mydata$V1 =='1' | mydata$V1 =='2' | mydata$V1 =='3' |  mydata$V1 =='4' |
                    mydata$V1 =='5' | mydata$V1 =='6' | mydata$V1 =='7' |  mydata$V1 =='8' |
                    mydata$V1 =='9' | mydata$V1 =='10' | mydata$V1 =='11' |  mydata$V1 =='12' |
                    mydata$V1 =='13' | mydata$V1 =='14' | mydata$V1 =='15' |  mydata$V1 =='16' |
                    mydata$V1 =='17' | mydata$V1 =='18' | mydata$V1 =='19' |  
                    mydata$V1 =='X' | mydata$V1 =='Y' )  , ]

head(mydata2)
dim(mydata2)

sort(table(mydata2$V1))



######################################################################################################
rm(mydata)
head(mydata2)
dim(mydata2)
# keep reads with R1 and R2

i1 <- which( duplicated(mydata2$V4) ==FALSE )
length(i1)
i2 <- which( duplicated(mydata2$V4) ==TRUE)
length(i2)

y = merge(x=mydata2[i1,], y=mydata2[i2,], by="V4", all=FALSE)
head(y)
dim(y)

#Combine R1 and R2
rm(mydata2)

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

head(mydata3)
dim(mydata3)


mydata2 <- mydata3
######################################################################################################


fwrite(x=mydata2, file = "test.bed", row.names=FALSE, col.names=FALSE, sep="\t")

