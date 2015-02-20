

## Load data
dat <- read.table('gene_count_data.txt',header=T,row.names=1)
len <- read.table('isoform_lengths.txt',header=F,row.names=1)



## keep columns
dat <- dat[,c(1,2)]


## merge read counts and lengths
input.dat <- merge(len,dat,by='row.names')


## Produce normalized data
output.tpm <- rnaseqnorm(input.dat[,2:dim(input.dat)[2]],input.dat[,1])




## normalization function
 
rnaseqnorm <- function(reads,lengths) {
    mat <- as.matrix(reads)

    lengths <- lengths / 1000 # convert to kb

    reads.rpkm <- mat
    reads.tpm <- mat

    libsum <- apply(mat,2,sum)

    for( j in 1:dim(mat)[2] ) {
        w <- cbind(mat[,j],lengths)

        v <- apply(w,1, function(x) x[1]/x[2] )

        # 75th percentile
        #quan <- quantile(v,0.75)
        quan <- quantile(v,1)

        altsum <- sum( v[which(v <= quan)] )
        for ( i in 1:dim(mat)[1] ) {
            reads.rpkm[i,j] <- ( mat[i,j] / lengths[i] ) * 1000000 / libsum[j]
            reads.tpm[i,j] <- ( mat[i,j] / lengths[i] ) * 1000000 / altsum
        }
    }

#   print( reads.rpkm )
#   print( reads.tpm )

    #return( reads.rpkm )
    return( reads.tpm )

}





