# this is a program to sequentially open a raster file row by row

set.seed(0)

library(raster)

## FUNCTIONS ##

## we store our results in a list of c( VALUE, STOPINDEX)
#getStarts <- function( x ) {
#  for( j in 1:length( x ) ) {
#    # get first row
#    if( j == 1) {
#      k = 1
#      a <- list( c(x[j],j) )
#    # for any other row that differs in category
#    } else if ( a[[k]][1] != x[j] )  {
#      k = k + 1 
#      a[k] <- list( c(x[j],j) )
#    }
#  }
#  return(a)
#} 
#
#
#
#
#
#
### example to use as input
x <- matrix( ceiling(4 * runif(12)), nrow=4, ncol=3 )
r <- raster(x)
#
#i <- 1
#oneRow <- getValues(r,i)
#
#b.old <- getStarts( oneRow )
#
#i = i + 1
#oneRow <- getValues(r,i)
#
#b.current <- getStarts( oneRow )


nrow_r <- nrow(r)
k <- 1

Q <- c()
# remember we are in a row major form 
j <- ceiling( k / nrow(r) ) 
i <-  ((k -1) %% nrow(r) ) + 1

v.current <- v[1] # our current target value
k.current <- k

# get the max in left and right
k <- k.current
while( j >= 1) {
      if( v[k] == v.current )  { 
        k.min.current <- k 
      } else { 
        break 
      }
      
      #check down
      if( i + 1 <= ncol_r ) 
        if( v[k + nrow_r] == v.current ) {
          Q <- c(Q,k + nrow_r)
        }

      #check up
      if( i - 1 <= ncol_r ) 
        if( v[k - nrow_r] == v.current ) {
          Q <- c(Q,k - nrow_r)
        }

      # update values
      k <-  k - 1
      j <- ceiling( k / nrow_r )    #column 
      i <-  ((k -1) %% nrow_r ) + 1 #row
}

k <- k.current
while( j <= nrow_r ) { 
      if (v[k] == v.current) & (j <= nrow_r ) { 
        k.max.current <- k 
      } else { 
        break
      }

      #check down
      if( i + 1 <= ncol_r ) 
        if( v[k + nrow_r] == v.current ) {
          Q <- c(Q,k + nrow_r)
        }

      #check up
      if( i - 1 <= ncol_r ) 
        if( v[k - nrow_r] == v.current ) {
          Q <- c(Q,k - nrow_r)
        }



      k <-  k + 1
      j <- ceiling( k / nrow_r )    #column 
      i <-  ((k -1) %% nrow_r ) + 1 #row
}

# now let's go up and down
