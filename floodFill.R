# this is a program to sequentially open a raster file row by row


library(raster)

## FUNCTIONS ##

### example to use as input
#set.seed(0)
#x <- matrix( ceiling(4 * runif(12)), nrow=4, ncol=3 )
#x[2,2] <- 4
#r <- raster(x)


floodFill.raster <- function(r) {

  # constants
  nrow_r <- nrow(r)
  ncol_r <- ncol(r)
  ncell_r <- ncell(r)
  
  v <- as.matrix(getValues(r))   # values of the raster image
  v.obj <- matrix(0,nrow=ncell(r),1) # ownership values we are making
  
  # get the initial position and derived values
  startPoint <- 1
  objNo <- 1
  
  while ( startPoint != 0 ) {
    
    Q <- startPoint
    v.current <- v[startPoint,] # our current target value
    
    # we iterate until the queue is empty
    while (length(Q) != 0) {
    
      k.current <- Q[1]  # get our initial position
      Q <- Q[-1]         # update the queue
    
      # set our initial object number
      
      ## get the max in left and right
    
      # first we need to reset k.min.current and k.min.current
      k.min.current <- 0
      k.max.current <- 0
  
      k.row.current <- ceiling(k.current / ncol_r)
      
      # left
      k <- k.current
      while( ceiling(k/ncol_r) == k.row.current    )  {
            if( identical( v[k,], v.current ) && (v.obj[k] == 0) )  { 
              k.min.current <- k 
            } else { 
              break 
            }
            
            # update values
            k <-  k - 1
      #      j <- ceiling( k / nrow_r )    #column 
      #      i <-  ((k -1) %% nrow_r ) + 1 #row
      }
      
      # right
      k <- k.current
      while( ceiling(k/ncol_r) == k.row.current    )  {
            if( identical( v[k,], v.current ) && (v.obj[k] == 0) ) { 
              k.max.current <- k 
            } else { 
              break
            }
      
            k <-  k + 1
      }
      
      # add our classified values if there are any to add
      if( k.min.current > 0 ) v.obj[k.min.current:k.max.current] <- objNo
   
    
      Q.high <- 
          (k.min.current:k.max.current + ncol_r)[
            ((k.min.current:k.max.current + ncol_r) <= ncell_r) &
            ((k.min.current:k.max.current) > 0)
          ]
    
      Q.low <- 
          (k.min.current:k.max.current - ncol_r)[
            ((k.min.current:k.max.current - ncol_r) >= 0) &
            ((k.min.current:k.max.current) > 0)
          ]
       
      #enqueue the new values
      Q <- c(Q,Q.high,Q.low)
    }
  
    # figure out the new start point
    startPoint <- 0 
    for(i in 1:ncell_r) 
      if( v.obj[i] == 0) {
        startPoint <- i
        break
      }
      
    # move to a new object
    objNo <- objNo + 1
  }

  return( v.obj )
}
