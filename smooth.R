
#source('~/Dropbox/Documents/src/simCrop.R')

#let's try to see what a single filter is like


#myFilter <- function(x) {
#
#  x.table <-table(x[!is.na(x)])
#  x.value <-names(x.table)
#  x.return <- x.value[which( x.table == max(x.table) )] 
#
#  if(length(x.return > 1) {
#    sample(x.return,size=1)
#  } else {
#    return(x.value)
#  }
#
#}
#
# This is what we want


# function to apply a modal response for each year in a spatial neighborhood defined by diameter s,
# based on classes defined within a temporal window, h.
localMode <- function( r, h, s) {

  a <- getValuesFocal(r,ngb=s)
  
  # a is a list where each element is pixels x focal
  # so what we want to do is to grab our first set and find what is the most popular
  a.return.values <- matrix(0, nrow=ncell(r), ncol=length(a))

  for(i in 1:ncell(r)) { 
  #for(i in 1:1000) { 

    a.hood <- c()
    for(j in 1:length(a)) {
    
      a.hood <- rbind(a.hood, a[[j]][i,] ) 
    
    }
  
    # run the filter over each year
    for( j in 1:length(a)) {

      # get the index of our response value
      h0 <- min(j,h+1)
  
      
      #only get years within the window 
      time.window <- max(1, j - h):min( length(a), j + h) 
  
      
      # so this mess really just gets the non NA neighborhood, and classes within the spatial domain together.
      # we transpose each of these classes so we can use aggregate to determine the most popular classes
      a.work <- data.frame(cbind(1, t(a.hood[time.window,!is.na(a.hood[1,])])))
      
      a.results <- aggregate( X1 ~ . , a.work, FUN=length)
      
      a.freq <- a.results$X1
      
      a.return <- a.results[which( a.freq == max(a.freq) ),1:length(time.window)] 
      
      if(nrow(a.return) > 1) {
        a.return <- a.return[sample(1:nrow(a.return),size=1),]
      }
  
      a.return.values[i,j] <- a.return[1,h0]
    }
  }

  return(a.return.values)

}
  


runTime <- proc.time()

h <- 3 
s <- 11

smoothValues <- localMode(r,h,s)

print("runtime:")
print(proc.time() - runTime)

x11()
r.new <- r
values(r.new) <- smoothValues
title(sprintf("h=%d, s=%d",h,s))


