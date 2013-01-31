
source('~/Dropbox/Documents/src/simCrop.R')

#let's try to see what a single filter is like


myFilter <- function(x) {

  x.table <-table(x[!is.na(x)])
  x.value <-names(x.table)
  x.return <- x.value[which( x.table == max(x.table) )] 

  if(length(x.return > 1) {
    sample(x.return,size=1)
  } else {
    return(x.value)
  }

}

rf1 <- focal(r1,w=3,fun=myFilter  )




