minDim <- function(x){
  n <- nrow(x)
  if( is.null(n) ) {
    return(1)
  } 
  
  if( ncol(x) < n ) {
    return( ncol(x) ) 
  } 

  return(n)
}
