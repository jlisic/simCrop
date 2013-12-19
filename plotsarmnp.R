

plotSARMNP <- function ( b, Beta=NULL, plotPoints=T,removeInitial=F ) {
  require(MASS)

  J <- length(b$crops) - 1
  Ncol <- ncol(b$cropValue)

  # get the values
  Z <- t(matrix( b$cropValue[,-1], nrow=J) )
 
  # get the transitions
  if( removeInitial ) {
    Y <- c( b$cropType[,c(-1,-2,-1*(Ncol+1))] )
  } else {
    Y <- c( b$cropType[,c(-1,-1*(Ncol+1))] )
  }

  if( length(Z) != J*length(Y)) stop(sprintf("values and J * prior states do not match (%d,%d,%d)", length(Z), J, length(Y)))

  if( !is.null(Beta)) {
   Beta <- t(matrix(Beta,nrow=J))
   n.beta <- length(Beta) / J
  }

  par(mfrow=c(1,J+1))

  for(k in 1:(J+1)) { 
    a.max <- max( abs(c(Z)) ) 
    x <- Z[Y==k,1]
    y <- Z[Y==k,2]

    if( !is.null(Beta) ) {
    plot(c(),
         xlim=c(-1 * a.max, a.max), 
         ylim=c(-1 * a.max, a.max),
         main=sprintf("J = %d, Beta = (%2.3f, %2.3f)", k, Beta[k,1], Beta[k,2] ), 
         xlab='x',
         ylab='y' 
         )
    } else {
    plot(c(),
         xlim=c(-1 * a.max, a.max), 
         ylim=c(-1 * a.max, a.max),
         main=sprintf("J = %d",k),
         xlab='x',
         ylab='y' 
         )
    } 

    # add polygons
    a <- matrix( 
           c(
             0         , 0,
             a.max     , a.max,
             -1 * a.max, a.max, 
             -1 * a.max, 0, 
             0         , 0
             ),byrow=T,ncol=2)
    polygon(a,col='lightblue')
    polygon(a[,c(2,1)],col='yellow')
   
    # plot x and y axis
    a <- matrix( c(-1 * a.max, 0, a.max,0),byrow=T,nrow=2)
    lines(a,col='red')
    lines(a[,c(2,1)],col='red')
  
    if( plotPoints ) points(x,y)


    myDensity <- kde2d(x,y)
    contour(myDensity,levels=(0:50)/50,add=T )

    # add beta
    if( !is.null(Beta)) {
      points(x=Beta[k,1], y=Beta[k,2], col="blue", pch=19, cex=3.0) 
    }


  }

}



