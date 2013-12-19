# the goal here is to plot with SARMNP

Beta <- c( 
           1.0,  1.0,
          -1.0, -1.0,
           0.5,  0.5
          )

x <- rnorm(10)
y <- rnorm(10)
d <- round(runif(10) + runif(10)) + 1

z <- c(t(cbind(x,y)))

J <- 2


#function ( Z, J, Beta ) {

  Z <- t(matrix(z,nrow=J))
  Beta <- t(matrix(Beta,nrow=J))

  n.beta <- length(Beta) / J

  par(mfrow=c(1,n.beta))

  for(k in 1:n.beta) { 

    Z <- matrix(z,ncol=J)
    a.max <- max( abs(z) ) 
    
    plot(c(),
         xlim=c(-1 * a.max, a.max), 
         ylim=c(-1 * a.max, a.max),
         main=sprintf("Beta = %d",k),
         xlab='x',
         ylab='y' 
         )
    
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
   
    points(x,y)

    # add beta
    points(x=Beta[k,1], y=Beta[k,2], col="blue", pch=19, cex=3.0) 
    text(x=Beta[k,1], y=Beta[k,2], col="blue", label=sprintf("Beta = (%2.3f, %2.3f)", Beta[k,1], Beta[k,2] ), cex=3.0,pos= ) 

  }

#}
