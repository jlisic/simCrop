# simple metropolis hastings implementation
library(mvtnorm)
library(Matrix)


mh <- function(foo,x0,iter,burnIn=50, mu = 0) {

  results <- 1:iter

  x <- x0
  # iterations for metropolis hastings algorithm
  for(i in 1:(iter + burnIn) ) {
    y <- rnorm(1,mean=x)

    u <- runif(1)

    rho <- min( foo(y) * dnorm(x,mean=y + mu) / ( foo(x) * dnorm(y,mean=x + mu) ), 1 )
    #rho <- min( foo(y) /  foo(x), 1 ) 

    if (u < rho) x <- y 

    if( i > burnIn) results[i - burnIn] <- x
  }

  return( results )
}

#  --- lambda ---

mh.lambda <- function(Z,W,x0,iter,burnIn=50, lambda.range) {
    
  n <- nrow(W)

  results <- 1:iter

  x <- x0
  sigma.inv.x <- diag(n) - x * W 
  det.sigma.x <- 1/det(sigma.inv.x) 

  # iterations for metropolis hastings algorithm
  for(i in 1:(iter + burnIn) ) {
    y <- runif(1, min=lambda.range[1],max=lambda.range[2]) 
    sigma.inv.y <- diag(n) - y * W 
    det.sigma.y <- 1/det(sigma.inv.y) 

    u <- runif(1)

    #rho <- min( foo(y) /  foo(x), 1 ) 
    rho2 <-  sqrt(det.sigma.y/det.sigma.x ) * exp( -1/2 * ( t(Z) %*% sigma.inv.y %*% Z - t(Z) %*% sigma.inv.x %*% Z ) )  
    rho <- min( rho2, 1 )

    if (u < rho) {
      x <- y
      sigma.inv.x <- sigma.inv.y
      det.sigma.x <- det.sigma.y 
    }


    if( i > burnIn) results[i - burnIn] <- x
  }

  return( results )
}


#  --- lambda ---

mh.lambda.sar <- function(Z,W,mu,x0,iter,burnIn=50, lambda.range) {
    
  n <- nrow(W)
  N <- nrow(Z)/n

  results <- 1:iter

  x <- x0
  lambda.inv.x <- diag(n) - x * W 
  det.lambda.x <- 1/det(lambda.inv.x) 
  Z.adj.x <- matrix(lambda.inv.x %*% matrix(Z,nrow=n),ncol=1) - mu 

  # iterations for metropolis hastings algorithm
  for(i in 1:(iter + burnIn) ) {
    y <- runif(1, min=lambda.range[1],max=lambda.range[2]) 
    lambda.inv.y <- diag(n) - y * W 
    det.lambda.y <- 1/det(lambda.inv.y) 
    Z.adj.y <- matrix(lambda.inv.y %*% matrix(Z,nrow=n),ncol=1) - mu 
 
    u <- runif(1)

    det.part <- ( (det.lambda.y) / (det.lambda.x) )^(N)

    #rho <- min( foo(y) /  foo(x), 1 ) 
    rho2 <- det.part * exp( -1/2 * ( t(Z.adj.y) %*% Z.adj.y  - t(Z.adj.x) %*%  Z.adj.x ) )  
    rho <- min( rho2, 1 )

    if (u < rho) {
      x <- y
      lambda.inv.x <- lambda.inv.y
      det.lambda.x <- det.lambda.y 
      Z.adj.x <- Z.adj.y
    }


    if( i > burnIn) results[i - burnIn] <- x
  }

  return( results )
}


# Z is n by t


lambda.propto.post <- function(lambda,Z , X, Beta, W) {

  lambda.support <- sort( 1/range(eigen(W)$values) )

  if( lambda.support[1] > lambda) return( 0 )
  if( lambda.support[2] < lambda) return( 0 )

  post <- 1

  for( i in 1:ncol(Z) ) {
    post <- post * dmvnorm(Z[,i],mean= W %*% Z[,i] + X %*% Beta, sigma= solve( diag(ncol(W)) - lambda * W ))
  }

  return(post)
}




# -- Test -- 
if ( F ) {

L <- 100 
W <- as.matrix(bandSparse(L,L,k=list(-1*sqrt(L),-1,1,sqrt(L))))*1

sigma <- solve(diag(L) - 0.25 * W)
sigma.eigen <- eigen(sigma)
sigma.half <- sigma.eigen$vectors %*% diag(sqrt(sigma.eigen$values)) 

deviates.mean <- c()
lambda.ls <- c()


for( j in 1:1000) {

Z <- sigma.half %*% rnorm(L)
X <- matrix(1,nrow=L,ncol=1)
B <- 0

lambda.range <- sort( 1/range(eigen(W)$values) )

# this is a function used to get a look at what the posterior distribution looks like
#l <- seq( lambda.range[1], lambda.range[2], length.out=102)
#l <- l[2:101]
#dl <- 1:100
#
#for( i in 1:100 ) {
#  print(i)
#  dl[i] <- lambda.propto.post( l[i], Z, X, B, W)
#}


# test the mh on the posterior

  Y <- Z - mean(Z)

  deviates <- mh.lambda(Y,W,0,1,100,lambda.range )

  deviates.mean[j] <- mean(deviates)


  lambda.ls[j] <- t(Y) %*% W %*% Y / (t(Y) %*% W %*% W %*% Y)
}

}
