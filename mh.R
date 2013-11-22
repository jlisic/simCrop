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

mh.lambda.sar <- function(Z,W,mu,tau,x0,iter,burnIn=50, rho.range) {
    
  n <- nrow(W)
  K <- nrow(Z)/n

  results <- 1:iter

  x <- x0
  Lambda.x <- diag(n) - x * W
  
  Z <- matrix(Z,nrow=n)
  mu <- matrix(mu,nrow=n)
  

  Z.adj.x <- Lambda.x %*% Z - mu 

  # iterations for metropolis hastings algorithm
  for(i in 1:(iter + burnIn) ) {
    y <- runif(1, min=rho.range[1],max=rho.range[2]) 
    Lambda.y <- diag(n) - y * W
    Z.adj.y <- Lambda.y %*% Z - mu 
 
    u <- runif(1)

    #rho <- min( foo(y) /  foo(x), 1 ) 
    quadratic.form <-  sum( Z.adj.y*Z.adj.y  - Z.adj.x*Z.adj.x )   
    
    if( quadratic.form < 0 ) { 
      acceptNewProb <- 1
    } else {
      rho2 <- exp( -1/2 * tau * sum( Z.adj.y*Z.adj.y  - Z.adj.x*Z.adj.x ) )  
      
      if( is.nan(rho2) ) { 
        acceptNewProb <- 0
        acceptNewProb <- min( rho2, 1 )
      } else {
        acceptNewProb <- min( rho2, 1 )
      }
    }

    if (u < acceptNewProb ) {
      x <- y
      Lambda.x <- Lambda.y
      Z.adj.x <- Z.adj.y
    }


    if( i > burnIn) results[i - burnIn] <- x
  }

  return( results )
}

# this function performs mh sampler on either rho1 or rho2, 
mh.q.sar <- function(rho,Z,W,Beta,q.value,iter,burnIn=50, rho.range) {
    
  n <- nrow(W)
  N <- nrow(Z)
  K <- N/n
   
  rho1 <- rho[1]
  rho2 <- rho[2]

  M <- matrix(1, nrow=K,ncol=K) 

  Lambda1 <- diag(n) - rho1 * W
  Lambda2 <- diag(n) - rho2 * W

  S1 <- kronecker( diag(K), solve(Lambda1) )

  Sigma1 <- t(S1) %*% S1
  Sigma2 <- kronecker( M, solve( t(Lambda2) %*% Lambda2 ) )

  Sigma.x <- q.value*Sigma1 + (1-q.value)*Sigma2
  Sigma.inv.x <- solve(Sigma.x) 
  Sigma.det.x <- det(Sigma.x)

  Q.x <- sqrt(q.value) * S1 %*% X 
    
  Z.adj.x <-  Z - Q.x %*%Beta 

  results <- 1:iter

  # set the initial value for what we want
  x <- q.value 

  # iterations for metropolis hastings algorithm
  for(i in 1:(iter + burnIn) ) {
    y <- runif(1) 
  
    Sigma.y <- y*Sigma1 + (1-y)*Sigma2
    Sigma.inv.y <- solve(Sigma.y) 
    Sigma.det.y <- det(Sigma.y)
    Q.y <- sqrt(y) * S1 %*% X 
    Z.adj.y <-  Z - Q.y %*%Beta 
    
    u <- runif(1)

    det.part <- sqrt( (Sigma.det.x) / (Sigma.det.y) )

    if (is.finite(det.part) ) {
      #rho <- min( foo(y) /  foo(x), 1 ) 
      alpha <- det.part * exp( -1/2 * ( t(Z.adj.y) %*% Sigma.inv.y %*%  Z.adj.y  - t(Z.adj.x) %*% Sigma.inv.x %*%  Z.adj.x ) )  
      alpha <- min( alpha, 1 )
    } else {
      alpha <- 0 
    }

    det.part <<- det.part

    if (u < alpha) {
      x <- y
      Sigma.inv.x <- Sigma.inv.y
      Sigma.det.x <- Sigma.det.y 
      Z.adj.x <- Z.adj.y
    }


    if( i > burnIn) results[i - burnIn] <- x
  }

  return( results )
}



# this function performs mh sampler on either rho1 or rho2, 
mh.lambda.sar2 <- function(rho.which,rho,Z,W,Beta,q.value,iter,burnIn=50, rho.range) {
    
  n <- nrow(W)
  N <- nrow(Z)
  K <- N/n
   
  rho1 <- rho[1]
  rho2 <- rho[2]

  M <- matrix(1, nrow=K,ncol=K) 

  Lambda1.x <- diag(n) - rho1 * W
  Lambda2.x <- diag(n) - rho2 * W

  S1.x <- kronecker( diag(K), sqrt(q.value) * solve(Lambda1.x) )

  Sigma1.x <- t(S1.x) %*% S1.x
  Sigma2.x <- (1-q.value) * kronecker( M, solve( t(Lambda2.x) %*% Lambda2.x ) )

  Sigma.x <- Sigma1.x + Sigma2.x
  Sigma.inv.x <- solve(Sigma.x) 
  Sigma.det.x <- det(Sigma.x)

  Q.x <- S1.x %*% X 
    
  Z.adj.x <-  Z - Q.x %*%Beta 
  Z.adj.y <- Z.adj.x # quick hack to help out the rho2 case 

  results <- 1:iter

  # set the initial value for what we want
  x <- rho[rho.which]

  # iterations for metropolis hastings algorithm
  for(i in 1:(iter + burnIn) ) {
    y <- runif(1, min=rho.range[1],max=rho.range[2]) 
    
    if(rho.which == 1) {
      Lambda1.y <- diag(n) - y * W
      S1.y <- kronecker( diag(K), sqrt(q.value) * solve(Lambda1.y) )
      Sigma1.y <- t(S1.y) %*% S1.y
      Sigma.y <- Sigma1.y + Sigma2.x
      Q.y <- S1.y %*% X 
      Z.adj.y <-  Z - Q.y %*%Beta 

    } else {
      Lambda2.y <- diag(n) - y * W 
      Sigma2.y <- (1-q.value) * kronecker( M, solve( t(Lambda2.y) %*% Lambda2.y) )
      Sigma.y <- Sigma1.x + Sigma2.y
    }
   

    Sigma.inv.y <- solve(Sigma.y) 
    Sigma.det.y <- det(Sigma.y) 
    
    u <- runif(1)

    det.part <- sqrt( (Sigma.det.x) / (Sigma.det.y) )

    #rho <- min( foo(y) /  foo(x), 1 ) 
    alpha <- det.part * exp( -1/2 * ( t(Z.adj.y) %*% Sigma.inv.y %*%  Z.adj.y  - t(Z.adj.x) %*% Sigma.inv.x %*%  Z.adj.x ) )  
    alpha <- min( alpha, 1 )

    if (u < alpha) {
      x <- y
      Sigma.inv.x <- Sigma.inv.y
      Sigma.det.x <- Sigma.det.y 
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
