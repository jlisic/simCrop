########################################## LIBRARY ###############################################
library(msm)  # for rtnorm
library(mvtnorm)
library(Matrix)

########################################## INCLUDES ###############################################
source('mh.R')

########################################## FUNCTION ###############################################


## side effect addition e.g. ++
`%+=%` = function(e1,e2) eval.parent(substitute(e1 <- e1 + e2))
x = 1
x %+=% 2 ; x


## function that simulates a sample of size m, using a transition matrix theta
simSample <- function( m, theta) {
  currentState <- c(1,0)
  rotate <- c()
  for(i in 1:m) {
  
    currentProb <-currentState %*% theta  
    trial <- rbinom(1,size=1,prob=1-currentProb[1]) # probability of soybeans
    currentState <- c(0,0)
    currentState[trial + 1] <-1
    rotate[i] <- which( currentState == 1) - 1
  }

  return(rotate)
}


## create function here ## 
probitGibbs <- function(y,X,Beta,Beta0,Sigma0,iter) {

  #init some values
  m <- length(y)
  Beta.save <- c() 
  z <- 1:m

  # the usual X'X 
  XX <- t(X) %*% X

  # inverse of the prior variance
  B.star.inv <- solve(Sigma0^2)
  B <- (B.star.inv + XX )^-1

  # the mcmc loop
  for(i in 1:iter) {
    
    # generate deviates for the latent variables
    for( j in 1:m) {
      if( y[j] == 1) {
        z[j] <- rtnorm( 1, mean=X[j,] %*% Beta, sd=1 , lower=0) 
      } else {
        z[j] <- rtnorm( 1, mean=X[j,] %*% Beta, sd=1 , upper=0) 
      }
    }
  
    # this is Albert / Chib Beta sqiggle 
    Beta.post <- B %*% (B.star.inv %*% Beta0 + t(X) %*% z) 
  
    # generate deviates for beta/mu
    Beta.save[i] <- rnorm(1,Beta.post, B )
    Beta <- Beta.save[i]
  }
  
  return(Beta.save)
}


########################################## CONSTANTS ###############################################


## simulated sample size 
m <- 200

## here are the true parameters
## to:   corn soy
theta <- matrix( c(.20, .80, # from: corn     
                   .98, .02) #       soy
         , byrow=T, ncol=2)


########################################## SCRIPT ###############################################

if( F ) {
#initial state
set.seed(0)

# I assume there exists a 

Beta0.corn <- 0 
Beta0.soy <- 0 

Sigma0.corn <- 10
Sigma0.soy <- 10

# note we need to find some better way to associate x with y, it doesn't matter right now 
# since x is just ones
rotate <- simSample( m, theta) 
rotate.mat <- matrix(0,nrow=2,ncol=2)
for( i in 1:( length(rotate) -1 ) ) rotate.mat[ rotate[i] + 1, rotate[i+1] +1] %+=% 1 

y.corn <- matrix( c(
                    rep(0,times=rotate.mat[1,1]), 
                    rep(1,times=rotate.mat[1,2])
                   ),
                    ncol=1
                   ) 

X.corn <- matrix(1,nrow=length(y.corn),ncol=1)

Beta.init.corn <- solve(t(X.corn) %*% X.corn) %*% t(X.corn) %*% y.corn 


# run the results for corn
Beta.Gibbs <- probitGibbs(y.corn, X.corn, Beta.init.corn, Beta0.corn, Sigma0.corn, 1000)
}




## create function here ## 
# Y vector of categorical responses in row major form, repeating for each year, length = (number of years) x fieldSize
# X matrix of covariates  in row major form, repeating for each year, length = (number of years) x fieldSize
# W matrix in row major form of spatial neighborhoods, dim is fieldSize x fieldSize
# fieldSize, number of observations in a given year
#  
probitGibbsSpatial <- function(Y,X,W,fieldSize,Beta.init,lambda.init,Beta0,Sigma0,iter) {

  # set initial conditions
  Beta <- Beta.init
  lambda <- lambda.init 

  #init some values
  m <- nrow(Y)     # number of observations over all years
  K <- m/fieldSize # the number of years
  p <- ncol(X)     # number of covariates

  Beta.save <- c()
  lambda.save <- c()

  W.big <- kronecker(diag(K),W) 
  Sigma.inv <- diag(m) - lambda * W.big

  lambda.range <- sort( 1/range(eigen(W)$values) )

  Z <- matrix(0,nrow=m,ncol=1) 

  # the usual X'X
  XX <-  t(X) %*% Sigma.inv %*% X
 
  # inverse of the prior variance
  B.star.inv <- solve(Sigma0^2)
  B <- (B.star.inv + XX )^-1

  # the mcmc loop
  for(i in 1:iter) {
 
    means <- X %*% Beta + lambda * W.big %*% (Z - X %*% Beta )

    # generate deviates for the latent variables
    for( j in 1:m) {
      if( Y[j] == 1) {
        Z[j] <- rtnorm( 1, mean=means[j], sd=1 , lower=0) 
      } else if( Y[j] == 0) {
        Z[j] <- rtnorm( 1, mean=means[j], sd=1 , upper=0) 
      } else {
        Z[j] <- rnorm( 1, mean=means[j], sd=1) 
      }
      
      means <- X %*% Beta + lambda * W.big %*% (Z - X %*% Beta )
    }

    if( T ) { 
    #XZ <- t(X) %*% (Z - lambda * W.big %*% (Z - X %*% Beta) ) 
    XZ <- t(X) %*% Sigma.inv %*%  Z  
  
    # this is Albert / Chib Beta sqiggle 
    Beta.post.location <- B %*% (B.star.inv %*% Beta0 + XZ )  
    
    # generate deviates for beta/mu
    Beta.save[i] <- rnorm(1,Beta.post.location, B )
    Beta <- Beta.save[i]
    }

    # generate lambda deviate
    if( T ) {
    lambda.save[i] <- mh.lambda(Z - X %*% Beta,W,0,1,100,lambda.range )
    lambda <- lambda.save[i] 
    }
  }
  
  return( list( Beta = Beta.save, lambda = lambda.save) )
}

test1 <- T
test2 <- F


if ( test1 ) {
# creat a problem

  lambda <- -.20
  Beta <- .2

L <- 25 
W <- as.matrix(bandSparse(L,L,k=list(-1*sqrt(L),-1,1,sqrt(L))))*1

# create a sequare root matrix for our simulation data with lambda = 0.20
sigma <- solve(diag(L) - lambda * W)
sigma.eigen <- eigen(sigma)
sigma.half <- sigma.eigen$vectors %*% diag(sqrt(sigma.eigen$values)) 

aa <-c()
bb <-c()
cc <- c()

for( i in 1:50000 ) {

  # setting up Z with appropriate variance and mean, pnorm(.5) = 0.6914
  Z <- sigma.half %*% rnorm(L) + Beta 
  Y <- matrix(as.numeric(Z > 0),ncol=1)
  
  X <- matrix(1,nrow=L,ncol=1)
  
  # set initial beta
  Beta.init <- matrix(0.2,nrow=1)
  lambda.init <- -.20
  
  Beta0 <- 0 
  Sigma0 <- 20 
  
  iter <- 1
  
  a <- probitGibbsSpatial(Y,X,W,L,Beta.init,lambda.init,Beta0,Sigma0,iter)
  aa[i] <- mean(a$Beta)
  bb[i] <- mean(a$lambda)
  cc[i] <- mean(Z) 

}

}


if ( test2 ) {

  lambda.corn <- .20
  Beta.corn <- -.5
  lambda.soy <- .20
  Beta.soy <- -.4

  # probability that a field initially starts off as corn
  corn.prob <- .80

  # creat a problem
  L <- 100 
  W <- as.matrix(bandSparse(L,L,k=list(-1*sqrt(L),-1,1,sqrt(L))))*1

  # identify parcels with corn 
  is.corn <- matrix(0,nrow=100,ncol=1)
  is.corn[sample(1:100,p=corn.prob),] <- 1

  # create a sequare root matrix for our simulation data with lambda = 0.20
  sigma.corn <- solve(diag(L) - .20 * W)
  sigma.eigen.corn <- eigen(sigma)
  sigma.half.corn <- sigma.eigen$vectors %*% diag(sqrt(sigma.eigen$values)) 
 
  Beta.corn <-c()
  lambda.corn <-c()
  Zobs.corn <- c()
  
  sigma.soy <- solve(diag(L) - .20 * W)
  sigma.eigen.soy <- eigen(sigma)
  sigma.half.soy <- sigma.eigen$vectors %*% diag(sqrt(sigma.eigen$values)) 
  
  Beta.soy <-c()
  lambda.soy <-c()
  Zobs.soy <- c()
  
  for( i in 1:100 ) {
  
    # setting up Z with appropriate variance and mean, pnorm(.5) = 0.6914
    Z.corn <- sigma.half.corn %*% rnorm(L) + Beta.corn 
    Y.corn <- matrix(as.numeric(Z.corn > 0),ncol=1)
    X.corn <- matrix(1,nrow=L,ncol=1)
   
    # setting up Z with appropriate variance and mean, pnorm(.5) = 0.6914
    Z.soy <- sigma.half.soy %*% rnorm(L) + Beta.soy 
    Y.soy <- matrix(as.numeric(Z.soy > 0),ncol=1)
    X.soy <- matrix(1,nrow=L,ncol=1)
   
    # set initial conditions for observations
    Y.corn[is.corn != 1,] <- -1 
    Y.soy[is.corn == 1,] <- -1 

    # set initial beta
    Beta.init <- matrix(0.1,nrow=1)
    lambda.init <- 0 
    
    Beta0 <- 0.5 
    Sigma0 <- 10
    
    iter <- 1
    
    sim.corn <- probitGibbsSpatial(Y.corn,X.corn,W,L.corn,Beta.init,lambda.init,Beta0,Sigma0,iter)
    sim.soy <- probitGibbsSpatial(Y.soy,X.soy,W,L.soy,Beta.init,lambda.init,Beta0,Sigma0,iter)
  
    Beta.corn[i] <- mean(a$Beta)
    lambda.corn[i] <- mean(a$lambda)
    zobs.corn[i] <- mean(Z) 
    
    Beta.soy[i] <- mean(a$Beta)
    lambda.soy[i] <- mean(a$lambda)
    zobs.soy[i] <- mean(Z) 
  }

}

