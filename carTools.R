########################################## carTools.R  ###############################################
# A set of tools to simulate CAR models


########################################## LIBRARY ###############################################
library(msm)  # for rtnorm
library(mvtnorm)
library(Matrix)





########################################## INCLUDES ###############################################
source('mh.R')
source('simCrop.R')



########################################## FUNCTION ###############################################


## a function to check if the value for rho ensures that I-rho *W is positive definite
# if a value of rho is not provided the allowable range to ensure positive definiteness is 
# provided
carTools.checkRho <- function(W, rho) {

  if( !missing(rho) ) {
    # for matrix valued rho
    if( !is.null(nrow(rho)) ) {
       if(sum(eigen(diag(nrow(W)) - rho %*% W)$values <= 0) != 0) {
        return(F) 
       } else {
         return(T)
       }
     } else {
    # for non-matrix valued rho
       if(sum(eigen(diag(nrow(W)) - rho * W)$values <= 0) != 0) {
        return(F) 
       } else {
         return(T)
       }
     }
  }
    
  # if rho is missing we return the scalar range of rho
  W.eigen.values <- eigen(W)$values

  return( sort( 1/range(W.eigen.values) ) )
}


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


## probit Gibbs function ## 
carTools.probitGibbs <- function(y,X,Beta,Beta0,Sigma0,iter) {

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


## probit Gibbs function ## 
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


carTools.deviates <- function( rho, W, X, Beta) {
  n <- nrow(X)
  if( is.null(n) ) n <- length(X)

  if( !is.null(nrow(rho)) ) { 
    L <- chol(diag(n) - rho %*% W) 
  } else {
    L <- chol(diag(n) - rho * W ) 
  }

  if(is.null(nrow(Beta))) Beta <- matrix(Beta,ncol=1)
  if(is.null(nrow(X)))    X <- matrix(X,ncol=1)

  return(  L %*% ( X %*% Beta + t(L) %*% L %*% rnorm(n) )  )
}


carTools.generateCropTypes <- function(a, p, rho, X, Beta) {

  if( !missing(p) ) {
    return( simCrop.generateCropTypes(a.neighbors,p) )
  }

  if( is.null(a$neighbors) ) {
    # figure out which QQS are neighbors 
    print("No neighbors exist, assuming rook (shared edge)")
    a <- simCrop.getNeighbors(a)
  }

  # create the distance matrix for a
  W <- simCrop.createRookDist(a)
  
  # W is sorted by object, so we need to sort our input by object
  myObjects <- a$map[,'object']
  myObjects.sortIndex <- sort( myObjects, index.return=T)
  myObjects.sort <-myObjects.sortIndex$x
  myObjects.sortIndex <-myObjects.sortIndex$ix
  # this gives us a way to un-sort the result
  myObjects.unsortIndex <- sort(myObjects.sortIndex, index.return=T)$ix

  # handle missing X and Beta
  n <- length(myObjects)
  if( missing(X) ) {
    X <- matrix(1,nrow=n,ncol=1)
    Beta <- list()
    for(i in 1:length(rho) ) Beta[[i]] <- rep(0,times=length(rho)) 
  }

  X.sort <- X[ myObjects.sortIndex,]
  Y <- matrix(1,nrow=n,ncol=1)

  for( i in 1:length(rho) ) {

    priorState <- a$cropType[,'x'] == i 

    if( sum(priorState) != 0 )  {
      Y <- carTools.deviates( rho[[i]], W, X.sort, Beta[[i]] )
      Y[priorState] <- (Y[myObjects.unsortIndex])[priorState] 
    }
   
  }
  # now we need to 
  a$cropType <- cbind( a$cropType, 1 + (Y>0) )

  return(a)
}


