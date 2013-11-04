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


sarTools.probitGibbsSpatial <- function( a, Beta.init, lambda.init, beta0,Sigma0, iter ) {


  myObjects <- a$cropType[,'myObjects']
  myObjects.sort <- sort( myObjects, index.return=T)$ix
  
  # take care of Y
  Y <- a$cropType[,-1]
  Y <- Y[myObjects.sort,]
  
  # take care of X
  X <- matrix(1,nrow=nrow(Y)*(ncol(Y) - 1),ncol=1)
    
  W <- simCrop.createRookDist(a)
  
  fieldSize <- nrow(W)

  result <- list() 
  for( i in 1:length(Beta.init) ) { 

    Y.change <- Y
    Y.change[,-1][!(Y.change==i)[,-ncol(Y.change)]] <- 0 
    Y.change <- Y.change[,-1] - 1

    Y.change <- matrix(Y.change,ncol=1)

    result[[i]] <- sarTools.probitGibbsSpatialRun(
                                      Y.change,
                                      X,
                                      W,
                                      Beta.init[[i]],
                                      lambda.init[[i]],
                                      Beta0[[i]],
                                      Sigma0[[i]],
                                      iter,
                                      10) 
    print( sprintf("Finished Running %d",i) )
  } 
  return( result)
}


## probit Gibbs function ## 
# Y vector of categorical responses in row major form, repeating for each year, length = (number of years) x fieldSize
# X matrix of covariates  in row major form, repeating for each year, length = (number of years) x fieldSize
# W matrix in row major form of spatial neighborhoods, dim is fieldSize x fieldSize
# fieldSize, number of observations in a given year
#  
sarTools.probitGibbsSpatialRun <- function(Y,X,W,Beta.init,rho.init,Beta0,Sigma0,iter,m) {

  # set initial conditions
  Beta <- Beta.init
  rho <- rho.init 
  #init some values
  n <- nrow(Y)     # number of observations
  K <- n / nrow(W)
  p <- ncol(X)     # number of covariates

  Beta.save <- c()
  rho.save <- c()

  W.big <- kronecker(diag(K),W) 
  rho.range <- sort( 1/range(eigen(W)$values) )

  Z <- matrix(0,nrow=n,ncol=1) 
  trunc.point <- Z

  # inverse of the prior variance
  T.inv <- diag(ncol(X)) / Sigma0^2
  B.star.inv <- solve( t(X) %*% X + T.inv ) 

  # the mcmc loop
  for(i in 1:iter) {
    print( sprintf("Lambda Update %d ", i) )
    last.time <- proc.time()

    # Lambda update
    Lambda.inv <- diag(n) - rho * W.big 
    Sigma.inv <- t(Lambda.inv) %*% Lambda.inv
  
    # the usual X'X
    
    print( proc.time() - last.time)  
    print( sprintf("Z Generation %d ", i) )
    last.time <- proc.time()
     

    # generate deviates for the latent variables
    for( k in 1:m) {
      for( j in 1:n) {
        trunc.point[j] <- -1/sqrt(Sigma.inv[j,j]) * ( X[j] %*% Beta + (Sigma.inv[j,-j]/Sigma.inv[j,j]) %*% Z[-j]) 
  
        if( Y[j] == 1) {
          Z[j] <- rtnorm( 1, lower=trunc.point[j], sd=1 ) 
        } else if( Y[j] == 0) {
          Z[j] <- rtnorm( 1, upper=trunc.point[j], sd=1 ) 
        } else {
          Z[j] <- rnorm( 1 ) 
        }
        
      }
    }
    print( proc.time() - last.time)  
    opt.result <- optim( 0, sarCheck, X=X, Y=Z, Beta=Beta, W=W, lower=rho.range[1] + 0.0001, upper=rho.range[2] - 0.0001,method="Brent" )
    print(opt.result)
    if( F ) { 
      print( sprintf("Beta Generation %d ", i) )
      last.time <- proc.time()

      B <- B.star.inv %*% (t(X) %*% Lambda.inv %*% Z  + T.inv %*% Beta0)  
      
      # generate deviates for beta/mu
      Beta.save[i] <- rnorm(1,B, B.star.inv )
      Beta <- Beta.save[i]
      print( proc.time() - last.time)  
    } else {
      Beta.save[1] <- Beta
    }

    # generate lambda deviate
    if( T ) {
      print( sprintf("Rho Generation %d ", i) )
      rho.save[i] <- mh.lambda.sar(Z,W, X%*% Beta,0,1,100,rho.range )
      rho <- rho.save[i] 
      print( proc.time() - last.time)  
    }
  }
  
  return( list( Beta = Beta.save, rho = rho.save) )
}


sarTools.deviates <- function( rho, W, X, Beta) {
  n <- nrow(X)
  if( is.null(n) ) n <- length(X)


  if( !is.null(nrow(rho)) ) { 
    Lambda <- (diag(n) - rho %*% W)
  } else {
    Lambda <- (diag(n) - rho * W)
  }
 
  Lambda.inv <- solve(Lambda)

  print( sprintf( "%d - %f",i,rho) )

  if(is.null(nrow(Beta))) Beta <- matrix(Beta,ncol=1)
  if(is.null(nrow(X)))    X <- matrix(X,ncol=1)

  return( solve(Lambda) %*% X %*% Beta + Lambda %*% rnorm(n) ) 
}


sarTools.generateCropTypes <- function(a, p, rho, X, Beta) {

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
      Y <- sarTools.deviates( rho[[i]], W, X.sort, Beta[[i]] )
      Y[priorState] <- (Y[myObjects.unsortIndex])[priorState] 
    }
   
  }
  # now we need to 
  a$cropType <- cbind( a$cropType, 1 + (Y>0) )

  return(a)
}



sarCheck <- function( rho, Y, X, Beta, W ) {

  n <- nrow(W)

  Lambda <- diag(n) - rho * W
  Z <- Lambda %*% Y - X %*% Beta  

  return( -1 * (log(det(Lambda))    -1/2 * t(Z) %*% Z )) 

}





