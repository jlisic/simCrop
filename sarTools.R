########################################## carTools.R  ###############################################
# A set of tools to simulate CAR models


########################################## LIBRARY ###############################################
library(msm)  # for rtnorm
library(mvtnorm)
library(Matrix)
library(mgcv)




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


sarTools.probitGibbsSpatial <- function( a, Beta.init, lambda.init, beta0,Sigma0, iter, m ) {

  myObjects <- a$cropType[,'myObjects']
  myObjects.sort <- sort( myObjects, index.return=T)$ix
  priorYears <- ncol( a$cropType) - 2  

  # take care of Y
  Y <- matrix(a$cropType[,c(-1,-2)],ncol=priorYears)
  
  Y.sort <- matrix(1:length(Y),ncol=priorYears)
  Y.sort <- c(Y.sort[myObjects.sort,])
  Y <- Y[Y.sort]
  Y <- matrix(Y,ncol=1)

  # take care of X
  X <- matrix(sarTools.priorStateDesignMatrix(a),ncol=priorYears)
  X <- X[ Y.sort, ]

  # take care of W  
  W <- simCrop.createRookDist(a)
  fieldSize <- nrow(W)

  X <<- X
  Y <<- Y

  result <- sarTools.probitGibbsSpatialRun(
                                      Y,
                                      a$crops,
                                      X,
                                      W,
                                      Beta.init,
                                      lambda.init,
                                      Beta0,
                                      Sigma0,
                                      iter,
                                      m) 
  return( result)
}

sarTools.probitGibbsSpatial2 <- function( a, Beta.init, rho.init, beta0,Sigma0, iter, m,thinning,burnin ) {

  require('spatialprobit')

  myObjects <- a$cropType[,'myObjects']
  myObjects.sort <- sort( myObjects, index.return=T)$ix
  priorYears <- ncol( a$cropType) - 2  

  # take care of Y
  Y <- matrix(a$cropType[,c(-1,-2)],ncol=priorYears)
  
  Y.sort <- matrix(1:length(Y),ncol=priorYears)
  Y.sort <- c(Y.sort[myObjects.sort,])
  Y <- Y[Y.sort]
  Y <- matrix(abs(Y-2),ncol=1) 
  Y <<- Y

  # take care of X
  X <- sarTools.priorStateDesignMatrix(a)
  X <- X[ Y.sort, ]
  X <<- X

  # take care of W  
  W <- simCrop.createRookDist(a)
  W.big <- kronecker(diag(priorYears),W) 
  fieldSize <- nrow(W)


  result <- sar_probit_mcmc(  
                                      y=Y,
                                      X=X,
                                      W=as(W.big,"CsparseMatrix"),
                                      ndraw=iter,
                                      burn.in=burnin,
                                      thinning=thinning,
                                      # prior - Not Set, using diffuse
                                      start= list( rho=rho.init, beta= Beta.init),
                                      m=m,
                                      showProgress=TRUE
                                     ) 
                                      
  return( result)
}


sarTools.gibbsSpatialRun <- function(Y,X,W,Beta.init,rho.init,Beta0,Sigma0,iter,m) {

  # set initial conditions
  Beta <- Beta.init
  rho <- rho.init 
  #init some values
  n <- nrow(Y)     # number of observations
  K <- n / nrow(W) # number of years
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

    # Lambda update
    Lambda.inv <- diag(n) - rho * W.big 
    Sigma.inv <- t(Lambda.inv) %*% Lambda.inv
 
    # because we don't have latents 
    Z <- Y 

    if( F ) { 

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
      rho.save[i] <- mh.lambda.sar(Z,W, X%*% Beta,0,1,100,rho.range )
      rho <- rho.save[i] 
    }
  }
  
  return( list( Beta = Beta.save, rho = rho.save) )
}


## probit Gibbs function ## 
# Y vector of categorical responses in row major form, repeating for each year, length = (number of years) x fieldSize
# X matrix of covariates  in row major form, repeating for each year, length = (number of years) x fieldSize
# W matrix in row major form of spatial neighborhoods, dim is fieldSize x fieldSize
# fieldSize, number of observations in a given year
#  
sarTools.probitGibbsSpatialRun <- function(Y,states,X,W,Beta.init,rho.init,Beta0,Sigma0,iter,m) {

  Beta0 <- matrix(Beta0,ncol=1)
  Sigma0 <- matrix(Sigma0,ncol=1)
  Beta.n <- nrow(Beta0)
  
  # set initial conditions
  Beta <- Beta.init
  rho <- rho.init 
  #init some values
  n <- nrow(Y)     # number of observations
  K <- n / nrow(W)
  p <- ncol(X)     # number of covariates

  Beta.save <- matrix(0,nrow=iter,ncol=length(Beta)) 
  rho.save <- matrix(0,nrow=iter,ncol=1) 

  W.big <- kronecker(diag(K),W) 
  rho.range <- sort( 1/range(eigen(W)$values) )

  Z <- matrix(0,nrow=n,ncol=1) 
  trunc.point <- Z

  # inverse of the prior variance
  T.inv <- diag( c(1/Sigma0^2) )
  T.inv <<- T.inv
  B.star.inv <- solve( t(X) %*% X + T.inv ) 
  B.star.inv.sqrt <- mroot(B.star.inv) # require mgcv

  # the mcmc loop
  for(i in 1:iter) {
    print(i)
    last.time <- proc.time()
    # Lambda update
    Lambda <- diag(n) - rho * W.big 
    Sigma.inv <- t(Lambda) %*% Lambda
  
    # generate deviates for the latent variables
    Z <- matrix(0,nrow=n,ncol=1) 
    for( k in 1:m) {
      for( j in 1:n) {

        Z.sd <- 1/sqrt(Sigma.inv[j,j])
        
        trunc.point[j] <-  (X[j,] %*% Beta - (Sigma.inv[j,-j]/Sigma.inv[j,j]) %*% Z[-j] )/ ( -1 * Z.sd ) 

        if( Y[j] == states[1]) {
          Z[j] <- rtnorm( 1, upper=trunc.point[j], sd=1 ) 
        } else if( Y[j] == states[2] ) {
          Z[j] <- rtnorm( 1, lower=trunc.point[j], sd=1 ) 
        } else {
          print( "Error:  Unknown State" )
        }
      }
    
      Z <- X %*% Beta + Z
  
      if( T ) { 
        B <- B.star.inv %*% (t(X) %*% Lambda %*% Z  + T.inv %*% Beta0)  
        
        # generate deviates for beta/mu
        Beta.save[i,] <- B + B.star.inv.sqrt %*% rnorm(Beta.n) 
  
        Beta <- matrix(Beta.save[i,],ncol=1)
        print(Beta)
      } else {
        Beta.save <- t(Beta)
      }

    }
    # generate lambda deviate
    if( T ) {
      rho.save[i,] <- mh.lambda.sar(Z,W, X%*% Beta,0,1,100,rho.range )
      rho <- rho.save[i,] 
    } else {
      rho.save <- t(rho)
    }
    print(proc.time() - last.time)
  }

  return( list( Beta = Beta.save, rho = rho.save) )
}


sarTools.deviates <- function( rho, W, X, Beta) {
  n <- nrow(X)
  if( is.null(n) ) n <- length(X)
  Lambda <- diag(n) - rho * W
  Lambda.inv <- solve(Lambda)
  print( sprintf( "Simulating with Beta=%f Rho=%f",Beta,rho) )
  if(is.null(nrow(Beta))) Beta <- matrix(Beta,ncol=1)
  if(is.null(nrow(X)))    X <- matrix(X,ncol=1)

  rhoRange <-carTools.checkRho(W) 
  epsilon <- rnorm(n)
  Y <- Lambda.inv %*% (X %*% Beta + epsilon) 

  Beta.hat <- solve( t(X) %*% X )  %*% t(X) %*% Lambda %*% Y 
  print( sprintf( "MLE Beta=%f Rho=%f",Beta.hat, 
    optim( rho, sarCheck, Y=Y, W=W, X=X, Beta=Beta, lower=rhoRange[1] + 0.0001, upper=rhoRange[2] - 0.0001,method="Brent" )$par ))

  return( Y ) 
}





sarTools.deviates.simple <- function( rho, W ) {
  n <- nrow(W)
  Lambda <- diag(n) - rho * W
  Lambda.inv <- solve(Lambda)
  print( sprintf( "Simulating with Beta=%f Rho=%f",0,rho) )
  return( Lambda.inv %*% matrix(rnorm(n),ncol=1) ) 
}


sarTools.priorStateDesignMatrix <- function(a,priorYear) {
  # get all the years of data we want
  # note the first column is for object id 'myObjects'
  # the last row is also not used
  X <- a$cropType
 
  if( missing(priorYear) ) { 
    # number of years
    n <- ncol(X) 
    # create a column major form listing of prior years
    py <- c(X[, c(-1,-n) ] )
  } else {
    # create a column major form listing of prior years
    py <- c(X[, priorYear+1 ] )
  }

  # create design matrix
  y <- rep(0,times=length(a$crops))
  return( t(sapply(py,function(x) {y[x] = 1; return(y)}) ) )
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

  # get the number of years
  years <- ncol(a$cropType) - 1 

  # W is sorted by object, so we need to sort our input by object
  myObjects <- a$map[,'object']
  myObjects.sortIndex <- sort( myObjects, index.return=T)
  myObjects.sort <-myObjects.sortIndex$x
  myObjects.sortIndex <-myObjects.sortIndex$ix

  # this gives us a way to un-sort the result
  myObjects.unsortIndex <- sort(myObjects.sortIndex, index.return=T)$ix

  priorState <- sarTools.priorStateDesignMatrix(a,priorYear=years)


  if( missing(X) ) {
    X <- priorState
  } else {
    X <- cbind(X,priorState)
  }
  X.sort <- X[ myObjects.sortIndex,]
  Y <- sarTools.deviates( rho, W, X.sort, Beta )
  Y <- Y[ myObjects.unsortIndex,]
   
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


sarCheck.simple <- function( rho, Y, W ) {
  n <- nrow(W)
  Lambda <- diag(n) - rho * W
  Z <- Lambda %*% Y   
  return( -1 * (log(det(Lambda))    -1/2 * t(Z) %*% Z )) 
}




