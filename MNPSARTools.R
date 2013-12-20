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
source('minDim.R')

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
# Y vector of categorical responses in row major form, repeating for each year, length = (number of years) x fieldSize
# X matrix of covariates  in row major form, repeating for each year, length = (number of years) x fieldSize
# W matrix in row major form of spatial neighborhoods, dim is fieldSize x fieldSize
sarTools.probitGibbsSpatial <- function( a, fun, ... ) {

  # number of alternatives
  J <- length(a$crops) - 1
  # number of years
  priorYears <- ncol( a$cropType) - 2  

  # take care of Y
  Y <- matrix( sortCrop(a$cropType[,c(-1,-2)],a), ncol=1) 
  Y <<- Y

  # take care of X
  X <- sarTools.priorStateDesignMatrix(a,sortResult=T)
  X <<- X
  X <- kronecker( X, diag(J) )
  X <<- X

  # take care of W  
  W <- simCrop.createRookDist(a)

  # run the function
  result <- fun( Y=Y, states=a$crops, X=X, W=W, ... )

  return( result)
}




# a sensible conditional SAR model
sarTools.probitGibbsSpatialRunConditional <- function(
# parameters (these are provided by sarTools.probitGibbsSpatial)
  Y,
  states,
  X,
  W,
#init parameters
  beta.init,
  rho.init,
  Z.init,
  alpha.init,
#hyper parameters
  Beta0,
  Sigma0,
  Wishart0.Z,
  Wishart0.alpha,
#runtime parameters
  iter,
  m,
  thinning,
  burnIn=0
  ) {

# required packages
  require(tmvtnorm)

  Beta0 <- matrix(Beta0,ncol=1)
  Sigma0 <- matrix(Sigma0,ncol=1)
  
  Beta.n <- nrow(Beta0)
  
  # set initial conditions
  Beta <- beta.init
  rho1 <- rho.init[1]
  rho2 <- rho.init[2] 
  Z <- Z.init
  V <- alpha.init

  #init some values
  NJ <- nrow(X)    # number of observations (times J*K )
  N <- length(Y)     # number of observations * K
  n <- nrow(W)     # number of observations within a year 
  K <- N / n       # number of years
  p <- ncol(X)     # number of covariates
  J <- NJ / N 
  nJ <- n * J 

  # take care of Y
  Y.lower <- rep(0,times=NJ) 
  Y.constraints <- bdiag(lapply( Y , function(x) DMat(J,x)))

  Y.constraints <<- Y.constraints

  # things to save
  Beta.save <- matrix(0,nrow=iter,ncol=length(Beta))  # Beta values
  rho.save <- matrix(0,nrow=iter,ncol=2)              # rho values
  Sigma1.inv.save <- matrix(0,nrow=iter,ncol=J^2)                # proportion parameter for variances
  Sigma2.inv.save <- matrix(0,nrow=iter,ncol=J^2)                # proportion parameter for variances
  V.save <- c()
  Z.save <- c()

  rho.range <- sort( 1/range(eigen(W)$values) )

  # create U matrix  1_K \otimes (diag(n) \otimes 1_J)
  U <- kronecker(matrix(1,ncol=1,nrow=K),diag(n*J))

  # static mu (zero) value for rho2
  mu2 <- matrix(0,ncol=1,nrow=n*J)

  # inverse of the prior variance
  S.inv <- diag( c(1/Sigma0^2) )

  # Create initial value for Sigma2.inv  
  Wishart0.Z.inv <- solve(Wishart0.Z[[2]])
  Wishart0.alpha.inv <- solve(Wishart0.alpha[[2]])
      
  # update rho step 
  Lambda1 <- kronecker( diag(n) - rho1 * W, diag(J))
  Lambda2 <- kronecker( diag(n) - rho2 * W, diag(J)) 
  Lambda1.K <- kronecker(diag(K), Lambda1)

  last.time <- proc.time()

#  Sigma1.inv <- diag(2)  
#  Sigma2.inv <- diag(2)  

  print(dim(Z))
  print(dim(X))
  print(dim(Lambda1.K))

  # the mcmc loop
  for(i in 1:(iter+burnIn) ) {

    if(i %% 10 == 0) {
      print(proc.time() - last.time)
      print(i)
      last.time <- proc.time()
    }

    for(l in 1:thinning) {
      # because occasionally we get location parameters that don't work out.
      retryCount <- 0
      while( TRUE ) {
      
        ## 3. generate Sigma2.inv deviate
        V.alt <- matrix(Lambda2 %*% V, nrow=J)
        Sigma2.inv <- rWishart(1 , n + Wishart0.alpha[[1]], solve( V.alt %*% t(V.alt) + Wishart0.alpha.inv ))[,,1] 
  
        ##tmp
        mu <- X %*% Beta + U%*%V
        Z.alt <- matrix(Lambda1.K %*%Z - mu, nrow=J)
        Sigma1.inv <- rWishart.control(N + Wishart0.Z[[1]], solve( Wishart0.Z.inv + Z.alt %*% t(Z.alt) ))  

        Sigma1.inv <<- Sigma1.inv

        # update the variance 
        S1.inv <- kronecker( diag(n), Sigma1.inv )
        S2.inv <- kronecker( diag(n), Sigma2.inv )
        S1.K.inv <- kronecker( diag(K), S1.inv ) 
        XX <- t(X) %*% S1.K.inv %*% X 
        UU <- t(U) %*% S1.K.inv %*%  U 
  
        ## 1. generate rho1 deviate
        rho1 <- mh.lambda.mnpsar(Z=Z,W=W,mu=mu,Sigma=Sigma1.inv,x0=rho1,iter=1,burnIn=100,rho.range=rho.range)
  
        ## 2. generate rho2 deviate
        rho2 <- mh.lambda.mnpsar(Z=V,W=W,mu=mu2,Sigma=Sigma2.inv,x0=rho2,iter=1,burnIn=100,rho.range=rho.range)
     
        # update rho step 
        Lambda1 <- kronecker( diag(n) - rho1 * W, diag(J))
        Lambda2 <- kronecker( diag(n) - rho2 * W, diag(J)) 
        
        Lambda1.inv <- solve(Lambda1)
  
        LH1L <- Lambda1 %*% S1.inv %*% Lambda1 
        LH2L <- Lambda2 %*% S2.inv %*% Lambda2
        
        Lambda1.K <- kronecker(diag(K), Lambda1)
        Lambda1.K.inv <- kronecker(diag(K), Lambda1.inv) 
        Sigma1.K.inv <- kronecker(diag(K), LH1L) 
  
        Sigma1.K.inv <<- Sigma1.K.inv 
        
        ############################################################################
        ## 5. generate deviates for the random effect latent variables
        # L = Lambda1, O - Sigma1.inv
        # M = Lambda2, A - Sigma2.inv
        # t(z - L.inv xb - L.inv ua)LOL(z - L.inv xb - L.inv ua) + aMAMa
        # V.var =  solve(  t(u)Ou + MAM   )
        # V.mu = V.var.inv t(u) O ( L z - xb )  
        ############################################################################
        
        Sigma2.cond <- solve(UU + LH2L )
        muV <- Sigma2.cond %*% t(U) %*% S1.K.inv %*% (Lambda1.K %*% Z - X%*%Beta)
        V <- matrix( rmvnorm(1, mean= muV, sigma=Sigma2.cond),  ncol=1)
  
        
  
  
        ############################################################################
        ## 4. generate new Beta
        # L = Lambda1, O - Sigma1.inv
        # B = Beta0
        # t(z - L.inv xb - L.inv ua)LOL(z - L.inv xb - L.inv ua) + t(B - b)S.inv(b-B) 
        # t(x) O x  + S.inv 
        # S.inv B + t(x) L.inv LOL( z - L.inv ua) = S.inv B + t(x) O ( L z - ua)
        ############################################################################
  
        
  
        # variance of the posterior distribution of Beta
        Beta.post.var <- solve( S.inv + XX ) 
       
        # mean of the posterior distribution of Beta
        Beta.post.mean <- Beta.post.var %*% ( t(X) %*% S1.K.inv %*% (Lambda1.K %*% Z - U%*%V)  + S.inv %*% Beta0)

        Beta <- matrix(rmvnorm(n=1,mean=Beta.post.mean,sigma=Beta.post.var),ncol=1)

        ############################################################################
        ## 6. generate deviates for the truncated latent variables

        Z.new <- c()
        for(k in 1:K) {
          interval.k <- (nJ*(k-1) + 1):(nJ*k)
          X.k <- X[interval.k,]
          U.k <- U[interval.k,]
          muZ <- Lambda1.inv %*% (X.k %*% Beta + U.k %*% V)
          Y.lower.k <- Y.lower[interval.k]
          Y.constraints.k <- Y.constraints[interval.k,interval.k]

          Z.new <- rbind(
              Z.new,
              matrix( rtmvnorm( n=1, mean=c(muZ),H=LH1L, lower=Y.lower.k, D=as.matrix(Y.constraints.k), algorithm="gibbs",burn.in.samples=m, ), ncol=1)
             )

        }

        Z <- Z.new

        if( !is.nan(sum(Z))) {
          break
        } else {
          if( i == 1 ) {
            Z <- Z.init
          } else {
            Z <- Z.save[,i-1]
          }
        }

        retryCount <- retryCount + 1
      }
      if(retryCount > 1) print( sprintf("retryCount = %d", retryCount))



    } # finish thinning 
    Beta.save[i - burnIn,] <- Beta  # save our result
    rho.save[i -  burnIn,1] <- rho1
    rho.save[i -  burnIn,2] <- rho2 
    Sigma1.inv.save[i -  burnIn,] <- Sigma1.inv 
    Sigma2.inv.save[i -  burnIn,] <- Sigma2.inv 
    V.save <- cbind(V.save,V)
    Z.save <- cbind(Z.save,Z)
  } # finish iteration

  return( list( Beta = Beta.save, Rho = rho.save,  Sigma1.inv = Sigma1.inv.save, Sigma2.inv = Sigma2.inv.save, V = V.save, Z = Z.save) )
}


  

# The form of X is 
#
# X_i,j
#
# i - choice
# j - position
#
# e.g.
# 111
# 211
# 121
# 221
#
sarTools.deviates <- function( rho, W, X, Beta, mu2, Sigma ) {
  n <- nrow(W)
  Lambda <- (diag(n) - rho * W)
  Lambda.inv <- solve(Lambda)

  print( sprintf( "Simulating with Beta:" ))
  print(Beta)
  print( sprintf("Rho=%f, Sigma:",rho) )
  print(Sigma)

  if(minDim(Beta) == 1) Beta <- matrix(Beta,ncol=1)
  if(minDim(X) == 1)    X <- matrix(X,ncol=1)
  if(minDim(Sigma) == 1) Sigma <- matrix(Sigma,ncol=1)
  if(missing(mu2)) mu2 <- matrix(0,nrow=nrow(X),ncol=1) 

  rhoRange <-carTools.checkRho(W)
  J <- minDim(Sigma)

  # if we only have tau we treat each observation individually 
  if( J == 1 ) { 
    epsilon <- rnorm(n)
    Z <- Lambda.inv %*% (X %*% Beta + mu2 + epsilon / Sigma ) 

  # if we do have sigma then we repeat the location parameter multiple times
  } else {
    # get the dim of Sigma (so we now how many choices) 
    Z <-  matrix( c( t( Lambda.inv %*% rmvnorm(n,sigma=Sigma) )),ncol=1)
    
    Z <-  kronecker(Lambda.inv, diag(J)) %*% (X %*% Beta + mu2  ) + Z 
  }

  # sanity checks
#  FullSigmaRoot <-  bdiag( rep( list( solve(Sigma) ), times=n) ) 
#  Beta.hat <- solve( t(X) %*%FullSigmaRoot %*% X )  %*% t(X) %*% FullSigmaRoot %*% kronecker(Lambda,diag(J)) %*% ( Z - kronecker(Lambda.inv, diag(J)) %*% mu2 )  
#  print( sprintf( "LSE Beta="))
#  print(Beta.hat)
#  print( sprintf( "MLE Rho="))
#
#  print(optim( rho, sarCheck, Y=Z, W=W, X=X, Beta=Beta, H=FullSigmaRoot, lower=rhoRange[1] + 0.0001, upper=rhoRange[2] - 0.0001,method="Brent" )$par )

  return( Z ) 
}








sarTools.priorStateDesignMatrix <- function(a,priorYear,sortResult=F) {
  # get all the years of data we want
  # note the first column is for object id 'myObjects'
  # the last row is also not used
  X <- a$cropType
  
  if( sortResult ) {
    X <- sortCrop( X, a)
  } 

 
  if( missing(priorYear) ) { 
    # number of years
    p <- ncol(X) 
    # create a column major form listing of prior years
    py <- c(X[, c(-1,-p) ] )
  } else {
    # create a column major form listing of prior years
    py <- c(X[, priorYear+1 ] )
  }

  # create design matrix
  y <- rep(0,times=length(a$crops))
  return( t(sapply(py,function(x) {y[x] = 1; return(y)}) ) )
}


sarTools.generateCropTypes <- function(a, p, rho, X, Beta, Sigma.list) {

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

  priorState <- sarTools.priorStateDesignMatrix(a,priorYear=years)

  J <- minDim(Sigma.list[[1]])

  # take care of X
  if( missing(X) ) {
    X <- priorState
  } else {
    X <- cbind(X,priorState)
  }
  X.sort <- sortCrop(X,a)
  X.sort <- kronecker(X.sort, diag(J) )

  X.sort <<- X.sort

  # take care of global error
  if( length(Sigma.list) == 2 ) {

    globalError <- sarTools.deviates( rho=rho[2], W=W,X=matrix(1,ncol=1,nrow=nrow(W)*J),Beta=0, Sigma=Sigma.list[[2]])
    V1 <<- globalError

    a$globalError <- cbind( rep( a$neighbors[,1], each=J), unsortCrop(globalError,a) )

    colnames(a$globalError) <- c('myObjects','error')
  }


  # if there is a global error available add it
  if(!is.null(a$globalError) ) {
    Y <- sarTools.deviates( 
                           rho[1], 
                           W, 
                           X.sort, 
                           Beta, 
                           mu2= sortCrop(a$globalError[,'error'], a), 
                           Sigma=Sigma.list[[1]] )
  } else {
    Y <- sarTools.deviates( rho[1], W, X.sort, Beta, Sigma=Sigma.list[[1]] )
  }
  Y1 <<- Y 
  
  Y <- unsortCrop( Y, a)
  # now we need to apply our rule to determine category
  a$cropType <- cbind( a$cropType, applyGroup( Y, J  ) )
  a$cropValue <- cbind( a$cropValue, Y )

  return(a)
}



sortCrop <- function( Z, a) {

  myObjects <- a$cropType[,'myObjects']
  object.sort <- sort(myObjects,index.return=T)$ix
  
  cols <- ncol(Z)
  if( is.null(cols)) cols <- 1

  rows <- length(Z) / cols 
  n <- length(object.sort)

  Z <- matrix( t(Z), ncol=n )[,object.sort] # it's sorted here 

  # now we need to get it back to the original shape
  Z <- t( matrix( Z, nrow=cols ))

  return( Z )
}

unsortCrop <- function( Z, a) {

  myObjects <- a$cropType[,'myObjects']
  object.sort <<- sort(myObjects,index.return=T)$ix
  object.unsort <<- sort(object.sort, index.return=T)$ix

  1:length(myObjects) 

  cols <- ncol(Z)
  if( is.null(cols)) cols <- 1
  
  rows <- length(Z) / cols 
  n <- length(object.sort)

  Z <- matrix( t(Z), ncol=n )[,object.unsort] # it's sorted here 

  # now we need to get it back to the original shape
  Z <- t( matrix( Z, nrow=cols ))

  return( Z )
}




DMat <- function( J, j ) {

  if( J <= 3 ) {

    DMatrix <-  list(
        matrix( 1,nrow=1),

        matrix( c(
                  -1,  0,
                   0, -1, 
                  
                   1,  0,
                   1, -1,

                  -1,  1,  
                   0,  1), ncol=2, byrow=T),

        matrix( c( 
                  -1,  0,  0,
                   0, -1,  0,
                   0,  0, -1,
                  
                   1, -1,  0,   
                   1,  0, -1,   
                   1,  0,  0,

                  -1,  1,  0,    
                   0,  1, -1,    
                   0,  1,  0,

                  -1,  0,  1,    
                   0, -1,  1,    
                   0,  0,  1), ncol=3, byrow=T)
        )


    if( missing(j) ) { 
      return( DMatrix[[J]] ) 
    } else {
      return( (DMatrix[[J]])[ ((j-1)*J +1):(J*j),] )
    }

  } 

  print("Not Yet Implemented")
}


# applyGroup
applyGroup <- function( X, J ) {
   # create a matrix of True and Falses for which category the latent variable falls in
   X.category <- matrix( colSums(matrix(matrix(DMat(J) %*% matrix( X, nrow=J) > 0, nrow=J*(J+1) ),nrow=J)) == J, nrow=J+1)  
   #X.category <<- X.category
   return( unlist( apply(X.category, 2, which) ))
}




# function to evaluate the pdf
sarCheck <- function( rho, Y, X, Beta, W, Sigma, H ) {

  n <- nrow(W)
# this is for completeness, but is VERY slow due to the solve call
  if( missing(H) ) {
    H <- solve(Sigma)
  }

  J <- length(Y)/nrow(W) 

  Lambda <- kronecker(diag(n) - rho * W, diag(J)) 

  Z <- Lambda %*% Y - X %*% Beta 

    
  return( as.numeric( -1 * ( .5 * log( det(Lambda %*% H %*% Lambda)  )    -1/2 * t(Z) %*% H %*%Z ))  )
}




# this function performs mh sampler on either rho1 or rho2,
logpdf.sar <- function(rho,Z,W,Beta,q.value) {

  n <- nrow(W)
  N <- nrow(Z)
  K <- N/n
  
  rho1 <- rho[1]
  rho2 <- rho[2]
  
  M <- matrix(1, nrow=K,ncol=K)
  
  Lambda1 <- diag(n) - rho1 * W
  Lambda2 <- diag(n) - rho2 * W
  
  S1 <- kronecker( diag(K), sqrt(q.value) * solve(Lambda1) )
  
  Sigma1.x <- t(S1) %*% S1
  Sigma2.x <- (1-q.value) * kronecker( M, solve( t(Lambda2) %*% Lambda2 ) )
  
  Sigma.x <- Sigma1 + Sigma2
  Sigma.inv <- solve(Sigma)
  Sigma.det <- det(Sigma)
  
  Q <- S1 %*% X
  
  Z.adj <-  Z - Q %*%Beta
  
  results <- - 1/2 * log(Sigma.det) - 1/2 * t(Z.adj) %*% Sigma.inv %*%  Z.adj
  
  return( results )

}


# generate a single wishart deviate under the constraint that the first diagonal element is 1
rWishart.control <- function(v, V) {
  p <- minDim(V) 

  L <- t(chol(V)) 

  A <- matrix(0,nrow=p,ncol=p)
  A[1,1] <- 1/L[1,1]
  
  i <- 2 
  while( i <= p ) {
  
    A[i,i] <- sqrt(rchisq(1, v + 1 - i)) 
  
    j <- 1 
    while( j < p ) {
      A[i,j] <-rnorm(1)
      j = j + 1
    }
  
    i = i + 1
  }
  
  return( L %*% A %*% t(A) %*% t(L) )
}
  



test.X <- matrix(c( 
  1, 0, 0,
  1, 0, 0,
  0, 1, 0,
  0, 1, 0,
  1, 0, 0,
  1, 0, 0,
  0, 0, 1,
  0, 0, 1
  ), ncol=3, byrow=T)


test.W <- matrix( c(
  0, 1, 1, 0,
  1, 0, 0, 1, 
  1, 0, 0, 1, 
  0, 1, 1, 0 ),ncol=4,byrow=T)



