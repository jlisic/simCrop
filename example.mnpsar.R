#library(RPostgreSQL)
source('~/src/simCrop/MNPSARTools.R')


parPlot <- function(x ,myTitle.hist, myTitle.acf, myFile) {
  pdf(myFile,width=7,height=.6*7 )
  par(mfrow=c(1,2))
  hist(x, main=myTitle.hist)
  acf(x,main=myTitle.acf)
  dev.off()
}

parPlot.crops <- function(a , myFile) {
  pdf(myFile)
  par(mfrow=c(1,1))
  #plot( raster(a$globalError[ sort(a$globalError[,1],index=T)$ix,2]) )
  simCrop.plotCropTypes(a,newDev=F)
  dev.off()
}




# In this example we are considering just two crops
# 1. corn
# 2. soybeans

set.seed(800)

# ratio of corn and soybeans
p <- c(.4, .4, .2)
J <- length(p) - 1
K <- 2 

# parameters 
rho <- c(-.15,.20)
Beta1.sim.corn  <- 1 
Beta2.sim.corn  <- 1 
Beta1.sim.soy   <- -1 
Beta2.sim.soy   <- -1 
Beta1.sim.other <- -0.5 
Beta2.sim.other <- -0.5 
Beta <- matrix( c( 
                  Beta1.sim.corn, 
                  Beta2.sim.corn, 
                  Beta1.sim.soy, 
                  Beta2.sim.soy, 
                  Beta1.sim.other,
                  Beta2.sim.other
                  ),ncol=1)

iter <- 500 
thinning <- 1 
burnIn <- 0 
m <- 20 

Sigma.Annual <- matrix( c(1,0,0,1), nrow=2,byrow=T)
Sigma.Environment <- matrix( c(1,0,0,1)/10, nrow=2,byrow=T)

# create a 2x2 section set of quarter-quarter sections (QQS)
a <- simCrop.partitionPLSS(1,1)

# add initial crop assignment
a.init <- simCrop.generateCropTypes(a,p)
a.init <- simCrop.getNeighbors(a.init)
W <- simCrop.createRookDist(a.init) 

# simulate 10 years of data

a.crops <- sarTools.generateCropTypes(a.init, rho=rho, Beta=Beta, Sigma.list= list(Sigma.Annual, Sigma.Environment) ) 
for(i in 2:K) {
  a.crops <- sarTools.generateCropTypes(a.crops, rho=rho, Beta=Beta, Sigma.list=list(Sigma.Annual) ) 
}


### useful for debugging

Z <- sortCrop(a.crops$cropValue, a.crops)
V <- sortCrop(a.crops$globalError[,'error'], a.crops)
n <- nrow(W)

sarTools.probitGibbsSpatial( a.crops, function(...) {} ) 

Lambda1 <- kronecker( diag(n) - rho[1] * W, diag(J))
Lambda2 <- kronecker( diag(n) - rho[2] * W, diag(J)) 
Lambda1.K <- kronecker( diag(K), Lambda1) 
Lambda1.K.inv <- kronecker( diag(K), solve(Lambda1)) 
  
FullSigmaRoot <-  kronecker(diag(K), bdiag( rep( list( solve(Sigma.Annual) ), times=n) )) 

Beta.hat <- solve( t(X) %*%FullSigmaRoot %*% X )  %*% 
  t(X) %*% FullSigmaRoot %*% Lambda1.K %*% 
  ( matrix(Z,ncol=1) - Lambda1.K.inv %*% kronecker( matrix(1,ncol=1,nrow=K), V) )  


##### inits
beta.init <- Beta 
rho.init <- c( -.15, .20)
Sigma.init <- list( Sigma.Annual, Sigma.Environment)

alpha.init <- V 
Z.init <-  matrix(Z,ncol=1)

##### hyper params
Beta0 <- c(0,0,0,0,0,0) 
Sigma0 <- rep(10,times=6) 
Wishart0 <- list(1,diag(2)) 



#
if ( T ) {
result <- sarTools.probitGibbsSpatial( 
  a.crops, 
  fun=sarTools.probitGibbsSpatialRunConditional,

  beta.init=beta.init,
  rho.init=rho.init,
  Z.init=Z.init,
  alpha.init=alpha.init,

##hyper parameters
  Beta0=Beta0,
  Sigma0=Sigma0,
  Wishart0.Z=Wishart0, # shape and rate
  Wishart0.alpha=Wishart0, # shape and rate
##runtime params
  iter=iter,
  m=m,         # thinning for Y
  thinning=thinning,  # thinning
  burnIn=burnIn     # burnIn
  )

print( colMeans(result$Beta))
print( colMeans(result$Rho))
print( colMeans(result$Sigma1))
print( colMeans(result$Sigma2))

par( mfrow=c(3,2) )
for( i in 1:length(Beta) ) {
  plot( result$Beta[,i],type='l',  main=sprintf("Beta %d",i ) )
}

}
#
