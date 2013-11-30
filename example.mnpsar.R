#library(RPostgreSQL)
source('MNPSARTools.R')


parPlot <- function(x ,myTitle.hist, myTitle.acf, myFile) {
  pdf(myFile,width=7,height=.6*7 )
  par(mfrow=c(1,2))
  hist(x, main=myTitle.hist)
  acf(x,main=myTitle.acf)
  dev.off()
}

parPlot.crops <- function(a , myFile) {
  pdf(myFile)
  par(mfrow=c(3,3))
  #plot( raster(a$globalError[ sort(a$globalError[,1],index=T)$ix,2]) )
  simCrop.plotCropTypes(a,newDev=F)
  dev.off()
}




# In this example we are considering just two crops
# 1. corn
# 2. soybeans

set.seed(400)

# ratio of corn and soybeans
p <- c(.4, .4, .2)
J <- length(p) - 1

# parameters 
rho <- c(-.15,.20)
Beta.sim.corn  <- -1 
Beta.sim.soy   <- 2 
Beta.sim.other <- 1 
Beta <- matrix( c( Beta.sim.corn, Beta.sim.soy, Beta.sim.other ),ncol=1)

iter <- 30 
thinning <- 20
burnIn <- 0 
m <- 10 

Sigma.Annual <- matrix( c(1,0,0,1), nrow=2,byrow=T)
Sigma.Environment <- matrix( c(1,0,0,1), nrow=2,byrow=T)

# create a 2x2 section set of quarter-quarter sections (QQS)
a <- simCrop.partitionPLSS(1,1)

# add initial crop assignment
a.init <- simCrop.generateCropTypes(a,p)
a.init <- simCrop.getNeighbors(a.init)
W <- simCrop.createRookDist(a.init) 

# simulate 10 years of data

a.crops <- sarTools.generateCropTypes(a.init, rho=rho, Beta=Beta, Sigma.list= list(Sigma.Annual, Sigma.Environment) ) 
for(i in 2:5) {
  a.crops <- sarTools.generateCropTypes(a.crops, rho=rho, Beta=Beta, Sigma.list=list(Sigma.Annual) ) 
}


## useful for debugging
myObjects <- a.crops$cropType[,'myObjects']
object.sort <- sort(myObjects,index.return=T)$ix

Z <- matrix( t(a.crops$cropValue), ncol=length(object.sort) )[,object.sort]
Z <- t( matrix( Z, ncol=ncol(a.crops$cropValue)*J ) )
Z <- matrix( Z, ncol=1)

V <- matrix( a.crops$globalError[,'error'], ncol=length(object.sort) )[,object.sort]
V <- matrix( V, ncol=1)

V2 <- rep(V,length=length(Z))


Z <- Z - V2

#### inits
##Beta.init <- c(0,0)
beta.init <- Beta 
rho.init <- c( -.15, .20)
Sigma.init <- list( Sigma.Annual, Sigma.Environment)
alpha.init <- matrix(0,nrow=nrow(W)*minDim(Sigma.init[[2]]),ncol=1)
Z.init <- matrix(0,nrow=nrow(W)*minDim(Sigma.init[[1]])*(ncol(a.crops$cropType) - 2),ncol=1)
#
#### hyper params
Beta0 <- c(0,0,0) 
Sigma0 <- c(2,2,2) 
Gamma0 <- c(1,1/2) # not likely to be used


#
if ( T ) {
result <- sarTools.probitGibbsSpatial( 
  a.crops, 
  fun=sarTools.probitGibbsSpatialRunConditional,

  beta.init=beta.init,
  rho.init=rho.init,
  Z.init=Z.init,
  alpha.init=alpha.init,
  Sigma.init=Sigma.init,

##hyper parameters
  Beta0=Beta0,
  Sigma0=Sigma0,
  Gamma0=Gamma0, # shape and rate
##runtime params
  iter=iter,
  m=m,         # thinning for Y
  thinning=thinning,  # thinning
  burnIn=burnIn     # burnIn
  )

print( colMeans(result$Beta))
#print( colMeans(result$Rho))
#print( colMeans(result$Tau))
}
#
