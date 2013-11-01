
source('carTools.R')

# In this example we are considering just two crops
# 1. corn
# 2. soybeans

# ratio of corn and soybeans
p <- c(.6, .4)

# create a 2x2 section set of quarter-quarter sections (QQS)
a <- simCrop.partitionPLSS(2,2)

# add initial crop assignment
a.init <- simCrop.generateCropTypes(a,p)
a.init <- simCrop.getNeighbors(a.init)
W <- simCrop.createRookDist(a.init) 

  # parameters for the CAR model
rhoFromCorn   <- min( carTools.checkRho( W ) ) + 0.0001
rhoFromSoy    <- max( carTools.checkRho( W ) )
BetaFromCorn  <-  0 
BetaFromSoy   <-  0 

rho <- list(rhoFromCorn,rhoFromSoy)
Beta <- list( BetaFromCorn, BetaFromSoy)

# add some simulated subsequent years
n <- length(a$map[,'object'])
X <- matrix(1,nrow=n,ncol=1)

#simulate 10 years of data
a.crops <- a.init
for(i in 1:10) {
  a.crops <- carTools.generateCropTypes(a.crops, rho=rho, X=X, Beta=Beta) 
}

## probit Gibbs function ## 
Beta.init.corn <- 0
Beta.init.soy <- 0
Beta.init <- list( Beta.init.corn, Beta.init.soy)

rho.init.corn <- 0
rho.init.soy <- 0
rho.init <- list( rho.init.corn, rho.init.soy )

Beta0.corn <- 0
Beta0.soy <- 0
Beta0 <- list( Beta0.corn, Beta0.soy)

Sigma0.corn <- 10
Sigma0.soy <- 10
Sigma0 <- list( Sigma0.corn, Sigma0.soy)

iter <- 250 

result <- carTools.probitGibbsSpatial(a.crops,Beta.init,rho.init,Beta0,Sigma0,iter)

