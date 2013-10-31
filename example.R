
source('carTools.R')

# In this example we are considering just two crops
# 1. corn
# 2. soybeans

# ratio of corn and soybeans
p <- c(1, 0)

# parameters for the CAR model
rhoFromCorn   <- min( carTools.checkRho( W ) ) +0.00001
rhoFromSoy    <- max( carTools.checkRho( W ) )
BetaFromCorn  <-  0 
BetaFromSoy   <-  1.00

rho <- list(rhoFromCorn,rhoFromSoy)
Beta <- list( BetaFromCorn, BetaFromSoy)

# create a 2x2 section set of quarter-quarter sections (QQS)
a <- simCrop.partitionPLSS(2,2)

# add initial crop assignment
a.init <- simCrop.generateCropTypes(a,p)

# add some simulated subsequent years
n <- length(a$map[,'object'])
X <- matrix(1,nrow=n,ncol=1)
a.crops <- carTools.generateCropTypes(a.init, rho=rho, X=X, Beta=Beta) 


