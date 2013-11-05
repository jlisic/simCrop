#[1] "Simulating with Beta=0.000000 Rho=0.150000"
#[1] "MLE Beta=0.000000 Rho=0.164480"
#[1] "Simulating with Beta=0.000000 Rho=0.100000"
#[1] "MLE Beta=0.000000 Rho=0.077363"
#



source('sarTools.R')

# In this example we are considering just two crops
# 1. corn
# 2. soybeans

# ratio of corn and soybeans
p <- c(.5, .5)

# create a 2x2 section set of quarter-quarter sections (QQS)
a <- simCrop.partitionPLSS(4,4)

# add initial crop assignment
a.init <- simCrop.generateCropTypes(a,p)
a.init <- simCrop.getNeighbors(a.init)
W <- simCrop.createRookDist(a.init) 

  # parameters for the CAR model
#rhoFromCorn   <- min( carTools.checkRho( W ) ) + 0.01
#rhoFromSoy    <- max( carTools.checkRho( W ) ) - 0.01

rhoFromCorn <- .15
rhoFromSoy <- .10

BetaFromCorn  <-  0 
BetaFromSoy   <-  0 

rho <- list(rhoFromCorn,rhoFromSoy)
Beta <- list( BetaFromCorn, BetaFromSoy)

# add some simulated subsequent years
n <- length(a$map[,'object'])
X <- matrix(1,nrow=n,ncol=1)

#simulate 10 years of data
a.crops <- a.init
for(i in 1:1) {
  a.crops <- sarTools.generateCropTypes(a.crops, rho=rho, X=X, Beta=Beta) 
}

rhoRange <-carTools.checkRho(W) 

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

iter <- 400

result <- sarTools.probitGibbsSpatial(a.crops,Beta.init,rho.init,Beta0,Sigma0,iter)

print( "corn")
print( mean( result[[1]]$rho) )
print( sd( result[[1]]$rho) )
print( mean( result[[1]]$Beta) )

print( "soy")
print( mean( result[[2]]$rho) )
print( sd( result[[2]]$rho) )
print( mean( result[[2]]$Beta) )
