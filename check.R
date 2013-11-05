
source('sarTools.R')

# In this example we are considering just two crops
# 1. corn
# 2. soybeans

# ratio of corn and soybeans
p <- c(1, 0)

# create a 2x2 section set of quarter-quarter sections (QQS)
a <- simCrop.partitionPLSS(2,2)

# add initial crop assignment
a.init <- simCrop.generateCropTypes(a,p)
a.init <- simCrop.getNeighbors(a.init)
W <- simCrop.createRookDist(a.init) 

  # parameters for the CAR model
#rhoFromCorn   <- min( carTools.checkRho( W ) ) + 0.01
#rhoFromSoy    <- max( carTools.checkRho( W ) ) - 0.01
rhoFromCorn <- -.15
rhoFromSoy <- -0.05 
BetaFromCorn  <-  0 
BetaFromSoy   <-  0 

rho <- list(rhoFromCorn,rhoFromSoy)
Beta <- list( BetaFromCorn, BetaFromSoy)

# add some simulated subsequent years
n <- length(a$map[,'object'])
X <- matrix(1,nrow=n,ncol=1)

#simulate 10 years of data
#a.crops <- a.init
#for(i in 1:1) {
#  a.crops <- sarTools.generateCropTypes(a.crops, rho=rho, X=X, Beta=Beta) 
#}




rhoRange <-carTools.checkRho(W) 
Beta.result <- c()
rhoMean.result <- c()
rhoSD.result <- c()
probit.Beta.result <- c()
probit.rhoMean.result <- c()
probit.rhoSD.result <- c()

for( i in 1:30) {
  Y <- sarTools.deviates( rhoFromCorn, W, X, BetaFromCorn) 
  
  result <- sarTools.probitGibbsSpatialRun( matrix(Y>0,ncol=1), X,W, 0,0,0,3,100,10)
  Beta.result[i] <- mean(result$Beta)
  rhoSD.result[i] <- sd(result$Beta)
  rhoMean.result[i] <- mean(result$rho)
  rhoSD.result[i] <- sd(result$rho)
  
  result <- sarTools.gibbsSpatialRun( Y, X,W, 0,0,0,3,100,10)
  probit.Beta.result[i] <- mean(result$Beta)
  probit.rhoSD.result[i] <- sd(result$Beta)
  probit.rhoMean.result[i] <- mean(result$rho)
  probit.rhoSD.result[i] <- sd(result$rho)
  
  opt.result[i] <- optim( rhoFromCorn, sarCheck, Y=Y, W=W, X=X, Beta=BetaFromCorn, lower=rhoRange[1] + 0.0001, upper=rhoRange[2] - 0.0001,method="Brent" )$par
  
}

print( mean( Beta.result))
print( mean( rhoMean.result ))
print( mean( rhoSD.result ))

print( mean( probit.Beta.result))
print( mean( probit.rhoMean.result ))
print( mean( probit.rhoSD.result ))

print( mean( opt.result ) )

#opt.result <- optim( rhoFromCorn, sarCheck, Y=Y, W=W, X=X, Beta=BetaFromCorn, lower=rhoRange[1] + 0.0001, upper=rhoRange[2] - 0.0001,method="Brent" )$par
#print( opt.result )

#opt.result <- c()
#for( i in 1:10) {
#  print(i)
#Y <- sarTools.deviates( rhoFromCorn, W, X, BetaFromCorn) 
#opt.result[i] <- optim( rhoFromCorn, sarCheck, Y=Y, W=W, X=X, Beta=BetaFromCorn, lower=rhoRange[1] + 0.0001, upper=rhoRange[2] - 0.0001,method="Brent" )$par
#}
#
#print(mean(opt.result))
#


