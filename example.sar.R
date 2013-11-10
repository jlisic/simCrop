


source('sarTools.R')

# In this example we are considering just two crops
# 1. corn
# 2. soybeans

# ratio of corn and soybeans
p <- c(.5, .5)

# parameters 
rho <- .15
Beta.sim.corn  <-  -2 
Beta.sim.soy   <-   2 
Beta <- matrix( c( Beta.sim.corn, Beta.sim.soy),ncol=1)
iter <- 1000 

# create a 2x2 section set of quarter-quarter sections (QQS)
a <- simCrop.partitionPLSS(4,4)

# add initial crop assignment
a.init <- simCrop.generateCropTypes(a,p)
a.init <- simCrop.getNeighbors(a.init)
W <- simCrop.createRookDist(a.init) 

# simulate 10 years of data
a.crops <- a.init
for(i in 1:10) {
  a.crops <- sarTools.generateCropTypes(a.crops, rho=rho, Beta=Beta) 
}

#
rhoRange <-carTools.checkRho(W) 
print(sprintf("The range of possible values for rho is (%f, %f)",rhoRange[1], rhoRange[2]))


## probit Gibbs function ## 
Beta.init.corn <- -2 
Beta.init.soy <-   2
Beta.init <- matrix( c(Beta.init.corn, Beta.init.soy), ncol=1) 

rho.init <- .15 

Beta0.corn <- 0
Beta0.soy <- 0
Beta0 <- matrix( c(Beta0.corn, Beta0.soy), ncol=1)

Sigma0.corn <- 3 
Sigma0.soy <- 3 
Sigma0 <- matrix( c(Sigma0.corn, Sigma0.soy), ncol=1)

result <- sarTools.probitGibbsSpatial2(a.crops,Beta.init,rho.init,Beta0,Sigma0,iter,m=10,thinning=20,burnin=200)

#if ( T ) {
#iter <- 100 
#result <- sarTools.probitGibbsSpatial2(a.crops,Beta.init,rho.init,Beta0,Sigma0,iter,50)
#print( "Beta")
#print( colMeans( result$Beta ) )
#print( apply( result$Beta ,2,sd ) )
#
#print( "Rho")
#print( colMeans( result$rho ) )
#print( sd( result$rho) )

#}



