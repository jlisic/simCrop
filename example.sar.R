#library(RPostgreSQL)
source('sarTools.R')

# In this example we are considering just two crops
# 1. corn
# 2. soybeans


# ratio of corn and soybeans
p <- c(.6, .4)

# parameters 
rho.global <- .20 
rho <- -.15
Beta.sim.corn  <- -1 
Beta.sim.soy   <- 2 
Beta <- matrix( c( Beta.sim.corn, Beta.sim.soy),ncol=1)
iter <- 130
thinning <- 10
burnIn <- 30
m <- 10
q.value <- .5



# create a 2x2 section set of quarter-quarter sections (QQS)
a <- simCrop.partitionPLSS(1,1)

# add initial crop assignment
a.init <- simCrop.generateCropTypes(a,p)
a.init <- simCrop.getNeighbors(a.init)
W <- simCrop.createRookDist(a.init) 

# simulate 10 years of data
a.crops <- sarTools.generateCropTypes(a.init, rho=rho, Beta=Beta, rho.global=rho.global, q.value=q.value) 
for(i in 2:5) {
  a.crops <- sarTools.generateCropTypes(a.crops, rho=rho, Beta=Beta, q.value=q.value) 
}


# useful for debugging
object.sort <- sort(a.crops$cropType[,'myObjects'],index.return=T)$ix
Z <- c(a.crops$cropValue[object.sort,])

Beta.init <- c(0,0)
rho.init <- c( 0, 0)
q.init <- q.value 
beta0 <- c(0,0) 
Sigma0 <- c(3,3) 

# result
result <- sarTools.probitGibbsSpatial( 
  a.crops, 
  Beta.init,
  rho.init,
  q.init,
  beta0,
  Sigma0,
  iter,
  m,
  thinning,  # thinning
  burnIn   # burnIn
  )


print( colMeans(result$Beta))
print( colMeans(result$rho))

if( F ) {
#
rhoRange <-carTools.checkRho(W) 
print(sprintf("The range of possible values for rho is (%f, %f)",rhoRange[1], rhoRange[2]))


## probit Gibbs function ## 
Beta.init.corn <- 0 
Beta.init.soy <-  0 
Beta.init <- matrix( c(Beta.init.corn, Beta.init.soy), ncol=1) 

rho.init <- 0 

result.can <- sarTools.probitGibbsSpatial2(a.crops,Beta.init,rho.init,Beta0,Sigma0,iter,m=10,thinning=20,burnin=200)
}




#if( F) {
## load driver for postgresql
#drv <- dbDriver("PostgreSQL")
# 
#  con <- dbConnect( drv, host="10.0.1.2", dbname="edges", user="pgsql")
#
#  # write our results
#  dbWriteTable(con, name="results", value=x, overwrite=FALSE)
#
#  # get the result back
#  dbDisconnect(con)
#
#dbUnloadDriver(drv)
#}


