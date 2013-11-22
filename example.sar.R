#library(RPostgreSQL)
source('sarTools.R')


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
p <- c(.6, .4)

# parameters 
rho <- c(-.15,.20)
Beta.sim.corn  <- -1 
Beta.sim.soy   <- 2 
Beta <- matrix( c( Beta.sim.corn, Beta.sim.soy),ncol=1)

iter <-  3 
thinning <- 1
burnIn <- 0 
m <- 1
tau <- c(1,.5)


# create a 2x2 section set of quarter-quarter sections (QQS)
a <- simCrop.partitionPLSS(1,1)

# add initial crop assignment
a.init <- simCrop.generateCropTypes(a,p)
a.init <- simCrop.getNeighbors(a.init)
W <- simCrop.createRookDist(a.init) 

# simulate 10 years of data
a.crops <- sarTools.generateCropTypes(a.init, rho=rho, Beta=Beta, tau=tau) 
for(i in 2:2) {
  a.crops <- sarTools.generateCropTypes(a.crops, rho=rho, Beta=Beta, tau=tau[1]) 
}


# useful for debugging
object.sort <- sort(a.crops$cropType[,'myObjects'],index.return=T)$ix
Z <- c(a.crops$cropValue[object.sort,])

Beta.init <- c(-1,2)
rho.init <- c( -.15, .20)
alpha.sigma.init <- matrix(0,nrow(W),ncol=1)
tau.init <- .5 


Beta0 <- c(0,0) 
Sigma0 <- c(2,2) 
Gamma0 <- c(1,1/2)


result <- sarTools.probitGibbsSpatial( 
  a.crops, 
  fun=sarTools.probitGibbsSpatialRunConditional,

  Beta.init=Beta.init,
  rho.init=rho.init,
  alpha.sigma.init=alpha.sigma.init,
  tau.init=tau.init,

#hyper parameters
  Beta0=Beta0,
  Sigma0=Sigma0,
  Gamma0=Gamma0, # shape and rate
#runtime params
  iter=iter,
  m=m,         # thinning for Y
  thinning=thinning,  # thinning
  burnIn=burnIn     # burnIn
  )

print( colMeans(result$Beta))
print( colMeans(result$Rho))
print( colMeans(result$Tau))


if ( F ) {
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
  burnIn,   # burnIn
  method="omnibus"
  )

print( colMeans(result$Beta))
print( colMeans(result$rho))
print( colMeans(result$q.value))

parPlot( result$Beta[,1], "Beta, From Corn Histogram", "Beta, From Corn ACF", "~/src/dissertation/tSnippets/BetaCorn.pdf")
parPlot( result$Beta[,2], "Beta, From Soybeans Histogram", "Beta, From Soybeans ACF", "~/src/dissertation/tSnippets/BetaSoybeans.pdf")
parPlot( result$rho[,1], "rho1 Histogram", "rho1 ACF", "~/src/dissertation/tSnippets/rho1.pdf")
parPlot( result$rho[,2], "rho2 Histogram", "rho2 ACF", "~/src/dissertation/tSnippets/rho2.pdf")
parPlot( result$q.value, "q Histogram", "q ACF", "~/src/dissertation/tSnippets/qvalue.pdf")

}





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

## test for mh
if( F ) {
result <- sarTools.probitGibbsSpatial(a.crops, fun=function(...) {} ) 

thinning <- 25 
rhos <- c(0)
start.time <- proc.time()
for(i in 2:1000) {
  rhos[i] <- mh.lambda.sar(Z=Z,W=W,mu=X %*% Beta,tau=1,x0=rhos[i-1],iter=1, burnIn=thinning, rho.range=carTools.checkRho(W)) 
}
print( proc.time() - start.time)
plot(rhos)
}


