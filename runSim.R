#library(RPostgreSQL)
source('MNPSARTools.R')
source('plotsarmnp.R')
source('generate.crops.R')




set.seed(800)

# ratio of corn and soybeans
init.p <- c(.4, .4, .2)
years <-  10 
years.params <- 3

burnIn <- 20 
burnIn.params <- 497 

m.params <- years.params + burnIn.params

J <- length(p) - 1

# parameters 
Rho <- c(-.15,.20)

Beta1.sim.corn  <- .5 
Beta2.sim.corn  <- -0.5 
Beta1.sim.soy   <- -1.5 
Beta2.sim.soy   <- -1.5 
Beta1.sim.other <- -0.5 
Beta2.sim.other <- -0.5 
Beta <- matrix( c( Beta1.sim.corn, Beta2.sim.corn, Beta1.sim.soy, Beta2.sim.soy, Beta1.sim.other, Beta2.sim.other),ncol=1)

Sigma.Annual <- matrix( c(1,0,0,1), nrow=2,byrow=T)
Sigma.Environment <- matrix( c(1,0,0,1)/10, nrow=2,byrow=T)

Sigma.list <- list( Sigma.Annual, Sigma.Environment)

QQS.size <- c(1,1)




# control switches
createParameters <- T 
runRegionMCMC <- T
runParcelMCMC <- T
runSpatialProbit <- T 







################################# PARAMETERS ##########################################
# generate parameters

if( createParameters ) {

  params <- generateCropDeviates(
    years.params,
    init.p,
    Beta,
    Sigma.list,
    Rho,
    burnIn.params,
    QQS.size,
    densityPlot=F
    ) 

  # "true" parameters
  params.parcel <- unlist(params$parcelMCMC)

}






################################# TEST #############################################
# generate a test data set

test <- generateCropDeviates(
  years,
  init.p,
  Beta,
  Sigma.list,
  Rho,
  burnIn,
  QQS.size,
  densityPlot=F
  ) 



# Estimate for 5 years, 10 years, 20 years
test.crop <- test$crops



##################  MCMC total estimation ##############################
if ( runRegionMCMC ) {
  test.regionMCMC <- calcMCMLE( test.crop$cropType[,-1], unique.values=test.crop$crops ) 
  test.regionMCMC.est <- rep(c(test.regionMCMC),times=nrow(test.crop$cropType))
}

##################  MCMC Bayesian estimation ###########################
if ( runParcelMCMC ) {
  test.parcelMCMC <- calcMCMLE( test.crop$cropType[,-1] ,byRow=T,useBayes=T, unique.values=test.crop$crops) 
  test.parcelMCMC.est <- unlist(test.parcelMCMC)
}

###################  Spatial Probit  ###################################
if ( runSpatialProbit ) {

  ##### inits
  beta.init  <- c(0, 0, 0, 0, 0, 0) 
  rho.init   <- c( 0, 0)
  Sigma.init <- list( diag(2), diag(2) )
  
  alpha.init <- matrix( 0, nrow=nrow(test.crop$cropValue), ncol=1)
  Z.init     <- matrix( 0, nrow=nrow(test.crop$cropValue) * (years - 1), ncol=1)
  
  ### hyper params
  Beta0 <- c(0,0,0,0,0,0) 
  Sigma0 <- rep(10,times=6) 
  Wishart0 <- list(1,diag(2)) 
  
  iter <- 3000 
  thinning <- 0 
  burnIn <- 0 
  m <- 20 
 
  years.SARMNP <- 3
  burnIn.SARMNP <- 497 
  m.SARMNP <- years.SARMNP + burnIn.SARMNP


  ###run model 
  test.SARMNP <- sarTools.probitGibbsSpatial( 
    test.crop, 
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

  print( colMeans(test.SARMNP$Beta))
  print( colMeans(test.SARMNP$Rho))
  print( 1/colMeans(test.SARMNP$Sigma1))
  print( 1/colMeans(test.SARMNP$Sigma2))

#  x11()
#  par( mfrow=c(3,2) )
#  for( i in 1:length(Beta) ) {
#    plot( test.SARMNP$Beta[,i],type='l',  main=sprintf("Beta %d",i ) )
#  }


  test.SARMNP.Z <- unsortCrop( test.SARMNP$Z, test.crop) 
  test.SARMNP.Z.type <- apply( test.SARMNP.Z,2, applyGroup, J=J )
  test.SARMNP.MCMC <- calcMCMLE( test.SARMNP.Z.type, unique.values=test.crop$crops ) 

  test.SARMNP.est.gen <- generateCropDeviates(
    years.SARMNP,
    init.p,
    Beta,
    Sigma.list,
    Rho,
    burnIn.SARMNP,
    QQS.size,
    densityPlot=F
  ) 
  
  test.SARMNP.est <- unlist(test.SARMNP.est.gen$parcelMCMC)
  
}


MSE.SARMNP <-  sum( (test.SARMNP.est - params.parcel)^2/m.params ) 
MSE.regionMCMC   <-  sum( (test.regionMCMC.est - params.parcel)^2/m.params ) 
MSE.parcelMCMC   <-  sum( (test.parcelMCMC.est - params.parcel)^2/m.params ) 






