# function to generate a set of crops 

#library(RPostgreSQL)
source('MNPSARTools.R')

#input 
# number of years
# initial conditions
# Parameters
#  Beta
#  Sigma.list
#  Rho
# burnIn
# Number of Quarter Sections

set.seed(800)

# ratio of corn and soybeans
p <- c(.4, .4, .2)
years <- 20 

# parameters 
rho <- c(-.15,.20)

Beta1.sim.corn  <- 1 
Beta2.sim.corn  <- 1 
Beta1.sim.soy   <- -1 
Beta2.sim.soy   <- -1 
Beta1.sim.other <- -0.5 
Beta2.sim.other <- -0.5 
Beta <- matrix( c( Beta1.sim.corn, Beta2.sim.corn, Beta1.sim.soy, Beta2.sim.soy, Beta1.sim.other, Beta2.sim.other),ncol=1)

Sigma.Annual <- matrix( c(1,0,0,1), nrow=2,byrow=T)
Sigma.Environment <- matrix( c(1,0,0,1)/10, nrow=2,byrow=T)

Sigma.list <- list( Sigma.Annual, Sigma.Environment)

QQS.size <- c(1,1)

generateCropDeviates <- function(
  years,
  init.p,
  Beta,
  Sigma.list,
  Rho,
  burnIn,
  QQS.size
  ) {

  # create a 2x2 section set of quarter-quarter sections (QQS)
  a <- simCrop.partitionPLSS(QQS.size[1],QQS.size[2])
  
  # add initial crop assignment
  a.init <- simCrop.generateCropTypes(a,p)
  a.init <- simCrop.getNeighbors(a.init)
  W <- simCrop.createRookDist(a.init) 
  
  # simulate years of data
  a.crops <- sarTools.generateCropTypes(a.init, rho=rho, Beta=Beta, Sigma.list= list(Sigma.Annual, Sigma.Environment) ) 
  for(i in 2:burnIn) {
    a.crops <- sarTools.generateCropTypes(a.crops, rho=rho, Beta=Beta, Sigma.list=list(Sigma.Annual) ) 
  }

  # subset resulf
  
  return(a.crops)
}

