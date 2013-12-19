
#library(RPostgreSQL)
source('MNPSARTools.R')



calcMCMLE <- function( X, byRow=F ) {
  
  r <- ncol(X)
  p <- nrow(X)
  # error if X doesn't have at least 2 rows
  if( is.null(r) ) stop("Insufficient Number of Columns")
  if( r < 2 ) stop("Insufficient Number of Columns")

  # fun with recursion  
  if( byRow ) {
    X.cols <- lapply( 1:p, function(x) return( matrix( X[x,],nrow=1)))
    return( lapply(X.cols, calcMCMLE) ) 
  }
 
  unique.values <- sort(unique(c(X)))
  J <- length(unique.values)
  
  result <- matrix(0,nrow=J,ncol=J)

  # convert the values in X to there order in unique.values
  Y <- apply(X,c(1,2), function(x) { which( unique.values %in%  x) } )

  Y.table <- table( c(Y[,-r] + Y[,-1]*J ) ) 
  Y.table.names <- as.numeric(names(Y.table))

  for(i in 1:(J*J + J) ) {
    if( i %in% Y.table.names) {
      result[i - J] <- Y.table[ which( Y.table.names %in% i) ]
    }
  }
  
  result <- result / rowSums(result)
  colnames(result) <- unique.values 
  rownames(result) <- unique.values

  return(result)
}







# function to generate a set of crops 
#input 
# number of years
# initial conditions
# Parameters
#  Beta
#  Sigma.list
#  Rho
# burnIn
# Number of Quarter Sections
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
  for(i in 1:(burnIn+years-1) ) {
    a.crops <- sarTools.generateCropTypes(a.crops, rho=rho, Beta=Beta, Sigma.list=list(Sigma.Annual) ) 
  }


  # calculate full transition probabilities
  regionMCMC <- calcMCMLE( a.crops$cropType[,-1] ) 
  
  parcelMCMC <- calcMCMLE( a.crops$cropType[,-1] ,byRow=T) 

  # subset result
  # here we take a subset of 
  #   cropType
  #   cropValue
  keepRange <- (burnIn+3):(burnIn+2+years)

  a.crops$cropType <- a.crops$cropType[,c(1,keepRange)]
  a.crops$cropValue <- a.crops$cropValue[,keepRange - 2]
  
  return(list(crops=a.crops,
              regionMCMC=regionMCMC,
              parcelMCMC=parcelMCMC))
}





set.seed(800)

# ratio of corn and soybeans
p <- c(.4, .4, .2)
years <-  5 

# parameters 
rho <- c(-.15,.20)

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

QQS.size <- c(2,2)

burnIn <- 1000 

b <- generateCropDeviates(
  years,
  init.p,
  Beta,
  Sigma.list,
  Rho,
  burnIn,
  QQS.size
  ) 

print(b)
