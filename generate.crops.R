
#library(RPostgreSQL)
source('MNPSARTools.R')
source('plotsarmnp.R')


  table.nomiss <- function( x, cats ) {
    y <- table(x)
    y.cats <- rep(0,times=length(cats))
    y.cats[ as.numeric(names(y)) ] <- y 
    return( y.cats )
  }


calcMCMLE <- function( X, byRow=F, useBayes=F, unique.values ) {
  
  r <- ncol(X)
  p <- nrow(X)
  # error if X doesn't have at least 2 rows
  if( is.null(r) ) stop("Insufficient Number of Columns")
  if( r < 2 ) stop("Insufficient Number of Columns")



  # fun with recursion  
  if( byRow ) {
    X.cols <- lapply( 1:p, function(x) return( matrix( X[x,],nrow=1)))
    if( missing( unique.values) ) { 
      return( lapply(X.cols, calcMCMLE, useBayes=useBayes) ) 
    } else {
      return( lapply(X.cols, calcMCMLE, useBayes=useBayes, unique.values=unique.values) ) 
    }
  }

  # allow forcing values 
  if( missing( unique.values) ) { 
    unique.values <- sort(unique(c(X)))
  } 
    
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

  # use a non-informative prior
  if( useBayes ) { 
    result <- result + 1 / ( rowSums(result) + J ) 
  } else {
    result <- result / rowSums(result)
  }

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
  p,
  Beta,
  Sigma.list,
  Rho,
  burnIn,
  QQS.size,
  densityPlot,
  ratioPlot=F
  ) {

  # create a 2x2 section set of quarter-quarter sections (QQS)
  a <- simCrop.partitionPLSS(QQS.size[1],QQS.size[2])
  
  # add initial crop assignment
  a.init <- simCrop.generateCropTypes(a,p)
  a.init <- simCrop.getNeighbors(a.init)
  W <- simCrop.createRookDist(a.init) 
  
  # simulate years of data
  a.crops <- sarTools.generateCropTypes(a.init, rho=Rho, Beta=Beta, Sigma.list= list(Sigma.Annual, Sigma.Environment) ) 
  for(i in 1:(burnIn+years-1) ) {
    a.crops <- sarTools.generateCropTypes(a.crops, rho=Rho, Beta=Beta, Sigma.list=list(Sigma.Annual) ) 
  }


  # calculate full transition probabilities
  regionMCMC <- calcMCMLE( a.crops$cropType[,-1], unique.values=a.crops$crops ) 
  parcelMCMC <- calcMCMLE( a.crops$cropType[,-1] ,byRow=T,useBayes=T, unique.values=a.crops$crops) 

  #get crop ratios

  cropRatio <- t(apply(a.crops$cropType[,-1],2,table.nomiss, a.crops$crops))/nrow(a.crops$cropType)
  rownames(cropRatio) <- 0:(nrow(cropRatio) - 1)

  if( ratioPlot != F ) { 
    pdf( file=ratioPlot )
    plot(cropRatio[,2],type='o',col=1,ylim=c(0,1))
    lines(cropRatio[,3],type='o',col=2)
    lines(cropRatio[,1],type='o',col=3)
    legend(1, legend=c(1:3), cex=0.8, col=c(1:3), pch=21:22, lty=1:2)
    dev.off()
  }
  # create a density plot
  if(densityPlot) plotSARMNP( a.crops, Beta=Beta, plotPoints=F,removeInitial=T ) 

  # subset result
  # here we take a subset of 
  #   cropType
  #   cropValue
  keepRange <- (burnIn+3):(burnIn+2+years)

  a.crops$cropType <- a.crops$cropType[,c(1,keepRange)]
  a.crops$cropValue <- a.crops$cropValue[,keepRange - 2]
  
  return(list(crops=a.crops,
              regionMCMC=regionMCMC,
              parcelMCMC=parcelMCMC,
              cropRatio=cropRatio
              ))
}




