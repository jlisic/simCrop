
set.seed(1)

library(raster)
library(rgdal)
library(clv)


# simCrop.R

# Township - 23040
# Section - 640
# Half-Section - 320
# Quarter-Section - 160
# Half-Quarter Section - 80
# Quarter-Quarter Section - 40
# assume approximately 56mx56m pixels so that one quarter section contains 49 pixels


# function to convert quarter sections to half quarter sections
# x value
q2hq <- function (x ) { 
  return( c(x*2,x*2-1) )
 }


# function to convert half quarter sections to quarter quarter sections
# x value
# w width
hq2qq <- function( a, w ) {

  y <- ceiling(a/w)
  x <- ((a-1) %% w) +1 

  return( c( (y-1)*w*2 + x, (y-1)*w*2 + w + x ) )
}


# This function returns a mapping of PLSS quarter-quarter sections to quarter-quarter half-quarter and quarter sections
simCrop.partitionPLSS <- function( nsSection=1, weSection=1, landRatioOfQuarterQuarterSections=1 ,landRatioOfHalfQuarterSections=0){
#nsSection=1
#weSection=1
#landRatioOfQuarterQuarterSections=.50 
#landRatioOfHalfQuarterSections=.25

  # we first start off with establishing the number of quarter sections, and the number
  # of these sections that need to be converted into smaller units

  qsCount <- nsSection * weSection * 4 

  # get how many quarter sections to partition
  partitionQSCount <- floor( (landRatioOfHalfQuarterSections + landRatioOfQuarterQuarterSections) * qsCount ) 
  # figure out which quarter sections to partition
  partitionQS <- sample(1:qsCount, size=partitionQSCount)

  # create a ref map for QQS -> QS 
  useAsQS <- (1:qsCount)[!(1:qsCount %in% partitionQS)]
  useAsQS <- q2hq(useAsQS)
  useAsQS <- hq2qq(useAsQS, 4 * weSection)
 
  refAsQS <- rep(useAsQS[1:(length(useAsQS)/4)],times=4)

  # if we don't need to partition any thing we just skip the partitioning.
  if( partitionQSCount > 0 ) {

    #now turn partitionQS into HQS
    partitionHQS <- q2hq(partitionQS)
  
    #get our HQS
    partitionHQSCount <- floor( landRatioOfHalfQuarterSections / (landRatioOfHalfQuarterSections + landRatioOfQuarterQuarterSections) * partitionQSCount*2)
  
  
    useAsHQS <- sample(partitionHQS, size=partitionHQSCount)
    useAsQQS <- partitionHQS[!(partitionHQS %in% useAsHQS)]
  
    useAsHQS <- hq2qq(useAsHQS, 4 * weSection)
    useAsQQS <- hq2qq(useAsQQS, 4 * weSection)
    refAsHQS <- rep(useAsHQS[1:(length(useAsHQS)/2)],times=2)
    refAsQQS <- useAsQQS
  
    # map to return
    myMap <- c()
  
    if( length(useAsQS) > 0 ) {
      myMap <-  cbind(refAsQS, useAsQS) 
    }
    if( length(useAsHQS) > 0 ) {
      myMap <-  rbind( myMap ,cbind(refAsHQS, useAsHQS) )
    }
    if( length(useAsQQS) > 0 ) {
      myMap <-  rbind( myMap ,cbind(refAsQQS, useAsQQS) )
    }

  } else {
    myMap <- c()
  
    if( length(useAsQS) > 0 ) {
      myMap <-  cbind(refAsQS, useAsQS) 
    }

  }

  colnames(myMap) <- c('object','qs')

  return( list( map=myMap, nsSection=nsSection, weSection=weSection) ) 

}


# this is a function that gets the spatial neighbors
# for input this takes an object of type simCrop
simCrop.getNeighbors <- function( x ) { 

  we  <- x$weSection
  ns  <- x$nsSection
  qqs <- x$map[,'qs']

  # 1. since the qqs relationships are known for a specified QS we start with getting those

  y <- matrix(0,nrow=length(qqs),ncol=4)

  y.index <- ( qqs %% ( we * 4 ) ) != 1
  y[y.index,1] <- qqs[y.index] -1

  y.index <- ( qqs %% ( we * 4 ) ) != 0
  y[y.index,2] <- qqs[y.index] +1

  y.index <- ( qqs /  ( we * 4 ) ) > 1
  y[y.index,3] <- qqs[y.index] - we * 4
  
  y.index <- ( qqs /  ( we * 4 ) ) <= (ns * 4 - 1) 
  y[y.index,4] <- qqs[y.index] + we * 4

  # 2. now we add our results on to myMap
  
 
  myMap <- cbind( x$map, y ) 
  
  colnames(myMap) <- c('object','qs','w.neighbor', 'e.neighbor','s.neighbor','n.neighbor')
 
  # add our results to the matrix 
  x$neighbors <- myMap

  return(x)
}


#This function creates a rook distance matrix from a neighbors matrix
simCrop.createRookDist <- function( x, fun=max ) {

  myNeighbors <- x$neighbors

  qsIndex <- sort(myNeighbors[,'qs'],index.return=T)$ix

  myObjects <- x$neighbors[,'object']
  myObjects.sorted <- myObjects[qsIndex]


  W.init <- apply(myNeighbors,1, function(x) {
        W.init <- matrix(0,nrow=1,ncol=nrow(myNeighbors))
        W.init[ c(x[3:6][x[3:6] > 0 ])] <- 1
        return(W.init)
    }
  )

  W <- aggregate( t(W.init), list( myObjects), fun )
  rownames(W) <- W[,1] 
  W <- W[,-1]
  W <- aggregate( t(W), list( myObjects.sorted), fun )
  rownames(W) <- W[,1] 
  W <- as.matrix(W[,-1])

  diag(W) <- 0

  return( W )
}


#This function gets the indexes for an object 
# It is assumed that qtr sections are listed going from nw to se e.g ( for nsSection = 1, weSection =2)
# 1 2 3 4 5 6 7 8
# 9 10 11 12 13 14 15 16
# 17 18 19 20 21 22 23 24
# 25 26 27 28 29 30 31 32
# l is the length of the side fo a square acre in pixels, default is 7
simCrop.getPixelsParcel <- function( a , objectId, l=7) {

  a.qs <- a$map[a$map[,'object'] == objectId,'qs'] 
  a.base <- a$weSection*4;

  y <- ceiling(a.qs/a.base)
  x <- ((a.qs-1) %% a.base) +1 
  
  lOffset <- l * a.base 

  start = (y - 1)  * lOffset * l + (x-1) * l + 1
  stop = start + l - 1

  b <- apply( cbind( start, stop), 1, function(x){ x[1]:x[2] } ) 
  return( cbind( objectId, c(apply( b, c(1,2), function(x) lOffset*(0:(l-1)) + x) ) ))

}


# add on pixel indexes to the parcel object as $value
simCrop.pixelParcel <- function(a, l=7) {

  b <- lapply( unique(a$map[,'object']),function(x) simCrop.getPixelsParcel(a,x,l) )
  
  x <- c()

  for( i in 1:length(b) ) {
    x <- rbind(x,b[[i]]) 
  }

  x <- x[sort(x[,2],index.return=T)$ix,]

  a$value <- cbind(x,0,0)
  a$l <- l

  colnames(a$value) <- c('object','pixel','value','error')
  
  return(a)
}


# convert the pixelParcel object to a raster
simCrop.pixelParcelToRaster <- function(a,x) {
  r <- raster(nrows= a$l * a$nsSection * 4, ncols=a$l * a$weSection * 4)
  values(r) <- a$value[,x]

  return(r)
}


# if transitionMatrix = T then use the transition matrix, else
# use the obs x length(transitionMatrix) matrix K
simCrop.addRotation <- function(a,K,transitionMatrix=F) {
  
  myObjects <- unique(a$map[,'object']) 
  myObjects.n <- length(myObjects)

  if(transitionMatrix==T) {
    K <- matrix( rep(c(t(K)),times=myObjects.n),byrow=T,nrow=myObjects.n) 
  }

  a$K <- cbind(myObjects,K) 

  return(a)
}

simCrop.plotCropTypes <- function(a) {
  
  cropType <- a$cropType

  objectIndex <- sort(cropType[,'myObjects'], index.return=T)$ix

  crops <- cropType[objectIndex,-1]

  m <- sqrt(nrow(crops))

  for( i in 1:(ncol(crops) )) {
    x11()
    plot( raster( matrix(crops[,i],byrow=T,nrow=m))) 
  }


}


# program that generates land cover type for objects in a simCrop object (this will be used by updates and pixel generation)  
# this land cover generation will be appended to any prior generation
# if K is provided, and a prior value exists the update will  based on the markov chain with transition matrix (K)
simCrop.generateCropTypes <- function(a,p,K) {
  
  # get the unique objects
  myObjects <- unique(a$map[,'object']) 

  # K will not be provided for a new update
  # and even if it is no update can be done without an initial value 
  if( !missing(K) && !is.null(a$cropType) ) {
  
    x <- a$cropType

    lastCol <- ncol(x)

    # if not init update with the conditional distribution
    p <- K[x[,lastCol],]

    # for the case of a single object
    if( is.null( nrow(p) ) ) {
      p <- matrix(p,nrow=1)
    }

    cropType.n <- ncol(p)

    a$cropType <- cbind( x, apply(p,1,function(x){ sample(1:cropType.n,size=1,prob=x) } )  )
    a$crops <-1:cropType.n

  
  } else if (!missing(p) ) {
    # if it's new or no K is provided we just use the provided p 

    cropType.n <- length(p)
    myObjects.n <- length(myObjects) 
     
    # select a primary crop
    x <- sapply(1:myObjects.n, function(x) { sample(1:cropType.n,size=1,prob=p) } )

    # append or add new
    if( is.null( a$cropType) ) {
      a$cropType <- cbind(myObjects, x) 
    } else {
      a$cropType <- cbind( a$cropType ,x) 
    }

    a$crops <-1:cropType.n

  } else if( !is.null(a$K) && !is.null(a$cropType) ) {
  
    x <- a$cropType
    K <- a$K[,-1]
    cropType.n <- sqrt(ncol(K))  
    lastCol <- ncol(x)

    # if not init update with the conditional distribution
    p <- t( apply( cbind( 1:nrow(x), x[,lastCol]) ,1, function(x) K[x[1], ((x[2]-1)*cropType.n + 1):(x[2]*cropType.n)] ) )

    a$cropType <- cbind( x, apply(p,1,function(x){ sample(1:cropType.n,size=1,prob=x) } )  )

    a$crops <-1:cropType.n
  
  } else {
    print("Error: no p or K provided, not doing anything")
  }

  return(a) 
}


### DEPRECIATED ###
# program that generates land cover type for a parcel
generateCropTypes <- function(x,p,K) {

  n <- length(x)

  cropTypes <- length(p)

  # if not init update with the conditional distribution
  if (x[1] != 0) p <- K[x[1],]

  # select a primary crop
  cropType <- sample(1:cropTypes,size=1,prob=p)

  y <- rep(cropType,times=n)

  return(y) 
}


# function that generates errors within a parcel
generateErrors <- function(x,p,errorRate) {

  y <- rbinom(length(x),size=1,prob=1-errorRate)   

  cropType <- x[1]

  errorCropType <- sample( (1:length(p))[-1 * cropType], size=1, prob=p[-1 * cropType]) 

  y[ y == 1] <- cropType 
  y[ y == 0] <- errorCropType

  return(y)
}


# update a set of pixel partials based on a function v 
simCrop.pixelParcelUpdate <- function(b,p,K,errorRate) {

  bValue <- b$value

  bValue <- bValue[sort( bValue[,'object'], index.return=T )$ix,]

  #bValue[,'value'] <- unlist(lapply(table(bValue[,'object']), function(x) generateCropTypes(x,p) ))

  # first need to subset by ownership
  nOwn <- unique(bValue[,'own'])

  for( i in 1:length(nOwn)) {
    # get the subset for a particular ownership
    subSetOwn <- which( bValue[,'own'] == nOwn[i] ) 
    #calculate next iteration and write it to the data set
    bValue[subSetOwn,'value'] <- unlist(aggregate( bValue[subSetOwn,'value'], list(bValue[subSetOwn,'object']), function(x) generateCropTypes(x,p,K[[i]]))$x)
  }

  bValue[,'error'] <- unlist(aggregate( bValue[,'value'], list(bValue[,'object']), function(x) generateErrors(x,p,errorRate))$x)

  b$value <- bValue[sort( bValue[,'pixel'],index.return=T)$ix,]

  return(b)

}


simCrop.pixelParcelChangeError <- function(b,p,errorRate) {
  bValue <- b$value

  bValue <- bValue[sort( bValue[,'object'], index.return=T )$ix,]
  bValue[,'error'] <- unlist(aggregate( bValue[,'value'], list(bValue[,'object']), function(x) generateErrors(x,p,errorRate))$x)
  b$value <- bValue[sort( bValue[,'pixel'],index.return=T)$ix,]

  return(b)
}

# function to add owner 
simCrop.objectOwner <- function(b,p.own) {

  # get a copy of the b$value
  bValue <- b$value

  # get unique objects
  bValueObjectUnique <- unique(bValue[,'object'])

  own <- rmultinom(length(bValueObjectUnique),size=1,prob=p.own) 

  # get the names
  bNames <- colnames(bValue)
 
  # assign ownership 
  bValue <- cbind(bValue,
    apply( data.frame(bValue[,'object']), 1, function(v) own[which(bValueObjectUnique %in% v)] ) ) 

  colnames(bValue) <- c(bNames,'own')

  b$value <- bValue

  return(b)
}


