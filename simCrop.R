
set.seed(1)

library(raster)
library(rgdal)
library(clv)

source('~/src/simCrop/floodFill.R')
source('~/src/simCrop/smooth.R')

# simCrop.R

# Township - 23040
# Section - 640
# Half-Section - 320
# Quarter-Section - 160
# Half-Quarter Section - 80
# Quarter-Quarter Section - 40
# assume approximately 56mx56m pixels so that one quarter section contains 49 pixels

# cropscape projection
aeaString <-"+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=23 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs +towgs84=0,0,0"

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
simCrop.partitionPLSS <- function( nsSection=1, weSection=1, landRatioOfQuarterQuarterSections=.50 ,landRatioOfHalfQuarterSections=.25){
#nsSection=1
#weSection=1
#landRatioOfQuarterQuarterSections=.50 
#landRatioOfHalfQuarterSections=.25




  # we first start off with establishing the number of quarter sections, and the number
  # of these sections that need to be converted into smaller units

  qsCount <- nsSection * weSection * 4 

  partitionQSCount <- floor( (landRatioOfHalfQuarterSections + landRatioOfQuarterQuarterSections) * qsCount ) 
  partitionQS <- sample(1:qsCount, size=partitionQSCount)

  # create a ref map for QQS -> QS 
  useAsQS <- (1:qsCount)[!(1:qsCount %in% partitionQS)]
  useAsQS <- q2hq(useAsQS)
  useAsQS <- hq2qq(useAsQS, 4 * weSection)
  
  refAsQS <- rep(useAsQS[1:(length(useAsQS)/4)],times=4)


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

  myMap <- rbind( cbind(refAsQS, useAsQS), cbind(refAsHQS, useAsHQS), cbind(refAsQQS,useAsQQS) ) 
  colnames(myMap) <- c('object','qs')

  return( list( map=myMap, nsSection=nsSection, weSection=weSection) ) 

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


# program that generates land cover type for a parcel
generateCropTypes <- function(x,p,K) {

  n <- length(x)

  # if not init update with the conditional distribution
  if (x[1] != 0) p <- K[x[1],]

  # select a primary crop
  cropType <- sample(1:4,size=1,prob=p)

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


# Test


# initial probability of land use 
p <- c(.05, .53, .37, .05)

# p ownership
p.own <- c(.8, .2)

# Non-Ag, Corn , Soybeans, Other Ag

K1 <- matrix( c( .98, .005, .005, .01,
                .01, .05,  .93,  .01,
                .01, .97,  .02,  .01,
                .01, .05,  .03,  .92),byrow=T,ncol=4) 

K2 <- matrix( c( .98, .005, .005, .01,
                .01, .90,  .10,  .01,
                .01, .50,  .50,  .01,
                .01, .05,  .03,  .92),byrow=T,ncol=4) 

K <- list(K1,K2)




a <- simCrop.partitionPLSS(6,6)
b <- simCrop.pixelParcel(a)

# add on own to the value object in b
b <- simCrop.objectOwner(b,p.own)

# create a list of raster objects and simulated landscapes
for(i in 1:10){

  if( i == 1 ) {
    b.list       <- list( simCrop.pixelParcelUpdate(b,p,K,.01) )
    r.list       <- list( simCrop.pixelParcelToRaster(b.list[[i]],'value') )
    projection(r.list[[i]]) <- CRS(aeaString)
  } 
  else {
    b.list[[i]]        <- simCrop.pixelParcelUpdate(b.list[[i-1]],p,K,.01)
    r.list[[i]] <- simCrop.pixelParcelToRaster(b.list[[i]],'value')
    projection(r.list[[i]]) <- CRS(aeaString)
  } 

}



#generate errors

runSim <- function( errorVal ) {
  j <- 1 # set initial index for sim

  for( h in c(1, 2, 3, 10) ) {  # iterate over years
    for( s in c(3,7,11,15) ) {  # iterate over x by x matrices 

      print( sprintf("h = %d, s = %d", h, s))  
      for( i in 1:10 ){
        if( i == 1 ) {
          b.temp <- simCrop.pixelParcelChangeError(b.list[[i]],p,errorVal)
          r.error.list <- list( simCrop.pixelParcelToRaster(b.temp,'error') )
          projection(r.error.list[[i]]) <- CRS(aeaString)
        } 
        else {
          b.temp <- simCrop.pixelParcelChangeError(b.list[[i]],p,errorVal)
          r.error.list[[i]] <- simCrop.pixelParcelToRaster(b.temp,'error')
          projection(r.error.list[[i]]) <- CRS(aeaString)
        } 
      }
      
      
      
      
      
      # create raster objects
      r.error <- brick( unlist(r.error.list) )
      r <- brick( unlist(r.list) )
      
      
      
      # we now get a new object index since the prior index my violate our parcel definition 
      r.object <- floodFill.raster(r) 
      #print(sprintf("New Object No. %d from %d", max(r.object), length(unique(b.list[[1]]$value[,'object']))))
      
      
      # run the kernel smoother on the item with errors
      r.smooth <- r.error
      values(r.smooth) <- localMode(r.error,h,s)
      
      # get the mse
      r.mse <- mean((values(r.smooth) - values( r ))^2)
      
      # calculate rand and other stats
      r.smooth.object <- floodFill.raster(r.smooth)
      
      r.compare <- std.ext(as.vector(r.smooth.object), as.vector(r.object) )
      
      r.compare[5] <- r.mse
      r.compare[6] <- errorVal 
      r.compare[7] <- h 
      r.compare[8] <- s 
      
      names(r.compare) <- c( names(r.compare)[1:4], "mse", "error", "h", "s")


      if( j > 1) {
        r.compare.all[[j]] <- r.compare
      }
      else {
        r.compare.all <- list(r.compare)
      }

     j <- j + 1 
    }
  }

} 




