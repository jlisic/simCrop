
set.seed(1)

library(raster)
library(rgdal)
library(clv)

source('~/src/simCrop/floodFill.R')
source('~/src/simCrop/smooth.R')
source('~/src/simCrop/simCrop.R')


# Township - 23040
# Section - 640
# Half-Section - 320
# Quarter-Section - 160
# Half-Quarter Section - 80
# Quarter-Quarter Section - 40
# assume approximately 56mx56m pixels so that one quarter section contains 49 pixels

# cropscape projection
aeaString <-"+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=23 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs +towgs84=0,0,0"

#################### Simulation ########################### 


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

  return(r.compare)

} 




