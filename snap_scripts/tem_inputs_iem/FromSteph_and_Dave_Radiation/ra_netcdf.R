#n=
#      [,1] [,2]
# [1,]    1   31
# [2,]   32   59
# [3,]   60   90
# [4,]   91  120
# [5,]  121  151
# [6,]  152  181
# [7,]  182  212
# [8,]  213  243
# [9,]  244  273
#[10,]  274  304
#[11,]  305  334
#[12,]  335  365)

#Calculate the latitude in radians at each grid cell center
lat.rad = matrix(NA,nx,ny)
for (y in 1:ny) {
	lat.rad[,y] = lat[y]*(pi/180)
}
rm(y)
lat = as.vector(lat)
lon = as.vector(lon)

#Pre-allocate a matrix for RA
ra = array(NA,c(nx,ny,12))
for (m in 1:12) {  #month loop
	#Identify the days of the yearassociated with that month
    md <-n[m,1]:n[m,2] #n is a matrix identifying the first and last days of the year
    #Also pre-allocate a temporary array to hold daily values of Ra
    temp <- array(NA,c(length(md),nx,ny))
    for (t in 1:length(md))    {  #loop over the number of days in each month
       r <- matrix(NA,nx,ny) #Make a matrix for  each day
       #Run the function  for each day
  	   r <- calcRa(md[t],lat.rad,nx,ny)
 	     #Slot the daily matrix r into array temp and delete r
       temp[t,,] = r
 	     rm(r)
    }      #End looping over number of days per month
    #Average the daily Ra values into a single month value
    ra[,,m] <- apply(temp,c(2,3),mean, na.rm = TRUE)
    rm(temp)   #clean up the temporary matrix
} #End looping over months

png('JAN_Ra.png')
image.plot(lon,lat,ra[,,1])
world(ylim = c(-90,90),add=T,shift=T)
dev.off()