################################################################################
#This function calculates extraterrestrial solar radiation (radiation at the top
#of the earth's atmosphere, aka insolation) from day of the year and latitude.
#It is calculated separately for each day of the year.
#NB: Latitude should be in radians (not degrees).
#NB: Non-real values are produced for w (sunset hour angle) at times of the year
#and at latitudes where the sun never rises or sets. However, it makes sense to
#use just the real portion of the value, as it sets w to 0 in the winter and pi
#(~3.14) in the summer. 
#Source: Allen et al. (1998)
################################################################################

calcRa <- function(jd,lat,nrow,ncol) {  #what follows is part of the function

#Calculate the earth-sun distance, which is a function solely of Julian day.  It is a single value for each day
d <- 1+(0.033*cos(2*pi*jd/365))

#Calculate declination, a function of Julian day.  It is a single value for each day.
dc <- 0.409*sin(((2*pi/365)*jd)-1.39)

#Pre-allocate the sunset hour angle (w) matrix
w <- matrix(NA,nrow,ncol)
#Calculate the sunset hour angle, a function of latitude and declination. Note that at v. high latitudes, this function can produce non-real values.  
w[!is.na(lat)] <- as.matrix(Re(acos(as.complex(-1*tan(dc)*tan(lat[!is.na(lat)])))))

#Pre-allocate the Ra matrix
Ra <- matrix(NA,nrow,ncol)

#Calculate Ra using the above variables
#S is the solar constant and is =0.082 MJ/m2min
Ra[!is.na(lat)] <- as.matrix((24*60/pi)*d*0.082*(w[!is.na(lat)]*sin(lat[!is.na(lat)])*sin(dc)+cos(lat[!is.na(lat)])*cos(dc)*sin(w[!is.na(lat)])))
return(Ra)
} #this closes the function

