################################################################################
#This function calculates daylength for deriving PET via the Hamon (1963) 
#equation as described in Lu et al. (2005).
#It is calculated for each day of the year and then averaged over months.
#NB: Latitude should be in radians (not degrees).
#NB: Non-real values are produces for w (sunset hour angle) at times of the year
#and at latitudes where the sun never rises or sets. However, it makes sense to
#use just the real portion of the value, as it sets w to 0 in the winter and pi
#in the summer. This translates into N of 0 and 24, respectively

################################################################################

calcNr <- function(j) {  #what follows is part of the function

#Calculate declination, a function of Julian day.  It is a single value for each day.
dc <- 0.409*sin(((2*pi/365)*j)-1.39)


#Calculate the sunset hour angle, a function of latitude and declination. Note that at v. high latitudes, this function can produce non-real values.  
calcw = function(y) {Re(acos(as.complex(-1*tan(dc)*tan(y))))}
w = calc(lat,calcw)
rm(calcw,dc)

#Calculate N from the sunset hour angle
calcn = function(x) {(24/pi)*x}
n = calc(w,calcn)
rm(calcn,w)

return(n)

} #this closes the function