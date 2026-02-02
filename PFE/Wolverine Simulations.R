
####### Time Simulations ######


###### Spatial Simulation #####

nimData$lowerHabCoords
coord_test <- rbinomPPSingle(n=1,
               lowerCoords = nimData$lowerHabCoords,
               upperCoords = nimData$upperHabCoords)

plot(coord_test[1], coord_test[2])

library(necountries)

zone <- countries(c("Sweden","Norway"))

library(sf); library(dplyr)

library()
zone <- zone %>%
  st_as_sf()
zone
plot()+plot(coord_test[1], coord_test[2])

coord_test

rbinomMNormSourcePPSingle()

rb