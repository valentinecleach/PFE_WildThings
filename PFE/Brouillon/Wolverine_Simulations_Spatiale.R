library(nimble) 
library(raster) 

home <- getwd() 
setwd("~/work/PFE_WildThings/ScriptAndDataMCMC/Wolverine")
load("22.J_Fa1.RData") # Loads not only data but also other stuff
#load("22.J_Ma1.RData")
source("dbin_LESS_Cached_MultipleCovResponse.R")
source("pointProcess.R")

# Importations

library(necountries)
library(sf); library(dplyr)
library(ggplot2)

###### Spatial Simulation #####

coord_test <- rbinomPPSingle(n=1,
               lowerCoords = nimData$lowerHabCoords,
               upperCoords = nimData$upperHabCoords)

# Charger les pays
zone <- countries(c("Sweden", "Norway")) %>%
  st_as_sf()
zone_buff <- zone %>%
  st_buffer(dist=2.5)

hab_windows <- data.frame(
  x_low = nimData$lowerHabCoords[, "x"],
  y_low = nimData$lowerHabCoords[, "y"],
  x_up  = nimData$upperHabCoords[, "x"],
  y_up  = nimData$upperHabCoords[, "y"]
)

hab_centroids <- hab_windows %>%
  mutate(
    x = (x_low + x_up) / 2,
    y = (y_low + y_up) / 2
  ) %>%
  st_as_sf(coords = c("x", "y"), crs = st_crs(zone))

keep <- st_intersects(hab_centroids, zone_buff, sparse = FALSE)[,1]
hab_windows_fix <- hab_windows[keep, ]

LowerHabCoords_fix <- as.matrix(
  hab_windows_fix[, c("x_low", "y_low")]
)
UpperHabCoords_fix <- as.matrix(
  hab_windows_fix[, c("x_up", "y_up")]
)

# On prends des pts qui tombent dans notre zone
rbinomPPSingle_inside_zone <- function(
    lower, 
    upper, 
    polygon, 
    max_iter = 1000
    ){
  
  for (i in 1:max_iter) {
    # On simule un point
    pt <- rbinomPPSingle( 
      n=1,
      lowerCoords = lower,
      upperCoords = upper
      )
    # On regarde s'il tombe dans la zone
    pt_sf <- st_as_sf(
      data.frame(x=pt[1], y=pt[2]),
      coords = c("x","y"),
      crs = st_crs(polygon)
      )
    if (st_within(pt_sf, polygon, sparse = FALSE)[1,1]) {
      return(pt)
    }
  }
  stop("Could not sample inside polygon")
}

# Meme chose
rbinomMNormSourcePPSingle_inside_zone <- function(
    lower, 
    upper, 
    polygon, 
    sourceCoords,
    max_iter = 1000) {
  
  for (i in 1:max_iter) {
    pt <- rbinomMNormSourcePPSingle(
      n=1,
      lowerCoords = lower,
      upperCoords = upper,
      sourceCoords = sourceCoords,
      normSD = 1)
    pt_sf <- st_as_sf(data.frame(x=pt[1], y=pt[2]),
                      coords = c("x","y"),
                      crs = st_crs(polygon))
    if (st_within(pt_sf, polygon, sparse = FALSE)[1,1]) {
      return(pt)
    }
  }
  stop("Could not sample inside polygon")
}

zone_union <- st_union(zone)  # merge Sweden + Norway

sxy <- array(NA,
      dim = c(nimConstants$n.individuals,
              nimConstants$n.years,
              2))
# sxy dim : [i , t , coord]

sxy <- array(NA,
             dim = c(nimConstants$n.individuals,
                     nimConstants$n.years,
                     2))

for (j in 1:nimConstants$n.individuals){
  
  # Year 1: initial activity center
  sxy[j, 1, ] <- rbinomPPSingle_inside_zone(
    lower = LowerHabCoords_fix,
    upper = UpperHabCoords_fix,
    polygon = zone_union)
  
  # Movement through years
  for (i in 2:nimConstants$n.years){
    
    sxy[j, i, ] <- rbinomMNormSourcePPSingle_inside_zone(
      lower = LowerHabCoords_fix,
      upper = UpperHabCoords_fix,
      sourceCoords = sxy[j, i-1, ],
      polygon = zone_union)
  }
}

zone %>%
  ggplot()+
  geom_sf()+
  geom_point(aes(x=sxy[1,1,1],
                 y=sxy[1,1,2],# [i , t , coord]
                 col="year1"))+
  geom_point(aes(x=sxy[1,2,1],
                 y=sxy[1,2,2],
                 col="year2"))+
  geom_point(aes(x=sxy[1,3,1],
                 y=sxy[1,3,2],
                 col="year3"))+
  geom_point(aes(x=sxy[1,4,1],
                 y=sxy[1,4,2],
                 col="year4"))+
  geom_point(aes(x=sxy[1,5,1],
                 y=sxy[1,5,2],
                 col="year5"))+
  geom_point(aes(x=sxy[1,6,1],
                 y=sxy[1,6,2],
                 col="year6"))+
  geom_point(aes(x=sxy[1,7,1],
                 y=sxy[1,7,2],
                 col="year7"))+
  labs(title = "Individu1")+
  theme_bw()
