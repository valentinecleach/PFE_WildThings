
####### Time Simulations ######


###### Spatial Simulation #####

coord_test <- rbinomPPSingle(n=1,
               lowerCoords = nimData$lowerHabCoords,
               upperCoords = nimData$upperHabCoords)

library(necountries)
library(sf); library(dplyr)
library(ggplot2)

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

individu_t1 <- rbinomPPSingle(n=1,
                             lowerCoords = LowerHabCoords_fix,
                             upperCoords = UpperHabCoords_fix)

individu_t2 <- rbinomMNormSourcePPSingle(n=1,
                                         lowerCoords = LowerHabCoords_fix,
                                         upperCoords = UpperHabCoords_fix,
                                         sourceCoords = individu_t1,
                                         normSD = 1)

zone %>%
  ggplot()+
  geom_sf()+
  geom_point(aes(x=individu_t1[1],
                 y=individu_t1[2],
                 col="temps1"))+
  geom_point(aes(x=individu_t2[1],
                 y=individu_t2[2],
                 col="temps2"))+
  theme_bw()