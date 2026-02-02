
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

# Convertir les coordonnées en sf, avec le bon CRS
LowerHabCoords_fix <- nimData$lowerHabCoords %>%
  as.data.frame() %>%
  st_as_sf(crs = st_crs(zone), coords = c("x", "y"))

# Intersection avec la zone
LowerHabCoords_fix <- LowerHabCoords_fix %>%
  st_intersection(zone_buff)

# Extraire x et y de la géométrie, puis supprimer la géométrie
LowerHabCoords_fix <- LowerHabCoords_fix %>%
  mutate(
    x = st_coordinates(.)[, 1],
    y = st_coordinates(.)[, 2]
  ) %>%
  select(x, y) %>%
  distinct() %>%
  st_drop_geometry()%>%
  as.matrix()

# Idem avec Upper

UpperHabCoords_fix <- nimData$upperHabCoords %>%
  as.data.frame() %>%
  st_as_sf(crs = st_crs(zone), coords = c("x", "y")) %>%
  st_intersection(zone_buff) %>%
  mutate(
    x = st_coordinates(.)[, 1],
    y = st_coordinates(.)[, 2]
  ) %>%
  select(x, y) %>%
  distinct() %>%
  st_drop_geometry()%>%
  as.matrix()
####

# test
zone %>%
  ggplot()+
  geom_sf()+
  geom_point(aes(x=LowerHabCoords_fix[1,1],
                 y=LowerHabCoords_fix[1,2]))+
  labs(x="x",y="y",title="Test")+theme_bw()
# fin test

head(nimData$upperHabCoords);head(nimData$lowerHabCoords) 

head(UpperHabCoords_fix); head(LowerHabCoords_fix)

coord_test <- rbinomPPSingle(n=1,
                             lowerCoords = LowerHabCoords_fix,
                             upperCoords = UpperHabCoords_fix)


zone %>%
  ggplot()+
  geom_sf()+
  geom_point(aes(coord_test[1], coord_test[2]))

rbinomMNormSourcePPSingle()


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

zone %>%
  ggplot()+
  geom_sf()+
  geom_point(aes(x=coord_test[1],
                 y=coord_test[2]))+
  theme_bw()

individu_t2 <- rbinomMNormSourcePPSingle(n=1,
                                         lowerCoords = LowerHabCoords_fix,
                                         upperCoords = UpperHabCoords_fix,
                                         sourceCoords = individu_t1,
                                         normSD = 1
                                        )

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

