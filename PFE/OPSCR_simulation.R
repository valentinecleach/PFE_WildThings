##--------------------------------------------------------------------------- --
##
## Script name: RovQuant OPSCR simulation
##
## This R script performs:
## 1. the simulation of the data, habitat and detectors using the 'rovquantR' and 'nimbleSCR' packages
## 2. the data preparation for the OPSCR analysis with the 'nimbleSCR' package
## 3. the model fitting using 'nimble' and 'nimbleSCR'
## 4. the post-processing of the MCMC output
##
## Author: Pierre Dupont
## Email: pierre.dupont@nmbu.no
##
## Date Created: 16/09/2025
##
## Copyright (c) Applied Quantitative Ecology Group (AQEG), 2025
## Faculty of Environmental Sciences and Natural Resource Management (MINA)
## Norwegian University of Life Sciences (NMBU), Ås, Norway 
##
## Script to understand the use of OPSCR, not linked to a particular article.
##--------------------------------------------------------------------------- --
##
## Notes: 
## This is based on the development branch of 'rovquantR' v.0.1
##   
##--------------------------------------------------------------------------- --

## ------ CLEAR-UP ENVIRONMENT ------

rm(list = ls())
gc()


## ------ INSTALL 'rovquantR' FROM GITHUB ------

## Ctrl + Shift + F10 (to restart R session)
# devtools::install_github("PierreDupont/rovquantR@devel")


## ------ LOAD REQUIRED LIBRARIES ------

library(sf)
library(dplyr)
library(raster)
library(ggplot2)
library(nimble)
library(nimbleSCR)
library(rovquantR)
library(basicMCMCplots)

##------------------------------------------------------------------------------

## ------ I. SIMULATION SET-UP ------

## ------     1. SET SIMULATION PARAMETER VALUES ------

##-- Population specifications
data = list( M = 250, # nb max d'individus potentiel
             n.years = 5, # duree de la simulation
             tau = 1.2, # dispersion du mouvement spatial
             betaHab = -0.8, # habitat defavorable, habitat impact distrib indiv
             psi = 0.1, # proba initial indiv present
             gamma = c(0.05,0.06,0.06,0.07), # probas recrutement pour chaque annee
             phi = runif(4, 0.75, 0.85), # proba de survie anuelle (al entre 0.75-0.85)
             p0 = runif(5, 0.02,0.05), # proba de détection de base
             sigma = 0.4) # param spatial de la detection

##-- Habitat specifications
habitat = list( resolution = 20000, # taille de pixels
                buffer.size = 60000) # zone tampon

##-- Detectors specifications
detectors = list( resolution.sub = 2000, # resolution des sous-detecteurs (2km)
                  resolution = 10000, # resolution principale (10km)
                  maxDist = 80000) # distance max de detection (80km)



## ------     2. SET-UP HABITAT ------

##-- Load pre-defined habitat rasters and shapefiles of Scandinavia
data("REGIONS", envir = environment()) # objet spatial des regions scandinaves
data("habitatRasters", envir = environment()) # rasters d'habitat prédéfinis (??)

##-- Select a county to simulate data
##-- Here, we use Norrbotten county in Sweden
habitat$studyArea <- REGIONS %>% # dans les regions
  filter(county == "Norrbotten") %>% # Suède
  sf::st_union() %>% # fusion des polygones
  sf::st_as_sf() # objet spatial

##-- Disaggregate predefined habitat raster to the desired resolution
# raster = ensemble de pixels
habRaster <- raster::disaggregate( # augmente la resolution du raster
  x = habitatRasters[["Habitat"]], # raster d'habitat d'origine
  fact = raster::res(habitatRasters[["Habitat"]])/habitat$resolution)
# facteur de division pour nouvelle resolution

##-- Make habitat from predefined Scandinavian raster of suitable habitat
habitat <- makeHabitatFromRaster( # transforme un raster en habitat exploitable
  poly = habitat$studyArea, # recup zone d'etude
  habitat.r = habitatRasters[["Habitat"]], # recup raster d'habitat
  buffer = habitat$buffer, # recup zone tampon
  plot.check = FALSE) %>% # pas d'affichage de la figure
  append(habitat,.) # ajout l'objet a habitat

##-- Retrieve number of habitat windows 
isHab <- habitat$habitat.r[] == 1 # id pixels qui sont habitat
habitat$n.habwindows <- sum(isHab) # nb total de cellule habitat

##-- Get coordinates of habitat windows
habitat$habitat.df <- cbind.data.frame(
  "id" = 1:habitat$n.habwindows, # identifiant unique pour les cellules habitat
  "x" = raster::coordinates(habitat$habitat.r)[isHab,1], # coordonnes centres cellules
  "y" = raster::coordinates(habitat$habitat.r)[isHab,2])

##-- Make a spatial grid from polygon
habitat$grid <- raster::rasterToPolygons( # conversion raster en polygones spatiaux
  habitat$habitat.r,
  fun = function(x){x>0}) %>%
  sf::st_as_sf() %>%
  mutate(id = 1:nrow(.)) # ajout id a chaque cellule

##-- Generate a random covariate
habitat$covariate <- rnorm(habitat$n.habwindows) # genre covar environ al



## ------     3. SET-UP DETECTORS ------

##-- Generate raster of sub-detectors based on the study area
detectors$subdetectors.r <- raster::disaggregate( # raster de sous-detecteurs
  x = habitat$habitat.rWthBuffer, # raster d'habitat avec buffer autour de la zone
  fact = raster::res(habitat$habitat.r)[1]/detectors$resolution.sub) # subdivision

##-- Generate NGS detectors based on the raster of sub-detectors
detectors <- makeSearchGrid( # réseau de détecteurs spatiaux
  data = detectors$subdetectors.r, # raster de base utilisé pour placer les detecteurs
  resolution = detectors$detResolution, # distance entre detecteurs
  div = (detectors$resolution/detectors$resolution.sub)^2, # nb de sous-cell dans une cell principale
  plot = FALSE) %>% 
  append(detectors,.) # ajout à la liste detectors

##-- Extract numbers of detectors
detectors$n.detectors <- nrow(detectors$main.detector.sp) # nb total de detecteurs

##-- Format detector locations & number of trials per detector
detectors$detectors.df <- cbind.data.frame( # dataframe des detecteurs
  "id" = 1:detectors$n.detectors, # id unique pour chaque detect
  "x" = sf::st_coordinates(detectors$main.detector.sp)[ ,1], # coord x du detect
  "y" = sf::st_coordinates(detectors$main.detector.sp)[ ,2], # coord y
  "size" = as.vector(table(detectors$detector.sp$main.cell.id))) # nb de ss-detect dans chaque detect principal

##-- Make a spatial grid from polygon
detectors$grid <- raster::rasterToPolygons( # convertion raster en polygone
  x = raster::aggregate( x = detectors$subdetectors.r, # aggregation des ss-detect
                         fact = detectors$resolution/detectors$resolution.sub), # nb tot de ss-cell à regr
  fun = function(x){x>0}) %>% # garder les cellules actives
  sf::st_as_sf() %>% # modif en objet spatial
  mutate(id = 1:nrow(.)) # ajout id a chaque cell



## ------     4. RESCALE COORDINATES ------

##-- Rescale coordinates
scaledObjects <- scaleCoordsToHabitatGrid( # toutes les coord dans meme syst d'echelle
  coordsData = detectors$detectors.df[ ,c("x","y")], # coord detecteurs
  coordsHabitatGridCenter = habitat$habitat.df[ ,c("x","y")]) # corrdonnes centres cell habitat

##-- Get lower and upper habitat cell coordinates
habitat$lowerHabCoords <- scaledObjects$coordsHabitatGridCenterScaled - 0.5 # coord coin inf cell hab
habitat$upperHabCoords <- scaledObjects$coordsHabitatGridCenterScaled + 0.5 # coord coin sup cell hab

##-- Get local objects
localDetectors <- getLocalObjects( # optimiser calculs SCR (id detecteurs proches de chaque cell)
  habitatMask = habitat$habitat.mx, # matrice hab
  coords = scaledObjects$coordsDataScaled, # coord detect mises à l'echelle
  dmax = detectors$maxDist/habitat$resolution, # distance maximale de detection
  resizeFactor = 1, # facteur d'echelle
  plot.check = TRUE)



##------------------------------------------------------------------------------

## ------ II. SIMULATE OPSCR DATASET ------

## ------     1. DEFINE MODEL CODE ------

modelCode <- nimbleCode({ # modele probabiliste nimble
  
  ##------ SPATIAL PROCESS ------##  
  
  ## Standard deviation of movement distribution
  tau ~ dgamma(0.001, 0.001) # prior dispersion du mouvement (distrib gamma)
  
  ## Intensity of movement point process
  betaHab ~ dnorm(0,0.01) # prior effet de la covar habitat
  logHabIntensity[1:n.habwindows] <- betaHab * habCov[1:n.habwindows] # intensite spatiale (depend covar habitat)
  habIntensity[1:n.habwindows] <- exp(logHabIntensity[1:n.habwindows]) # transformation exponentielle
  logSumIntensity <- log(sum(habIntensity[1:n.habwindows])) # normalisation du processus spatial
  
  ## FIRST YEAR 
  for(i in 1:M){ # pour chaque indiv
    s[i, 1:2,1] ~ dbernppAC(
      lowerCoords = lowerHabCoords[1:n.habwindows,1:2],
      upperCoords = upperHabCoords[1:n.habwindows,1:2],
      logIntensities = logHabIntensity[1:n.habwindows],
      logSumIntensity = logSumIntensity,
      habitatGrid = habitatGrid[1:y.max,1:x.max],
      numGridRows = y.max,
      numGridCols = x.max) # localisation du CA la premiere annee
    
    ## T > 1 
    for(t in 2:n.years){ # boucle temporelle
      s[i,1:2,t] ~ dbernppACmovement_normal( # modele de mouvement spatial des indiv
        lowerCoords = lowerHabCoords[1:n.habwindows,1:2],
        upperCoords = upperHabCoords[1:n.habwindows,1:2],
        s = s[i,1:2,t-1],
        sd = tau,
        baseIntensities = habIntensity[1:n.habwindows],
        habitatGrid = habitatGrid[1:y.max,1:x.max],
        numGridRows = y.max,
        numGridCols = x.max,
        numWindows = n.habwindows)
    }#t
  }#i
  
  
  
  ##----- DEMOGRAPHIC PROCESS -----## 
  
  ## Initial inclusion prob.
  psi ~ dunif(0,1) # proba init d'inclusion (alive)
  
  ## Initial state prob. vector
  omeg1[1:3] <- c(1-psi,psi,0) # proba init des etats
  
  for(t in 1:(n.years-1)){ # chaque annee
    ## Survival prob.
    phi[t] ~ dunif(0,1) # proba de survie
    
    ## Recruitment prob.
    gamma[t] ~ dunif(0,1) # proba de recrutement
    
    ## Transition matrix
    omega[1,1:3,t] <- c(1-gamma[t],gamma[t],0)
    omega[2,1:3,t] <- c(0,phi[t],1-phi[t])
    omega[3,1:3,t] <- c(0,0,1)
  }#t
  
  ## Individual states ; here, we use a basic demographic model:
  ## z = 1 ==> not yet alive (unborn)
  ## z = 2 ==> alive
  ## z = 3 ==> dead
  for(i in 1:M){ # pour chaque indiv
    z[i,1] ~ dcat(omeg1[1:3]) # etat init
    for(t in 1:(n.years-1)){ # temps d'apres
      z[i,t+1] ~ dcat(omega[z[i,t],1:3,t]) # transition entre etats (modele de Markov)
    }#t
  }#i
  
  
  
  ##----- DETECTION PROCESS -----## 
  
  ## Scale parameter of the detection function
  ## Note that all spatial parameters are relative to the habitat resolution
  ## e.g. here, a scale parameter of 5 corresponds to 5*20 = 100km
  sigma ~ dunif(0,10) # param spatial de la detection
  
  for(t in 1:n.years){ # chaque annee
    
    ## Annual baseline detection prob.
    p0[t] ~ dunif(0,1) # proba de detection de base
    
    for (i in 1:M){ # chaque indiv
      
      ## Alive indicator
      isAlive[i,t] <- (z[i,t] == 2) # indicatrice si etat alive
      
      ## Individual detections
      y[i, 1:lengthYCombined,t] ~ dbinomLocal_normal(
        # modele de detection binomiale spatiale --> selon dist, detects proches, etat indiv
        size = size[1:n.detectors],
        p0 = p0[t],
        s = s[i,1:2,t],
        sigma = sigma,
        trapCoords = detCoords[1:n.detectors,1:2],
        localTrapsIndices = localDetsIndices[1:n.habwindows,1:numLocalIndicesMax],
        localTrapsNum = localDetsNum[1:n.habwindows],
        resizeFactor = resizeFactor,
        habitatGrid = habitatGrid[1:y.max,1:x.max],
        indicator = isAlive[i,t],
        lengthYCombined = lengthYCombined)
    }#i
    
    ## Population size
    N[t] <- sum(isAlive[1:M,t]) # taille annuelle de la population
  }#t
  
})



## ------     2. PREPARE LISTS OF DATA, CONSTANTS & INITIAL VALUES ------

##-- Constants = constantes du modele
simConstants <- list( M = data$M,
                      n.years = data$n.years,
                      n.detectors = detectors$n.detectors,                     
                      n.habwindows = habitat$n.habwindows,
                      y.max = dim(localDetectors$habitatGrid)[1],
                      x.max = dim(localDetectors$habitatGrid)[2],
                      resizeFactor = localDetectors$resizeFactor,
                      numLocalIndicesMax = localDetectors$numLocalIndicesMax,
                      localDetsIndices = localDetectors$localIndices,
                      localDetsNum = localDetectors$numLocalIndices,
                      lengthYCombined = localDetectors$numLocalIndicesMax * 2 + 1)


##-- Data = donnes observees (par simulation)
simData <- list( detCoords = scaledObjects$coordsDataScaled,
                 size = detectors$detectors.df$size,
                 habCov = habitat$covariate,
                 lowerHabCoords = habitat$lowerHabCoords,
                 upperHabCoords = habitat$upperHabCoords,
                 habitatGrid = localDetectors$habitatGrid)


##-- Initial Values = valeurs initiales des param
simInits <- list( tau = data$tau,
                  betaHab = data$betaHab,
                  psi = data$psi,
                  gamma = data$gamma,
                  phi = data$phi,
                  sigma = data$sigma,
                  p0 = data$p0) 



## ------     3. SIMULATE DETECTION HISTORY FROM THE NIMBLE MODEL OBJECT ------

##-- Build the nimble model object
simModel <- nimbleModel( code = modelCode,
                         constants = simConstants,
                         data = simData,
                         inits = simInits,
                         check = F,       
                         calculate = T) # creation modele NIMBLE

##-- Identify nodes in the model to simulate
nodesToSim <- simModel$getDependencies( c("tau","betaHab", # id variables a simuler
                                          "psi","gamma","phi",
                                          "p0","sigma"),
                                        self = FALSE,
                                        downstream = TRUE,
                                        returnScalarComponents = TRUE)

##-- Simulate those nodes 
simModel$simulate(nodesToSim, includeData = FALSE) # generation donnes simulees

##-- Check the simulated population size
N <- apply(simModel$z,2,function(x)sum(x==2)) # calcul de la vraie taille de pop
N



##------------------------------------------------------------------------------

## ------ III. FIT OPSCR MODEL ------

## ------     1. BUNDLE SIMULATED DATA ------

##-- Individual state matrix
true.z <- z.data <- z.inits <- simModel$z # etats reels des individus
whichDet <- simModel$y[ ,1, ] > 0 # id indiv detectes

##-- Set state to NA in the data for individuals not detected
z.data[!whichDet] <- NA # etat inconnus pour indiv non detectes
z.inits[whichDet] <- NA # pas d'initialisation si detectes

##-- Constants
nimConstants <- list( M = data$M,
                      n.years = data$n.years,
                      n.detectors = detectors$n.detectors,                     
                      n.habwindows = habitat$n.habwindows,
                      y.max = dim(localDetectors$habitatGrid)[1],
                      x.max = dim(localDetectors$habitatGrid)[2],
                      resizeFactor = localDetectors$resizeFactor,
                      numLocalIndicesMax = localDetectors$numLocalIndicesMax,
                      localDetsIndices = localDetectors$localIndices,
                      localDetsNum = localDetectors$numLocalIndices,
                      lengthYCombined = localDetectors$numLocalIndicesMax * 2 + 1)


##-- Data
nimData <- list( y = simModel$y,
                 z = z.data,
                 detCoords = scaledObjects$coordsDataScaled,
                 size = detectors$detectors.df$size,
                 habCov = habitat$covariate,
                 lowerHabCoords = habitat$lowerHabCoords,
                 upperHabCoords = habitat$upperHabCoords,
                 habitatGrid = localDetectors$habitatGrid)


##-- Initial Values
nimInits <- list( z = z.inits,
                  s = simModel$s,
                  tau = data$tau,
                  betaHab = data$betaHab,
                  psi = data$psi,
                  gamma = data$gamma,
                  phi = data$phi,
                  sigma = data$sigma,
                  p0 = data$p0) 



## ------     2. FIT THE MODEL ------

##-- CREATE & COMPILE THE NIMBLE MODEL
model <- nimbleModel( code = modelCode,
                      constants = nimConstants,
                      data = nimData,
                      inits = nimInits,
                      check = F,
                      calculate = F)
cmodel <- compileNimble(model) # compilation en c++
cmodel$calculate()
MCMCconf <- configureMCMC(model = model,
                          monitors = c( "tau","betaHab",
                                        "psi","gamma","phi",
                                        "p0","sigma",
                                       "N"),
                          monitors2 = c("s","z"),
                          control = list(reflective = TRUE),
                          thin = 1,
                          thin2 = 10) # config du MCMC --> param a suivre, thinning, options
MCMC <- buildMCMC(MCMCconf) # construction de l'algo MCMC
cMCMC <- compileNimble(MCMC, project = model, resetFunctions = TRUE) # compile

##-- RUN THE MCMC (can take a long time)
nimOutput <- runMCMC( mcmc = cMCMC, # execution MCMC
                      nburnin = 0,
                      niter = 1000, # 1000 iterations
                      nchains = 2, # 2 chaines
                      samplesAsCodaMCMC = TRUE)



## ------     3. OUTPUT PROCESSING ------

##-- Traceplots 

## Population size
chainsPlot( nimOutput$samples,
            var = c("N[1]","N[2]","N[3]","N[4]","N[5]"),
            line = N) # pour verifier visuellement la convergence

## Other parameters
chainsPlot( nimOutput$samples, 
            var = c("p0","sigma","psi","gamma","phi","tau"), 
            line = c(data$p0,data$sigma,data$psi,data$gamma,data$phi,data$tau))
# graphiques pour les parametres

##------------------------------------------------------------------------------
