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
data("habitatRasters", envir = environment()) #

##-- Select a county to simulate data
##-- Here, we use Norrbotten county in Sweden
habitat$studyArea <- REGIONS %>%
  filter(county == "Norrbotten") %>%
  sf::st_union() %>%
  sf::st_as_sf()

##-- Disaggregate predefined habitat raster to the desired resolution
habRaster <- raster::disaggregate(
  x = habitatRasters[["Habitat"]],
  fact = raster::res(habitatRasters[["Habitat"]])/habitat$resolution)

##-- Make habitat from predefined Scandinavian raster of suitable habitat
habitat <- makeHabitatFromRaster(
  poly = habitat$studyArea,
  habitat.r = habitatRasters[["Habitat"]],
  buffer = habitat$buffer,
  plot.check = FALSE) %>%
  append(habitat,.)

##-- Retrieve number of habitat windows 
isHab <- habitat$habitat.r[] == 1
habitat$n.habwindows <- sum(isHab)

##-- Get coordinates of habitat windows
habitat$habitat.df <- cbind.data.frame(
  "id" = 1:habitat$n.habwindows,
  "x" = raster::coordinates(habitat$habitat.r)[isHab,1],
  "y" = raster::coordinates(habitat$habitat.r)[isHab,2])

##-- Make a spatial grid from polygon
habitat$grid <- raster::rasterToPolygons(
  habitat$habitat.r,
  fun = function(x){x>0}) %>%
  sf::st_as_sf() %>%
  mutate(id = 1:nrow(.))

##-- Generate a random covariate
habitat$covariate <- rnorm(habitat$n.habwindows)



## ------     3. SET-UP DETECTORS ------

##-- Generate raster of sub-detectors based on the study area
detectors$subdetectors.r <- raster::disaggregate(
  x = habitat$habitat.rWthBuffer,
  fact = raster::res(habitat$habitat.r)[1]/detectors$resolution.sub)

##-- Generate NGS detectors based on the raster of sub-detectors
detectors <- makeSearchGrid( 
  data = detectors$subdetectors.r,
  resolution = detectors$detResolution,
  div = (detectors$resolution/detectors$resolution.sub)^2,
  plot = FALSE) %>%
  append(detectors,.)

##-- Extract numbers of detectors
detectors$n.detectors <- nrow(detectors$main.detector.sp)

##-- Format detector locations & number of trials per detector
detectors$detectors.df <- cbind.data.frame(
  "id" = 1:detectors$n.detectors,
  "x" = sf::st_coordinates(detectors$main.detector.sp)[ ,1],
  "y" = sf::st_coordinates(detectors$main.detector.sp)[ ,2],
  "size" = as.vector(table(detectors$detector.sp$main.cell.id)))

##-- Make a spatial grid from polygon
detectors$grid <- raster::rasterToPolygons(
  x = raster::aggregate( x = detectors$subdetectors.r,
                         fact = detectors$resolution/detectors$resolution.sub),
  fun = function(x){x>0}) %>%
  sf::st_as_sf() %>%
  mutate(id = 1:nrow(.))



## ------     4. RESCALE COORDINATES ------

##-- Rescale coordinates
scaledObjects <- scaleCoordsToHabitatGrid(
  coordsData = detectors$detectors.df[ ,c("x","y")],
  coordsHabitatGridCenter = habitat$habitat.df[ ,c("x","y")])

##-- Get lower and upper habitat cell coordinates
habitat$lowerHabCoords <- scaledObjects$coordsHabitatGridCenterScaled - 0.5
habitat$upperHabCoords <- scaledObjects$coordsHabitatGridCenterScaled + 0.5

##-- Get local objects
localDetectors <- getLocalObjects( 
  habitatMask = habitat$habitat.mx,
  coords = scaledObjects$coordsDataScaled,
  dmax = detectors$maxDist/habitat$resolution,
  resizeFactor = 1,
  plot.check = TRUE)



##------------------------------------------------------------------------------

## ------ II. SIMULATE OPSCR DATASET ------

## ------     1. DEFINE MODEL CODE ------

modelCode <- nimbleCode({
  
  ##------ SPATIAL PROCESS ------##  
  
  ## Standard deviation of movement distribution
  tau ~ dgamma(0.001, 0.001)
  
  
  ## Intensity of movement point process
  betaHab ~ dnorm(0,0.01)
  logHabIntensity[1:n.habwindows] <- betaHab * habCov[1:n.habwindows]
  habIntensity[1:n.habwindows] <- exp(logHabIntensity[1:n.habwindows])
  logSumIntensity <- log(sum(habIntensity[1:n.habwindows]))
  
  ## FIRST YEAR 
  for(i in 1:M){
    s[i, 1:2,1] ~ dbernppAC(
      lowerCoords = lowerHabCoords[1:n.habwindows,1:2],
      upperCoords = upperHabCoords[1:n.habwindows,1:2],
      logIntensities = logHabIntensity[1:n.habwindows],
      logSumIntensity = logSumIntensity,
      habitatGrid = habitatGrid[1:y.max,1:x.max],
      numGridRows = y.max,
      numGridCols = x.max)
    
    ## T > 1 
    for(t in 2:n.years){
      s[i,1:2,t] ~ dbernppACmovement_normal(
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
  psi ~ dunif(0,1)   
  
  ## Initial state prob. vector
  omeg1[1:3] <- c(1-psi,psi,0)                                                 
  
  for(t in 1:(n.years-1)){
    ## Survival prob.
    phi[t] ~ dunif(0,1)      
    
    ## Recruitment prob.
    gamma[t] ~ dunif(0,1)
    
    ## Transition matrix
    omega[1,1:3,t] <- c(1-gamma[t],gamma[t],0)
    omega[2,1:3,t] <- c(0,phi[t],1-phi[t])
    omega[3,1:3,t] <- c(0,0,1)
  }#t
  
  ## Individual states ; here, we use a basic demographic model:
  ## z = 1 ==> not yet alive (unborn)
  ## z = 2 ==> alive
  ## z = 3 ==> dead
  for(i in 1:M){
    z[i,1] ~ dcat(omeg1[1:3])
    for(t in 1:(n.years-1)){
      z[i,t+1] ~ dcat(omega[z[i,t],1:3,t])
    }#t
  }#i
  
  
  
  ##----- DETECTION PROCESS -----## 
  
  ## Scale parameter of the detection function
  ## Note that all spatial parameters are relative to the habitat resolution
  ## e.g. here, a scale parameter of 5 corresponds to 5*20 = 100km
  sigma ~ dunif(0,10)
  
  for(t in 1:n.years){
    
    ## Annual baseline detection prob.
    p0[t] ~ dunif(0,1)
    
    for (i in 1:M){
      
      ## Alive indicator
      isAlive[i,t] <- (z[i,t] == 2) 
      
      ## Individual detections
      y[i, 1:lengthYCombined,t] ~ dbinomLocal_normal(
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
    N[t] <- sum(isAlive[1:M,t])
  }#t
  
})



## ------     2. PREPARE LISTS OF DATA, CONSTANTS & INITIAL VALUES ------

##-- Constants
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


##-- Data
simData <- list( detCoords = scaledObjects$coordsDataScaled,
                 size = detectors$detectors.df$size,
                 habCov = habitat$covariate,
                 lowerHabCoords = habitat$lowerHabCoords,
                 upperHabCoords = habitat$upperHabCoords,
                 habitatGrid = localDetectors$habitatGrid)


##-- Initial Values
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
                         calculate = T)  

##-- Identify nodes in the model to simulate
nodesToSim <- simModel$getDependencies( c("tau","betaHab",
                                          "psi","gamma","phi",
                                          "p0","sigma"),
                                        self = FALSE,
                                        downstream = TRUE,
                                        returnScalarComponents = TRUE)

##-- Simulate those nodes 
simModel$simulate(nodesToSim, includeData = FALSE)

##-- Check the simulated population size
N <- apply(simModel$z,2,function(x)sum(x==2))
N



##------------------------------------------------------------------------------

## ------ III. FIT OPSCR MODEL ------

## ------     1. BUNDLE SIMULATED DATA ------

##-- Individual state matrix
true.z <- z.data <- z.inits <- simModel$z
whichDet <- simModel$y[ ,1, ] > 0

##-- Set state to NA in the data for individuals not detected
z.data[!whichDet] <- NA
z.inits[whichDet] <- NA

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
cmodel <- compileNimble(model)
cmodel$calculate()
MCMCconf <- configureMCMC(model = model,
                          monitors = c( "tau","betaHab",
                                        "psi","gamma","phi",
                                        "p0","sigma",
                                       "N"),
                          monitors2 = c("s","z"),
                          control = list(reflective = TRUE),
                          thin = 1,
                          thin2 = 10)
MCMC <- buildMCMC(MCMCconf)   
cMCMC <- compileNimble(MCMC, project = model, resetFunctions = TRUE) 

##-- RUN THE MCMC (can take a long time)
nimOutput <- runMCMC( mcmc = cMCMC,                 
                      nburnin = 0,
                      niter = 1000,
                      nchains = 2,
                      samplesAsCodaMCMC = TRUE)



## ------     3. OUTPUT PROCESSING ------

##-- Traceplots 

## Population size
chainsPlot( nimOutput$samples,
            var = c("N[1]","N[2]","N[3]","N[4]","N[5]"),
            line = N)

## Other parameters
chainsPlot( nimOutput$samples, 
            var = c("p0","sigma","psi","gamma","phi","tau"), 
            line = c(data$p0,data$sigma,data$psi,data$gamma,data$phi,data$tau))



##------------------------------------------------------------------------------
