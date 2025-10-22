rm(list=ls())

library(nimble) # Bayesien
library(raster) # pour lire les maps pixels

home <- getwd() 

##########
## BEAR ##
##########
## LOAD NIMBLE FUNCTIONS AND DATA 



setwd(file.path(home,"ScriptAndDataMCMC/Bear")) # Selectionnement de l'ours


## LOAD FEMALE DATA
load("17.F_12_18_INPUTChain1.RData") # Importer dans env les choses


print(modelCode)

## LOAD MALE DATA 
#load("17.M_12_18_INPUTChain1.RData")

## LOAD THE CUSTOM NIMBLE FUNCTIONS
source("dbin_LESSCachedAllSparseBear_v2.R") # Import de fonctions
source("pointProcess.R")

#RUN NIMBLE MODEL (demonstration with single chain and low number of iterations)
model <- nimbleModel( code = modelCode
                      , constants = nimConstants
                      , data = nimData
                      , inits = nimInits
                      , check = FALSE       
                      , calculate = FALSE)  
cmodel <- compileNimble(model)
cmodel$calculate()
MCMCconf <- configureMCMC(model = model, monitors = c(nimParams),
                          control = list(reflective = TRUE, adaptScaleOnly = TRUE),
                          useConjugacy = FALSE) 
MCMC <- buildMCMC(MCMCconf) ## NE SE LANCE PAS??
cMCMC <- compileNimble(MCMC, project = model, resetFunctions = TRUE)
Runtime <- system.time(myNimbleOutput <- runMCMC( mcmc = cMCMC
                                                  , nburnin = 0
                                                  , niter = 100
                                                  , nchains = 1
                                                  , samplesAsCodaMCMC = TRUE))

#--PLOT STORED ANNUAL DENSITY RASTERS (average posterior utilization density, see Methods) 

setwd(file.path(home,"DensityRasterMaps"))
load("DensityRasterBrickBear.RData")
plot(DensityRasterBrick)



##########
## WOLF ##
##########


setwd(file.path(home,"ScriptAndDataMCMC/Wolf"))


## LOAD FEMALE DATA
load("9.F1218Cached_INPUTChain1.RData")
## LOAD MALE DATA 
#load("9.M1218Cached_INPUTChain1.RData")

## LOAD THE CUSTOM NIMBLE FUNCTIONS
source("dbin_LESSCachedAllSparseWolf.R")
source("pointProcess.R")

#RUN NIMBLE MODEL (demonstration with single chain and low number of iterations)
model <- nimbleModel( code = modelCode
                      , constants = nimConstants
                      , data = nimData
                      , inits = nimInits
                      , check = FALSE       
                      , calculate = FALSE)  
cmodel <- compileNimble(model)
cmodel$calculate()
MCMCconf <- configureMCMC(model = model, monitors = c(nimParams),
                          control = list(reflective = TRUE, adaptScaleOnly = TRUE),
                          useConjugacy = FALSE) # ca ca marche oks 
MCMC <- buildMCMC(MCMCconf)
cMCMC <- compileNimble(MCMC, project = model, resetFunctions = TRUE)
Runtime <- system.time(myNimbleOutput <- runMCMC( mcmc = cMCMC
                                                  , nburnin = 0
                                                  , niter = 100
                                                  , nchains = 1
                                                  , samplesAsCodaMCMC = TRUE))
plot(myNimbleOutput)
#--PLOT STORED ANNUAL DENSITY RASTERS (average posterior utilization density, see Methods) 

setwd(file.path(home,"DensityRasterMaps"))
load("DensityRasterBrickWolf.RData")
plot(DensityRasterBrick)




###############
## WOLVERINE ##
###############

# To help load without writing whole path (doesn't always work)
setwd(file.path(home,"ScriptAndDataMCMC/Wolverine")) 

## LOAD FEMALE DATA
load("22.J_Fa1.RData") # Loads not only data but also other stuff

## LOAD MALE DATA 
#load("22.J_Ma1.RData")

## LOAD THE CUSTOM NIMBLE FUNCTIONS
source("dbin_LESS_Cached_MultipleCovResponse.R")
source("pointProcess.R")

#RUN NIMBLE MODEL (demonstration with single chain and low number of iterations)


model <- nimble::nimbleModel(
  code = modelCode, 
  constants = nimConstants, 
  data = nimData, 
  inits = nimInits, 
  check = FALSE, 
  calculate = FALSE
  )  

cmodel <- nimble::compileNimble(model)

cmodel$calculate()

MCMCconf <- nimble::configureMCMC(
  model = model, 
  monitors = c(nimParams),
  control = list(reflective = TRUE, 
                 adaptScaleOnly = TRUE), 
  useConjugacy = FALSE) 


MCMC <- nimble::buildMCMC(MCMCconf) # okay

cMCMC <- nimble::compileNimble(MCMC, 
                               project = model, 
                               resetFunctions = TRUE)

Runtime <- system.time(myNimbleOutput <- nimble::runMCMC( 
  mcmc = cMCMC, 
  nburnin = 0, 
  niter = 100, 
  nchains = 1, 
  samplesAsCodaMCMC = TRUE)
  )

#--PLOT STORED ANNUAL DENSITY RASTERS (average posterior utilization density, see Methods) 

setwd(file.path(home,"DensityRasterMaps"))
load("DensityRasterBrickWolverine.RData")
plot(DensityRasterBrick)