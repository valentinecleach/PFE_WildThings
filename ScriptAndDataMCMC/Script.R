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

print(modelCode)

#RUN NIMBLE MODEL (demonstration with single chain and low number of iterations)

model <- nimble::nimbleModel( # Processes BUGS model code and returns NIMBLE model
  code = modelCode, # Code detailled in Wolverine_Explanations.Rmd
  constants = nimConstants, # n.individuals, n.detected, n.years, etc..
  data = nimData, # detection, non detection on diff years + data on habitat etc 
  inits = nimInits, # diff omega, gamma, sigma etc.
  check = FALSE, 
  calculate = FALSE
  )  

nimParams  

summary(model$getNodeNames())  # All variables/parameters in model
# 96498 parameters in model

model$getVarNames()   # Just the variable names

cmodel <- nimble::compileNimble(model) # Compile the model we created above

cmodel$calculate() # -1193130 -> exp(cmodel$calculate) is very small

MCMCconf <- nimble::configureMCMC( # Creates a default MCMC config for model
  model = model, 
  monitors = c(nimParams),
  control = list(reflective = TRUE, 
                 adaptScaleOnly = TRUE), 
  useConjugacy = FALSE) 


MCMC <- nimble::buildMCMC(MCMCconf) # Creates an uncompiled executable MCMC object

cMCMC <- nimble::compileNimble(MCMC, # algorithm that samples from model
                               project = model, # the statistical model (data, priors, likelihood)
                               resetFunctions = TRUE # 
                               )
cMCMC$my_initializeModel

Runtime <- system.time(
  myNimbleOutput <- nimble::runMCMC( # Runs MCMC
    mcmc = cMCMC, # compiled nimble model
    nburnin = 0, # no burnin period (to reduce time)
    niter = 100, # 100 iterations (here low to reduce time)
    nchains = 1, # Only 1 chaine (to reduce time here too)
    samplesAsCodaMCMC = TRUE # coda mcmc returned instead of R matrix of samples
    )
  )

#--PLOT STORED ANNUAL DENSITY RASTERS (average posterior utilization density, see Methods) 

setwd(file.path(home,"DensityRasterMaps"))
load("DensityRasterBrickWolverine.RData")
plot(DensityRasterBrick)