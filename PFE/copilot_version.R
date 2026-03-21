##--------------------------------------------------------------------------- --
##--------------------------------------------------------------------------- --
## Script name: RovQuant OPSCR simulation (FIXED FOR STABILITY)
##--------------------------------------------------------------------------- --
##--------------------------------------------------------------------------- --

## ------ CLEAR-UP ENVIRONMENT ------

rm(list = ls())
gc()

## ------ LIBRARIES ------

packages <- c("sf", "dplyr", "raster", "ggplot2", "nimble", "nimbleSCR", 
              "rovquantR", "basicMCMCplots", "spatstat", "patchwork")

for(pkg in packages) {
  if(!require(pkg, character.only = TRUE)) {
    install.packages(pkg)
    library(pkg, character.only = TRUE)
  }
}

##--------------------------------------------------------------------------- --
## ------ I. SET-UP ------
##--------------------------------------------------------------------------- --

## REDUCED PARAMETERS FOR STABILITY
data = list(
  M = 15,
  n.years = 2,          # This means you need gamma and phi of length n.years-1 = 1
  tau = 1.2,
  betaHab = -0.8,
  psi = 0.1,
  gamma = 0.05,         # Single value for 1 transition
  phi = 0.80,           # Single value for 1 transition
  p0 = c(0.03, 0.03),   # Two values for 2 years
  sigma = 0.4
)

habitat = list(
  resolution = 50000,    # INCREASED (fewer cells)
  buffer.size = 30000    # REDUCED
)

detectors = list(
  resolution.sub = 10000,  # INCREASED (fewer sub-detectors)
  resolution = 50000,      # INCREASED
  maxDist = 50000          # REDUCED
)

## ------ 2 - Habitat ------
data("REGIONS", envir = environment())
data("habitatRasters", envir = environment())

habitat$studyArea <- REGIONS %>%
  filter(county == "Norrbotten") %>%
  sf::st_union() %>%
  sf::st_as_sf()

##-- Disaggregate or aggregate predefined habitat raster to the desired resolution
original_res <- raster::res(habitatRasters[["Habitat"]])[1]

if(habitat$resolution > original_res) {
  # If target resolution is larger, aggregate (coarsen)
  fact <- habitat$resolution / original_res
  habRaster <- raster::aggregate(
    x = habitatRasters[["Habitat"]],
    fact = fact)
} else {
  # If target resolution is smaller, disaggregate (refine)
  fact <- original_res / habitat$resolution
  habRaster <- raster::disaggregate(
    x = habitatRasters[["Habitat"]],
    fact = fact)
}

habitat <- makeHabitatFromRaster(
  poly = habitat$studyArea,
  habitat.r = habitatRasters[["Habitat"]],
  buffer = habitat$buffer,
  plot.check = FALSE) %>%
  append(habitat, .)

isHab <- habitat$habitat.r[] == 1
habitat$n.habwindows <- sum(isHab)

habitat$habitat.df <- cbind.data.frame(
  "id" = 1:habitat$n.habwindows,
  "x" = raster::coordinates(habitat$habitat.r)[isHab,1],
  "y" = raster::coordinates(habitat$habitat.r)[isHab,2])

habitat$grid <- raster::rasterToPolygons(habitat$habitat.r, fun = function(x){x>0}) %>%
  sf::st_as_sf() %>%
  mutate(id = 1:nrow(.))

habitat$covariate <- rnorm(habitat$n.habwindows)

## ------ 3 - Detectors ------
detectors$subdetectors.r <- raster::disaggregate(
  x = habitat$habitat.rWthBuffer,
  fact = raster::res(habitat$habitat.r)[1]/detectors$resolution.sub)

detectors <- makeSearchGrid(
  data = detectors$subdetectors.r,
  resolution = detectors$resolution,  # FIXED: was detResolution
  div = (detectors$resolution/detectors$resolution.sub)^2,
  plot = FALSE) %>%
  append(detectors, .)

detectors$n.detectors <- nrow(detectors$main.detector.sp)

detectors$detectors.df <- cbind.data.frame(
  "id" = 1:detectors$n.detectors,
  "x" = sf::st_coordinates(detectors$main.detector.sp)[ ,1],
  "y" = sf::st_coordinates(detectors$main.detector.sp)[ ,2],
  "size" = as.vector(table(detectors$detector.sp$main.cell.id)))

detectors$grid <- raster::rasterToPolygons(
  x = raster::aggregate(x = detectors$subdetectors.r,
                        fact = detectors$resolution/detectors$resolution.sub),
  fun = function(x){x>0}) %>%
  sf::st_as_sf() %>%
  mutate(id = 1:nrow(.))

## ------ 4 - Rescale Coords ------
scaledObjects <- scaleCoordsToHabitatGrid(
  coordsData = detectors$detectors.df[ ,c("x","y")],
  coordsHabitatGridCenter = habitat$habitat.df[ ,c("x","y")])

habitat$lowerHabCoords <- scaledObjects$coordsHabitatGridCenterScaled - 0.5
habitat$upperHabCoords <- scaledObjects$coordsHabitatGridCenterScaled + 0.5

localDetectors <- getLocalObjects(
  habitatMask = habitat$habitat.mx,
  coords = scaledObjects$coordsDataScaled,
  dmax = detectors$maxDist/habitat$resolution * 2,  # INCREASED: multiply by 2
  resizeFactor = 1,
  plot.check = FALSE)

cat("numLocalIndicesMax:", localDetectors$numLocalIndicesMax, "\n")

# DIAGNOSTIC: Check the structure
cat("localIndices class:", class(localDetectors$localIndices), "\n")
cat("localIndices dimensions:", dim(localDetectors$localIndices), "\n")
cat("localIndices first few rows:\n")
print(head(localDetectors$localIndices))
cat("\nnumLocalIndicesMax:", localDetectors$numLocalIndicesMax, "\n")
cat("numLocalIndices length:", length(localDetectors$numLocalIndices), "\n")

# ENSURE localIndices is a proper 2D matrix
if(!is.matrix(localDetectors$localIndices)) {
  cat("WARNING: localIndices is not a matrix, attempting to convert...\n")
  localDetectors$localIndices <- as.matrix(localDetectors$localIndices)
}

# Verify dimensions match expectations
expected_rows <- habitat$n.habwindows
expected_cols <- localDetectors$numLocalIndicesMax

cat("\nExpected dimensions:", expected_rows, "x", expected_cols, "\n")
cat("Actual dimensions:", nrow(localDetectors$localIndices), "x", ncol(localDetectors$localIndices), "\n")

if(nrow(localDetectors$localIndices) != expected_rows) {
  stop("Row mismatch: localIndices has ", nrow(localDetectors$localIndices), 
       " rows but expected ", expected_rows)
}

if(ncol(localDetectors$localIndices) != expected_cols) {
  stop("Column mismatch: localIndices has ", ncol(localDetectors$localIndices), 
       " columns but expected ", expected_cols)
}
##--------------------------------------------------------------------------- --
## ------ II. NIMBLE MODEL CODE ------
##--------------------------------------------------------------------------- --

modelCode <- nimbleCode({
  tau ~ dgamma(0.001, 0.001)
  betaHab ~ dnorm(0,0.01)
  logHabIntensity[1:n.habwindows] <- betaHab * habCov[1:n.habwindows]
  habIntensity[1:n.habwindows] <- exp(logHabIntensity[1:n.habwindows])
  logSumIntensity <- log(sum(habIntensity[1:n.habwindows]))
  
  for(i in 1:M){
    s[i, 1:2,1] ~ dbernppAC(
      lowerCoords = lowerHabCoords[1:n.habwindows,1:2],
      upperCoords = upperHabCoords[1:n.habwindows,1:2],
      logIntensities = logHabIntensity[1:n.habwindows],
      logSumIntensity = logSumIntensity,
      habitatGrid = habitatGrid[1:y.max,1:x.max],
      numGridRows = y.max,
      numGridCols = x.max)
    
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
    }
  }
  
  psi ~ dunif(0,1)
  omeg1[1:3] <- c(1-psi, psi, 0)
  
  for(t in 1:(n.years-1)){
    phi[t] ~ dunif(0,1)
    gamma[t] ~ dunif(0,1)
    omega[1,1:3,t] <- c(1-gamma[t], gamma[t], 0)
    omega[2,1:3,t] <- c(0, phi[t], 1-phi[t])
    omega[3,1:3,t] <- c(0, 0, 1)
  }
  
  for(i in 1:M){
    z[i,1] ~ dcat(omeg1[1:3])
    for(t in 1:(n.years-1)){
      z[i,t+1] ~ dcat(omega[z[i,t], 1:3, t])
    }
  }
  
  sigma ~ dunif(0,10)
  
  for(t in 1:n.years){
    p0[t] ~ dunif(0,1)
    for (i in 1:M){
      isAlive[i,t] <- (z[i,t] == 2)
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
    }
    N[t] <- sum(isAlive[1:M,t])
  }
})

##--------------------------------------------------------------------------- --
## ------ III. CONSTANTS, DATA, INITIAL VALUES ------
##--------------------------------------------------------------------------- --

simConstants <- list(
  M = data$M,
  n.years = data$n.years,
  n.detectors = detectors$n.detectors,                     
  n.habwindows = habitat$n.habwindows,
  y.max = dim(localDetectors$habitatGrid)[1],
  x.max = dim(localDetectors$habitatGrid)[2],
  resizeFactor = localDetectors$resizeFactor,
  numLocalIndicesMax = localDetectors$numLocalIndicesMax,
  localDetsIndices = localDetectors$localIndices,  # Now guaranteed to be 2D matrix
  localDetsNum = localDetectors$numLocalIndices,
  lengthYCombined = localDetectors$numLocalIndicesMax * 2 + 1)

# Verify dimensions
cat("localDetsIndices dimensions:", dim(simConstants$localDetsIndices), "\n")
cat("localDetsNum length:", length(simConstants$localDetsNum), "\n")
cat("numLocalIndicesMax:", simConstants$numLocalIndicesMax, "\n")

# If localDetsIndices is not 2D, reshape it
if(is.null(dim(simConstants$localDetsIndices))) {
  stop("localDetsIndices is not a matrix! Check nimbleSCR::getLocalObjects() output.")
}

if(nrow(simConstants$localDetsIndices) != habitat$n.habwindows) {
  stop("localDetsIndices row count does not match n.habwindows!")
}

simData <- list(
  detCoords = scaledObjects$coordsDataScaled,
  size = detectors$detectors.df$size,
  habCov = habitat$covariate,
  lowerHabCoords = habitat$lowerHabCoords,
  upperHabCoords = habitat$upperHabCoords,
  habitatGrid = localDetectors$habitatGrid)

simInits <- list(
  tau = data$tau,
  betaHab = data$betaHab,
  psi = data$psi,
  gamma = data$gamma,   # Will be length 1
  phi = data$phi,       # Will be length 1
  sigma = data$sigma,
  p0 = data$p0)         # Will be length 2

simModel <- nimbleModel(code = modelCode,
                        constants = simConstants,
                        data = simData,
                        inits = simInits,
                        check = F,       
                        calculate = FALSE)  # DON'T calculate during init

true.z <- simModel$z
z.data <- simModel$z
z.inits <- simModel$z
whichDet <- simModel$y[ ,1, ] > 0
whichDet[is.na(whichDet)] <- FALSE

for(i in 1:data$M) {
  if(!any(whichDet[i, ])) {
    z.data[i, ] <- NA
  }
}

z.inits[is.na(z.inits[ , 1]), 1] <- 2

nimConstants <- simConstants
nimData <- list(y = simModel$y, z = z.data,
                detCoords = scaledObjects$coordsDataScaled,
                size = detectors$detectors.df$size,
                habCov = habitat$covariate,
                lowerHabCoords = habitat$lowerHabCoords,
                upperHabCoords = habitat$upperHabCoords,
                habitatGrid = localDetectors$habitatGrid)

nimInits <- list(z = z.inits, s = simModel$s, tau = data$tau, betaHab = data$betaHab,
                 psi = data$psi, gamma = data$gamma, phi = data$phi,
                 sigma = data$sigma, p0 = data$p0)

##--------------------------------------------------------------------------- --
## ------ IV. SIMULATE DATASET ------
##--------------------------------------------------------------------------- --

hab_poly <- raster::rasterToPolygons(habitat$habitat.rWthBuffer, 
                                     fun=function(x){x==1}, 
                                     dissolve=TRUE) %>%
  sf::st_as_sf()

hab_poly_valid <- sf::st_make_valid(hab_poly)
win <- spatstat.geom::as.owin(sf::st_geometry(hab_poly_valid))

generate_centers <- function(type, M, window, hab_raster) {
  if(type == "IPP") {
    r_mat <- as.matrix(hab_raster)
    dens_im <- spatstat.geom::im(mat = r_mat,
                                 xcol = seq(raster::xmin(hab_raster), raster::xmax(hab_raster), length.out = ncol(r_mat)),
                                 yrow = seq(raster::ymin(hab_raster), raster::ymax(hab_raster), length.out = nrow(r_mat)))
    dens_im <- dens_im / max(dens_im, na.rm = TRUE)
    pts <- spatstat.random::rpoint(n = M, f = dens_im, win = window)
  } else if(type == "Neyman-Scott") {
    pts_raw <- spatstat.random::rMatClust(kappa = 0.00005, scale = 10000, mu = 10, win = window)
    pts <- pts_raw[sample(1:pts_raw$n, M)]
  }
  
  coords <- data.frame(x = pts$x, y = pts$y)
  scaled <- nimbleSCR::scaleCoordsToHabitatGrid(coordsData = coords,
                                                coordsHabitatGridCenter = habitat$habitat.df[,c("x","y")])
  return(list(scaled = scaled$coordsDataScaled, raw = coords))
}

run_simulation <- function(scenario) {
  cat("\n=== Starting", scenario, "scenario ===\n")
  
  centers <- generate_centers(scenario, data$M, win, habitat$habitat.r)
  s_init <- as.matrix(centers$scaled)
  
  nimInits$s <- array(NA_real_, dim = c(data$M, 2, data$n.years))
  for(t in 1:data$n.years) nimInits$s[,,t] <- s_init
  
  model <- nimbleModel(code = modelCode, constants = nimConstants, data = nimData, inits = nimInits, check = F, calculate = FALSE)
  cmodel <- compileNimble(model)
  
  MCMCconf <- configureMCMC(model, monitors = c("N", "sigma", "p0"))
  MCMC <- buildMCMC(MCMCconf)
  cMCMC <- compileNimble(MCMC, project = model)
  
  # REDUCED ITERATIONS FOR STABILITY
  nimOutput <- runMCMC(cMCMC, niter = 50, nchains = 1, nburnin = 10)
  
  cat("=== Completed", scenario, "===\n")
  return(list(output = nimOutput, coords = centers$raw))
}

res_IPP <- run_simulation("IPP")
res_NS  <- run_simulation("Neyman-Scott")

# Filter points to keep only those inside the window
inside_IPP <- spatstat.geom::inside.owin(res_IPP$coords$x, res_IPP$coords$y, win)
coords_IPP <- res_IPP$coords[inside_IPP, ]

inside_NS <- spatstat.geom::inside.owin(res_NS$coords$x, res_NS$coords$y, win)
coords_NS <- res_IPP$coords[inside_NS, ]

##--------------------------------------------------------------------------- --
## ------ V. GRAPHICS ------
##--------------------------------------------------------------------------- --

p1 <- ggplot() +
  geom_sf(data = hab_poly, fill = "grey90", color = NA) +
  geom_point(data = coords_IPP, aes(x = x, y = y),
             color = "blue", size = 1.2) +
  ggtitle("IPP") +
  theme_minimal()

p2 <- ggplot() +
  geom_sf(data = hab_poly, fill = "grey90", color = NA) +
  geom_point(data = coords_NS, aes(x = x, y = y),
             color = "red", size = 1.2) +
  ggtitle("Neyman–Scott") +
  theme_minimal()

print(p1 + p2)

##---------------------------------------------------------------------------- -
##---------------------------------------------------------------------------- -

## ------ VI. ANALYSE STATIALE (fonction K Ripley) ------

pp_IPP <- ppp(coords_IPP$x, coords_IPP$x, window = win)
pp_NS  <- ppp(coords_NS$x, coords_NS$y, window = win)

K_IPP <- Kest(pp_IPP)
K_NS  <- Kest(pp_NS)

plot(K_IPP$r, K_IPP$iso, type="l", col="blue", lwd=2,
     xlab="Distance", ylab="K(r)",
     main="Fonction K : IPP vs Neyman–Scott")

lines(K_NS$r, K_NS$iso, col="red", lwd=2)

legend("topleft",
       legend=c("IPP (Poisson)", "Neyman–Scott (cluster)"),
       col=c("blue","red"), lwd=2)
