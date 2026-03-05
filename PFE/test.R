library(raster)
library(sp)

# --- Paramètres de simulation ---
set.seed(123)
area_ext <- c(0, 100, 0, 100) # Étendue X et Y
n_years <- 3
phi <- 0.8                   # Survie annuelle
move_sd <- 2                 # Mouvement des centres d'activité (AC)
n_parents <- 10              # Nombre de clusters (meutes/zones)
lambda_offspring <- 15       # Nombre moyen d'individus par cluster
sigma_cluster <- 5           # Rayon du cluster (Neyman-Scott)

# --- Étape 1 : Création de la population initiale (Année 1) ---
parents_x <- runif(n_parents, area_ext[1], area_ext[2])
parents_y <- runif(n_parents, area_ext[3], area_ext[4])

pop_list <- list()
# Simulation Neyman-Scott
for(p in 1:n_parents){
  n_ind <- rpois(1, lambda_offspring)
  # Les enfants sont distribués autour des parents (Processus de Thomas)
  x <- rnorm(n_ind, parents_x[p], sigma_cluster)
  y <- rnorm(n_ind, parents_y[p], sigma_cluster)
  pop_list[[p]] <- data.frame(id = paste0("p", p, "_i", 1:n_ind), 
                              x = x, y = y, year = 1, alive = 1)
}
pop <- do.call(rbind, pop_list)

# --- Étape 2 : Évolution temporelle (Années 2 à T) ---
full_pop <- pop
for(t in 2:n_years){
  prev_year <- full_pop[full_pop$year == (t-1) & full_pop$alive == 1, ]
  
  # Survie
  prev_year$alive <- rbinom(nrow(prev_year), 1, phi)
  
  # Mouvement des survivants
  survivors <- prev_year[prev_year$alive == 1, ]
  survivors$x <- survivors$x + rnorm(nrow(survivors), 0, move_sd)
  survivors$y <- survivors$y + rnorm(nrow(survivors), 0, move_sd)
  survivors$year <- t
  
  full_pop <- rbind(full_pop, survivors)
}

# Nettoyage (rester dans les limites)
full_pop$x <- pmax(area_ext[1], pmin(area_ext[2], full_pop$x))
full_pop$y <- pmax(area_ext[3], pmin(area_ext[4], full_pop$y))



################################################################################


# Supposons une grille de pièges (detectors)
traps <- expand.grid(x = seq(10, 90, by=10), y = seq(10, 90, by=10))
sigma_det <- 3  # Rayon de détection de l'animal
p0 <- 0.2        # Probabilité de détection de base

# Boucle simplifiée de détection
detections <- list()
for(i in 1:nrow(full_pop)) {
  if(full_pop$alive[i] == 1) {
    # Distance aux pièges
    dists <- sqrt((full_pop$x[i] - traps$x)^2 + (full_pop$y[i] - traps$y)^2)
    prob <- p0 * exp(-dists^2 / (2 * sigma_det^2))
    det_count <- rbinom(nrow(traps), 1, prob)
    if(sum(det_count) > 0) {
      detections[[length(detections)+1]] <- data.frame(id=full_pop$id[i], year=full_pop$year[i], trap=which(det_count>0))
    }
  }
}
obs_data <- do.call(rbind, detections)