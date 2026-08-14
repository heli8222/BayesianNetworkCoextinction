suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(igraph))



# Functions ---------------------------------------------------------------------------

# isolates(): function to detect identity of single isolated species in adjacency matrix
# input: S x S matrix for positive integer S
# output: vector of length a for non-negative integer a
isolates <- function(adjacency_matrix) {
  # initialize a vector to store indices of isolated species:
  isolated_species <- numeric(length = 0)
  # loop through all species to identify isolated species:
  for (i in 1:nrow(adjacency_matrix)) {
    # add species to vector if it eats no one and no one eats it:
    if (sum(adjacency_matrix[i,]) == 0 & sum(adjacency_matrix[,i]) == 0) {
      isolated_species <- append(isolated_species, i)
    }
  }
  isolated_species
}


# identicals(): detect identity of single trophically identical species in adj. matrix
# input: S x S matrix for positive integer S
# output: vector of length a for non-negative integer a
identicals <- function(adjacency_matrix) {
  # initialize a vector to store indices of identical species:
  identical_species <- numeric(length = 0)
  # loop through all species-pairs to identify identical species:
  for (i in 1:nrow(adjacency_matrix)) {
    # initialize a vector to store all species that species i eats:
    i_eats <- which(adjacency_matrix[i,] == 1)
    # initialize a vector to store all species that eat species i:
    i_eaten <- which(adjacency_matrix[,i] == 1)
    for (j in 1:nrow(adjacency_matrix)) {
      if (j != i){
        # initialize a vector to store all species that species j eats:
        j_eats <- which(adjacency_matrix[j,] == 1)
        # initialize a vector to store all species that eats species j:
        j_eaten <- which(adjacency_matrix[,j] == 1)
        # add species j to vector if it has the same eats and eaten vectors as species i:
        if (identical(i_eats, j_eats) & identical(i_eaten, j_eaten)) {
          identical_species <- append(identical_species, j)
        }
      }
    }
  }
  identical_species
}



# Generate web lists from adjacency matrices ------------------------------------------

stratified_web_df <- read_rds("data/stratified_webs.rds")
niche_web_df <- read_rds("data/niche_webs.rds")
empirical_web_df <- read_rds("data/empirical_webs.rds")


stratified_webs <- lapply(stratified_web_df$adj_mat, t) %>%
  lapply(graph_from_adjacency_matrix, mode = "directed")
niche_webs <- lapply(niche_web_df$adj_mat, t) %>%
  lapply(graph_from_adjacency_matrix, mode = "directed")
empirical_webs <- lapply(empirical_web_df$adj_mat, t) %>%
  lapply(graph_from_adjacency_matrix, mode = "directed")



# Check that webs are connected -------------------------------------------------------

stratified_connected <- sapply(stratified_webs, is_connected, mode = "weak")
niche_connected <- sapply(niche_webs, is_connected, mode = "weak")
empirical_connected <- sapply(empirical_webs, is_connected, mode = "weak")


all(stratified_connected)
all(niche_connected)
all(empirical_connected)



# Check that webs do not have isolated species ----------------------------------------

stratified_isolates <- lapply(stratified_web_df$adj_mat, isolates) %>%
  sapply(length)
niche_isolates <- lapply(niche_web_df$adj_mat, isolates) %>%
  sapply(length)
empirical_isolates <- lapply(empirical_web_df$adj_mat, isolates) %>%
  sapply(length)


sum(stratified_isolates)
sum(niche_isolates)
sum(empirical_isolates)



# Check that webs do not have trophically identical species ---------------------------

stratified_identicals <- lapply(stratified_web_df$adj_mat, identicals) %>%
  sapply(length)
niche_identicals <- lapply(niche_web_df$adj_mat, identicals) %>%
  sapply(length)
empirical_identicals <- lapply(empirical_web_df$adj_mat, identicals) %>%
  sapply(length)


sum(stratified_identicals)
sum(niche_identicals)
sum(empirical_identicals)



# Check that niche webs fall within connectance tolerance -----------------------------

niche_connectance <- sapply(niche_web_df$adj_mat,
                            function(x) { sum(x) / (nrow(x) * (nrow(x) - 1)) })
niche_connectance_difference <- abs(niche_connectance - niche_web_df$C) / niche_web_df$C
niche_within_tolerance <- as.integer(ifelse(niche_connectance_difference < 0.05, 1, 0))


empirical_connectance <- sapply(empirical_web_df$adj_mat,
                                function(x) { sum(x) / (nrow(x) * (nrow(x) - 1)) })
empirical_connectance_difference <- abs(empirical_connectance - empirical_web_df$C) / empirical_web_df$C
empirical_within_tolerance <- as.integer(ifelse(empirical_connectance_difference < 0.05, 1, 0))


all(niche_within_tolerance)
all(empirical_within_tolerance)



# Check that webs are not identical to each other -------------------------------------

stratified_identical_webs <- matrix(0,
                                    nrow = length(stratified_webs),
                                    ncol = length(stratified_webs))
for ( i in 1:(nrow(stratified_identical_webs) - 1) ) {
  for ( j in (i + 1):nrow(stratified_identical_webs) ) {
    stratified_identical_webs[i, j] <- ifelse(identical(stratified_web_df$adj_mat[[i]],
                                                        stratified_web_df$adj_mat[[j]]),
                                              1, 0)
  }
}
sum(stratified_identical_webs)


niche_identical_webs <- matrix(0,
                               nrow =length(niche_webs),
                               ncol = length(niche_webs))
for ( i in 1:(nrow(niche_identical_webs) - 1) ) {
  for ( j in (i + 1):nrow(niche_identical_webs) ) {
    niche_identical_webs[i, j] <- ifelse(identical(niche_web_df$adj_mat[[i]],
                                                   niche_web_df$adj_mat[[j]]),
                                         1, 0)
  }
}
sum(niche_identical_webs)


empirical_identical_webs <- matrix(0,
                                   nrow = length(empirical_webs),
                                   ncol = length(empirical_webs))
for ( i in 1:(nrow(empirical_identical_webs) - 1) ) {
  for ( j in (i + 1):nrow(empirical_identical_webs) ) {
    empirical_identical_webs[i, j] <- ifelse(identical(empirical_web_df$adj_mat[[i]],
                                                       empirical_web_df$adj_mat[[j]]),
                                             1, 0)
  }
}
sum(empirical_identical_webs)