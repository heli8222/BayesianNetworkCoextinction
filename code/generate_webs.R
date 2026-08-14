# Random generator seed can be passed as a command line argument:
clargs <- commandArgs(trailingOnly = TRUE)
if (length(clargs) > 0 && !is.na(as.integer(clargs[1]))) set.seed(as.integer(clargs[1]))


suppressPackageStartupMessages(library(tidyverse))



# Functions ---------------------------------------------------------------------------

# Stratified Model
create_strat_model <- function(trophic_levels, species_per_level, mixing_parameter) {
  # Mixing matrix, with prob. of interaction between species in level i and level j:
  mixing_matrix <- matrix(0, nrow = trophic_levels, ncol = trophic_levels)
  mixing_matrix[row(mixing_matrix) == col(mixing_matrix) + 1] <- 1
  # Kronecker product with matrix of 1s, to expand to full size with 1 at every potential
  # interaction and 0 elsewhere:
  adj_matrix <- mixing_matrix %x% matrix(1, species_per_level, species_per_level)
  # Multiply every entry by 0 (with probability 1 - mixing_parameter) or 1 (with
  # probability mixing_parameter) randomly:
  s <- trophic_levels * species_per_level
  adj_matrix * matrix(sample(0:1, size = s*s, replace = TRUE,
                             prob = c(1 - mixing_parameter, mixing_parameter)), s, s)
}


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


# make_edg(): function to construct an edge list from an adjacency matrix
# input: S x S matrix for positive integer S
# output: a x 2 matrix for positive integer a
make_edg <- function(adjacency_matrix) {
  n <- sum(adjacency_matrix) # store the number of edges in an adjacency matrix
  edge_list <- matrix(NA, nrow = n, ncol = 2) # initialize matrix of network edges
  colnames(edge_list) <- c("consumer", "resource") # name matrix columns
  row_index <- 1 # initialize counter to track rows of the edge list
  for (i in 1:nrow(adjacency_matrix)) { # loop through all rows of adjacency matrix
    for (j in 1:ncol(adjacency_matrix)) { # loop through all columns of adjacency matrix
      if (adjacency_matrix[i,j] == 1) { # add i,j to edge list if species i eats j
        edge_list[row_index, 1] <- i
        edge_list[row_index, 2] <- j
        row_index <- row_index + 1 # increment row_index to the next row of edge list
      }
    }
  }
  as_tibble(edge_list)
}


# decycle(): function to decycle an adjacency matrix
# input: S x S adjacency matrix
# output: S x S matrix with cycles removed
decycle <- function(adjacency_matrix) {
  DFSCOLOR <- numeric(0)
  DFSBACKEDGE <- numeric(0)
  ORDERVERTICES <- numeric(0)
  DFSVisit <- function(adjacency_matrix, i) {
    DFSCOLOR[i] <<- 1
    for (j in 1:nrow(adjacency_matrix)) {
      if(adjacency_matrix[i,j] != 0) {
        if (DFSCOLOR[j] == 0) {
          DFSVisit(adjacency_matrix,j)
        } else {
          if(DFSCOLOR[j] == 1) {
            DFSBACKEDGE[i,j] <<- 1 # It's a back edge: list for removal
          }
        }
      }
    }
    DFSCOLOR[i] <<- 2
    ORDERVERTICES <<- c(i, ORDERVERTICES)
  }
  run_DFS <- function(adjacency_matrix) {
    S <- nrow(adjacency_matrix)
    DFSCOLOR <<- rep(0, S)
    DFSBACKEDGE <<- matrix(0, S, S)
    ORDERVERTICES <<- numeric(0)
    for (i in 1:S) if (DFSCOLOR[i] == 0) DFSVisit(adjacency_matrix, i)
    adjacency_matrix - DFSBACKEDGE
  }
  run_DFS(adjacency_matrix)
}


find_strat_web <- function(trophic_levels, species_per_level, mixing_parameter) {
  repeat {
    candidate_web <- create_strat_model(trophic_levels,
                                        species_per_level,
                                        mixing_parameter)
    candidate_graph <- igraph::graph_from_adjacency_matrix(t(candidate_web),
                                                           mode = "directed")
    if (length(isolates(candidate_web)) == 0 &
        length(identicals(candidate_web)) == 0 &
        igraph::is_connected(candidate_graph, mode = "weak")) break
  }
  candidate_web
}


find_niche_web <- function(S, C, tolerance, C_factor) {
  repeat {
    candidate_web <- decycle(t(ATNr::create_niche_model(S, C_factor * C)))
    candidate_rows <- nrow(candidate_web)
    candidate_graph <- igraph::graph_from_adjacency_matrix(t(candidate_web),
                                                           mode = "directed")
    candidate_web_conn <- sum(candidate_web) / (candidate_rows * (candidate_rows - 1))
    connectance_difference <- abs(candidate_web_conn - C) / C
    if (length(isolates(candidate_web)) == 0 &
        length(identicals(candidate_web)) == 0 &
        igraph::is_connected(candidate_graph, mode = "weak") &
        connectance_difference < tolerance) break
  }
  candidate_web
}


strat_maker <- function(pars) {
  if (pars$web_ID %% 10 == 0) print(pars$web_ID) # Simple progress indicator
  find_strat_web(pars$trophic_levels, pars$species_per_level, pars$mixing_parameter)
}


niche_maker <- function(pars) {
  if (pars$web_ID %% 10 == 0) print(pars$web_ID) # Simple progress indicator
  find_niche_web(pars$S, pars$C, pars$tolerance, pars$C_factor)
}


generate_webs <- function(param_table, web_maker) {
  param_table %>%
    # Package each set of parameters into a list:
    mutate(pars = pmap(., list)) %>%
    # Create adjacency matrices, with the parameters in each row of the table:
    mutate(adj_mat = map(pars, web_maker)) %>%
    # Extract edge & node lists from adjacency matrices:
    mutate(edge_list = map(adj_mat, make_edg)) %>%
    # Remove the `pars` column that is no longer needed:
    select(-pars)
}


process_empirical_web <- function(node_list_path, edge_list_path) {
  node_list <- read_csv(node_list_path)
  edge_list <- read_csv(edge_list_path)
  adj_mat <- matrix(0, nrow(node_list), nrow(node_list))
  rownames(adj_mat) <- node_list$SpeciesID
  colnames(adj_mat) <- node_list$SpeciesID
  for (i in 1:nrow(edge_list)) {
    adj_mat[as.character(edge_list$ConsumerSpeciesID[i]),
            as.character(edge_list$ResourceSpeciesID[i])] <- 1
  }
  ES_vector <- rownames(adj_mat)[as.numeric(rownames(adj_mat)) > 599]
  species_vector <- rownames(adj_mat)[!(rownames(adj_mat) %in% ES_vector)]
  adj_mat <- adj_mat[species_vector, species_vector]
  isolated_species_to_remove <- isolates(adj_mat) # get vector of isolated species
  while (length(isolated_species_to_remove) > 0) {
    # remove isolated species from adjacency matrix:
    adj_mat <- adj_mat[-isolated_species_to_remove[1], -isolated_species_to_remove[1]]
    # update vector of isolated species:
    isolated_species_to_remove <- which(colSums(adj_mat) + rowSums(adj_mat) == 0)
  }
  # vector of trophically identical species:
  identical_species_to_remove <- identicals(adj_mat)
  # continue removing while there are identical species:
  while(length(identical_species_to_remove) > 0) {
    # remove first trophically identical species from adjacency matrix:
    adj_mat <- adj_mat[-identical_species_to_remove[1], -identical_species_to_remove[1]]
    # update vector of trophically identical species:
    identical_species_to_remove <- identicals(adj_mat)
  }
  adj_mat <- decycle(adj_mat)
  adj_mat
}


add_ES <- function(adj_mat, edge_list_path) {
  edge_list <- read_csv(edge_list_path) %>%
    filter(ConsumerSpeciesID > 599)
  ES_nodes <- unique(edge_list$ConsumerSpeciesID)
  new_adj_mat <- cbind(adj_mat, matrix(0, nrow = nrow(adj_mat),
                                       ncol = length(ES_nodes),
                                       dimnames = list(rownames(adj_mat), ES_nodes)))
  new_adj_mat <- rbind(new_adj_mat, matrix(0, nrow = length(ES_nodes), ncol = ncol(new_adj_mat),
                                           dimnames = list(ES_nodes, colnames(new_adj_mat))))
  for (i in 1:nrow(edge_list)) {
    if (as.character(edge_list$ResourceSpeciesID[i]) %in% colnames(new_adj_mat)) {
      new_adj_mat[as.character(edge_list$ConsumerSpeciesID[i]),
                  as.character(edge_list$ResourceSpeciesID[i])] <- 1
    }
  }
  new_adj_mat
}



# Food web generation -----------------------------------------------------------------

# Generate stratified food webs:
expand_grid(replicate = 1:50,
            trophic_levels = 2:10,
            species_per_level = 10,
            mixing_parameter = 0.5) %>%
  rowid_to_column("web_ID") %>%
  generate_webs(strat_maker) %>%
  write_rds("data/stratified_webs.rds", compress = "xz")


# Generate niche food webs:
expand_grid(replicate = 1:50,
            S = c(25, 50, 100, 150, 200),
            C = c(0.05, 0.1, 0.15, 0.2),
            tolerance = 0.05) %>%
  # Set connectance overshoot factor to find webs with the right connectance faster:
  mutate(C_factor = case_match(S, 25~1.02, 50~1.02, 100~1.03, 150~1.04, 200~1.05)) %>%
  rowid_to_column("web_ID") %>%
  generate_webs(niche_maker) %>%
  write_rds("data/niche_webs.rds", compress = "xz")


# Generate empirical food webs:
empirical_webs <- tibble(
  node_list_path = Sys.glob("data/empirical_webs/*_nodes.csv"),
  edge_list_path = Sys.glob("data/empirical_webs/*_edges.csv")
) %>%
  rowid_to_column("web_ID") %>%
  mutate(replicate = 0) %>%
  mutate(tolerance = 0.05) %>%
  mutate(C_factor = 1) %>%
  mutate(adj_mat = NA) %>%
  mutate(edge_list = NA)
for (i in 1:nrow(empirical_webs)) {
  empirical_webs$adj_mat[i] <- list(process_empirical_web(empirical_webs$node_list_path[i],
                                                          empirical_webs$edge_list_path[i]))
  empirical_webs$edge_list[i] <- list(make_edg(empirical_webs$adj_mat[[i]]))
}

empirical_webs <- empirical_webs %>%
  mutate(S = sapply(adj_mat, nrow)) %>%
  mutate(C = sapply(adj_mat, function(x) sum(x)/(nrow(x)*(nrow(x) - 1)))) %>%
  select(web_ID, replicate, S, C, tolerance, C_factor, adj_mat, edge_list)
empirical_niche_webs <- expand_grid(replicate = 1:50,
                                    SC = paste(empirical_webs$S, empirical_webs$C, sep = "_"),
                                    tolerance = 0.05,
                                    C_factor = 1.02) %>%
  separate(SC, into = c("S", "C"), sep = "_") %>%
  mutate(S = as.numeric(S)) %>%
  mutate(C = as.numeric(C)) %>%
  rowid_to_column("web_ID") %>%
  generate_webs(niche_maker)
write_rds(rbind(empirical_webs, empirical_niche_webs) %>%
            select(-web_ID) %>%
            rowid_to_column("web_ID"),
          "data/empirical_webs.rds", compress = "xz")


# Generate empirical food webs with ES nodes included:
empirical_webs_with_ES <- empirical_webs %>%
  mutate(edge_list_path = Sys.glob("data/empirical_webs/*_edges.csv"),
         temp_adj = adj_mat,
         temp_edge = edge_list)
for (i in 1:nrow(empirical_webs_with_ES)) {
  empirical_webs_with_ES$adj_mat[i] <- list(add_ES(empirical_webs_with_ES$temp_adj[[i]],
                                                   empirical_webs_with_ES$edge_list_path[i]))
  empirical_webs_with_ES$edge_list[i] <- list(make_edg(empirical_webs_with_ES$adj_mat[[i]]))
}
empirical_webs_with_ES <- empirical_webs_with_ES %>%
  select(web_ID, replicate, S, C, tolerance, C_factor, adj_mat, edge_list)
write_rds(empirical_webs_with_ES %>%
            select(-web_ID) %>%
            rowid_to_column("web_ID"),
          "data/empirical_webs_ES.rds", compress = "xz")