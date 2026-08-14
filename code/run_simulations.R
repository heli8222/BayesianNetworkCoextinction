# Number of Monte Carlo iterations can be passed as a command line argument:
clargs <- commandArgs(trailingOnly = TRUE)
if (length(clargs) > 0 && !is.na(as.integer(clargs[1]))) {
  nreps <- as.integer(clargs[1])
} else {
  nreps <- 1e3L
}


suppressPackageStartupMessages(library(tidyverse))



# Functions ---------------------------------------------------------------------------

# Source C file with a single function, persistence_mc(). Performs Monte Carlo
# simulation of Bayesian network, returning species' marginal persistence probabilities.
# Input:
# - A: adjacency matrix of the food web
# - Pb: vector of baseline extinction probabilities
# - alpha: first parameter of the Beta distribution
# - beta: second parameter of the Beta distribution
# - nreps: number of iterations for Monte Carlo sim of Bayesian network
# Output:
# - vector of marginal persistence probabilities
Rcpp::sourceCpp("code/bayesian_nw_monte_carlo.cpp")


# Create adjacency matrix from a data frame of edges
# Input:
# - edgelist: Data frame with two columns of taxon names, where consumers
#             are in column 1 and resources in column 2. Assumes that the
#             corresponding web is acyclic.
# Output:
# - The adjacency matrix, with A[i,j] = 1 if species i eats j and 0 otherwise.
#   Species are topologically sorted to make the matrix lower triangular.
#   The rows and columns of the matrix are labeled by taxon names.
create_adj_matrix <- function(edgelist) {
  # Create adjacency matrix A; links point from resource to consumer, so
  # the two columns of the edge list data frame are flipped:
  A <- igraph::graph_from_data_frame(d = edgelist[,2:1]) %>% # Create igraph graph
    igraph::as_adjacency_matrix() %>% # Convert to igraph adjacency matrix
    as.matrix() %>% # Convert to a regular matrix
    t() # Transpose result, so A[i,j] is 1 if i eats j and 0 otherwise
  A <- A[order(rowSums(A)),order(rowSums(A))] # Make sure producers are the first rows
  # Find a sorting of A's rows and columns to make A lower triangular:
  s <- igraph::graph_from_adjacency_matrix(A) %>% igraph::topo_sort(mode = "in")
  A[s,s] # Return sorted, lower triangular adjacency matrix A
}


# Obtain table of persistence probabilities
# Input:
# - edges: data frame with two columns of taxon names, where consumers
#          are in column 1 and resources in column 2
# - threats: data frame with its 1st column corresponding to taxon names,
#            and the rest to specific threats (and the baseline extinction
#            probabilities for each taxon listed under them)
# - alpha: 1st parameter of beta distribution, governing consumer response
#          to prey loss
# - beta: 2nd parameter of beta distribution, governing consumer response
#         to prey loss
# - nreps: number of iterations for Monte Carlo sim of the Bayesian network
# Output:
# - Data frame with taxon names, threats, and persistence probabilities
persistence_table <- function(edges, threats, alpha, beta, nreps) {
  A <- create_adj_matrix(unique(edges)) # sorted, lower triangular adj. matrix A
  threats <- threats[match(rownames(A), threats$taxon),] # match row order of A
  threatlist <- threats %>% # get vector of threat names
    select(-taxon) %>% # all columns except `taxon`
    colnames() # extract names of remaining columns
  persistence_vector <- tibble() # persistence table
  for (threat in threatlist) { # for each threat in list of threats:
    Pb <- threats %>% # obtain baseline extinction probs. from table of threats
      pull(threat) %>% # column `threat` as a vector
      replace_na(0) # replace NAs (for taxa not in the threat table) with 0s
    persistence_df <- persistence_mc(A, Pb, alpha, beta, nreps)
    persistence_vector <- tibble( # create data frame with:
      taxon = colnames(A), # taxon names
      threat = threat, # threat names
      persist = colMeans(persistence_df) # persistence
    ) %>%
      bind_rows(persistence_vector) # attach to main persistence table
  }
  list(persistence_df, persistence_vector)
}


# Obtain persistence probabilities given a web and threat file,
# plus model parameters alpha, beta, and nreps
# Input:
# - web: file with web to analyze (with path & extension), in csv format;
#        the column names must be `consumer` and `resource`
# - threatfile: file containing threat information (with path & extension),
#               in csv format; the first column must be "taxon", the rest are
#               threat names
# - alpha: 1st parameter of beta distribution, governing consumer response
#          to prey loss
# - beta: 2nd parameter of beta distribution, governing consumer response
#         to prey loss
# - nreps: number of iterations for Monte Carlo sim of the Bayesian network
# Output:
# - Data frame with taxon names, threats, persistence probabilities,
#   and parameters alpha, beta, and nreps
threat_analysis <- function(web, threatfile, alpha, beta, nreps) {
  edges <- web
  threats <- threatfile %>% # only spp in subweb:
    filter(taxon %in% unique(c(edges$consumer, edges$resource)))
  persistence_table(edges, threats, alpha, beta, nreps)
}


param_combos <- function(web_table) {
  web_table %>%
    drop_na() %>%
    select(web_ID) %>%
    distinct() %>%
    expand_grid(
      threat_type = c(
        "deterministic", # All species have equal baseline probability of extinction (BP)
        "stochastic", # Each species' BP is uniform in: target probability +/- 0.025
        "basal-only" # All basal species have equal BP; all non-basal species have 0 BP
      ),
      threat_magnitude = 1:9 / 10,
      function_ID = 1:3 # Controls combination of alpha and beta parameters
    ) %>%
    # Join with a table containing the three legal alpha-beta combinations:
    left_join(tibble(function_ID = 1:3, alpha = c(1, 1, 2), beta = c(1, 2, 1)),
              by = join_by(function_ID)) %>%
    select(-function_ID) # Remove this now-obsolete column
}


# Helper function to generate threat profile data frames for different threat types:
generate_threat_profile <- function(adj_mat, threat_type, threat_magnitude) {
  case_when(
    threat_type == "deterministic" ~ tibble(
      taxon = 1:nrow(adj_mat),
      threat = threat_magnitude
    ),
    threat_type == "stochastic" ~ tibble(
      taxon = 1:nrow(adj_mat),
      threat = runif(n = nrow(adj_mat),
                     min = threat_magnitude - 0.025,
                     max = threat_magnitude + 0.025)
    ),
    threat_type == "basal-only" ~ tibble(
      taxon = 1:nrow(adj_mat),
      threat = (rowSums(adj_mat) == 0) * threat_magnitude
    )
  )
}


# Helper function to simulate extinctions on a web:
simulate_extinction <- function(sim_ID, threat_type, threat_magnitude, alpha, beta,
                                adj_mat, edge_list, nreps) {
  if (sim_ID %% 100 == 0) print(sim_ID) # Simple progress indicator
  threat_profile <- generate_threat_profile(adj_mat, threat_type, threat_magnitude)
  threat_analysis(edge_list, threat_profile, alpha, beta, nreps)[[2]] %>%
    select(-threat) %>%
    left_join(mutate(threat_profile, taxon = as.character(taxon)), by=join_by(taxon)) %>%
    summarise(direct_extinctions=sum(threat), indirect_extinctions=sum(1 - persist)) %>%
    mutate(excess_extinctions = indirect_extinctions / direct_extinctions)
}


# Helper function to simulate ES extinctions on a web:
simulate_extinction_ES <- function(sim_ID, web_ID, threat_type, threat_magnitude, alpha, beta,
                                   adj_mat, edge_list, nreps) {
  if (sim_ID %% 100 == 0) print(sim_ID) # Simple progress indicator
  threat_profile <- generate_threat_profile(adj_mat, threat_type, threat_magnitude)
  threat_analysis(edge_list, threat_profile, alpha, beta, nreps)[[2]] %>%
    select(-threat) %>%
    left_join(mutate(threat_profile, taxon = as.character(taxon)), by=join_by(taxon)) %>%
    filter(taxon %in% ES_node_table$node_ID[ES_node_table$web_ID == web_ID]) %>%
    mutate(excess_extinctions = (1 - persist) / threat)
}


# Helper function to perform simulations on all webs in parameter table:
process_webs <- function(web_table) {
  left_join(param_combos(web_table), web_table, by = join_by(web_ID)) %>%
    rowid_to_column("sim_ID") %>% # For the progress indicator in `simulate_extinction`
    mutate(extinction = pmap(list(sim_ID, threat_type, threat_magnitude, alpha, beta,
                                  adj_mat, edge_list), simulate_extinction,
                             nreps = nreps)) %>%
    unnest(extinction) %>%
    select(-sim_ID) # This was only needed for the progress indicator; remove
}


# Helper function to perform ES simulations on all webs in parameter table:
process_webs_ES <- function(web_table) {
  left_join(param_combos(web_table), web_table, by = join_by(web_ID)) %>%
    rowid_to_column("sim_ID") %>% # For the progress indicator in `simulate_extinction`
    mutate(extinction = pmap(list(sim_ID, web_ID, threat_type, threat_magnitude, alpha, beta,
                                  adj_mat, edge_list), simulate_extinction_ES,
                             nreps = nreps)) %>%
    unnest(extinction) %>%
    left_join(ES_node_table, by = join_by(web_ID, taxon == node_ID)) %>%
    select(-sim_ID) # This was only needed for the progress indicator; remove
}



# Simulate extinctions ----------------------------------------------------------------

# Simulate stratified webs and save results:
read_rds("data/stratified_webs.rds") %>%
  process_webs() %>%
  write_rds("data/stratified_outcome.rds", compress = "xz") %>%
  select(-c(adj_mat, edge_list)) %>%
  write_rds("data/stratified_outcome_small.rds", compress = "xz")


# Simulate niche webs and save results:
read_rds("data/niche_webs.rds") %>%
  process_webs() %>%
  write_rds("data/niche_outcome.rds", compress = "xz") %>%
  select(-c(adj_mat, edge_list)) %>%
  write_rds("data/niche_outcome_small.rds", compress = "xz")


# Simulate empirical webs and save results:
read_rds("data/empirical_webs.rds") %>%
  process_webs() %>%
  write_rds("data/empirical_outcome.rds", compress = "xz") %>%
  select(-c(adj_mat, edge_list)) %>%
  write_rds("data/empirical_outcome_small.rds", compress = "xz")


# Simulate ES on empirical webs and save results:
ES_node_table <- read_rds("data/ES_node_ID_table.rds")
read_rds("data/empirical_webs_ES.rds") %>%
  process_webs_ES() %>%
  write_rds("data/empirical_ES_outcome.rds", compress = "xz") %>%
  select(-c(adj_mat, edge_list)) %>%
  write_rds("data/empirical_ES_outcome_small.rds", compress = "xz")