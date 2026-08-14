suppressPackageStartupMessages(library(tidyverse))



# Functions ---------------------------------------------------------------------------

# Calculate prey-averaged trophic level for each species from an adjacency matrix
# input: S x S adjacency matrix
# output: length S vector of trophic levels
trophic_levels <- function(adjacency_matrix) {
  A <- adjacency_matrix
  A <- A / rowSums(A)
  A[is.nan(A)] <- 0
  S <- nrow(A)
  solve(diag(S) - A, rep(1, S))
}


infer_trophic_level <- function(outcome_data) {
  rename(outcome_data, x = threat_magnitude, y = excess_extinctions) %>%
    nls(y ~ (m*x - (1-x)*(1-(1-x)^m))/(m*x^2), data = ., start = list(m = 0.01)) %>%
    summary() %>%
    `$`(coefficients) %>%
    `[`(1)
}


aggregate_trophic_levels <- function(data, size, conn, slope, int) {
  expand.grid(m = int + slope*data$actual_m[data$S == size & data$C == conn],
              threat_magnitude = seq(0.1, 0.9, 0.1)) %>%
    mutate(y = (m*threat_magnitude - (1-threat_magnitude)*(1-(1-threat_magnitude)^m))/(m*threat_magnitude^2)) %>%
    summarise(excess_extinctions = mean(y), .by = c(threat_magnitude)) %>%
    infer_trophic_level()
}


# Calculate omnivory index for a community from an adjancency matrix
# input: S x S adjacency matrix
# output: non-negative real number
get_omnivory <- function(adjacency_matrix) {
  TL_vector <- trophic_levels(adjacency_matrix)
  matrix <- t(t(adjacency_matrix)*TL_vector)
  matrix[matrix == 0] <- NA
  sd_vector <- matrixStats::rowSds(matrix, na.rm = T)
  sd_vector[is.na(sd_vector)] <- 0
  return(mean(sd_vector))
}



# Obtain parameters -------------------------------------------------------------------

# Calculate actual and best-fit parameters for stratified webs:
read_rds("data/stratified_outcome.rds") %>%
  filter(alpha == 1 & beta == 1 & threat_type == "deterministic") %>%
  nest(data = !web_ID & !trophic_levels & !alpha & !beta & !threat_type) %>%
  mutate(actual_m = map_dbl(data, function(x) max(trophic_levels(x$adj_mat[[1]])))) %>%
  mutate(inferred_m = map_dbl(data, infer_trophic_level)) %>%
  write_rds("data/stratified_param_fit.rds", compress = "xz")


# Calculate actual and best-fit parameters for niche webs:
read_rds("data/niche_outcome.rds") %>%
  filter(alpha == 1 & beta == 1 & threat_type == "deterministic") %>%
  nest(data = !web_ID & !S & !C & !alpha & !beta & !threat_type) %>%
  mutate(actual_m = map_dbl(data, function(x) max(trophic_levels(x$adj_mat[[1]])))) %>%
  mutate(inferred_m = map_dbl(data, infer_trophic_level)) %>%
  write_rds("data/niche_param_fit.rds", compress = "xz")


# Calculate adjusted parameters for niche webs:
model <- read_rds("data/niche_param_fit.rds") %>%
  lm(inferred_m ~ actual_m, data = .)
read_rds("data/niche_param_fit.rds") %>%
  group_by(S, C) %>%
  summarise(mean_actual_m = aggregate_trophic_levels(data = ., size = S, conn = C,
                                                     slope = 1, int = 0),
            mean_corrected_m = aggregate_trophic_levels(data = ., size = S, conn = C,
                                                        slope = model$coefficients[[2]],
                                                        int = model$coefficients[[1]])) %>%
  ungroup() %>%
  write_rds("data/niche_param_adjusted.rds", compress = "xz")


# Calculate actual and best-fit parameters for empirical webs:
read_rds("data/empirical_outcome.rds") %>%
  filter(alpha == 1 & beta == 1 & threat_type == "deterministic") %>%
  nest(data = !web_ID & !S & !C & !alpha & !beta & !threat_type) %>%
  mutate(actual_m = map_dbl(data, function(x) max(trophic_levels(x$adj_mat[[1]])))) %>%
  mutate(inferred_m = map_dbl(data, infer_trophic_level)) %>%
  write_rds("data/empirical_param_fit.rds", compress = "xz")


# Calculate omnivory indices for empirical webs:
read_rds("data/empirical_webs.rds") %>%
  filter(replicate == 0) %>%
  mutate(omnivory = sapply(adj_mat, get_omnivory)) %>%
  select(c(web_ID, omnivory)) %>%
  cbind({read_rds("data/empirical_param_fit.rds") %>%
      select(-c(data)) %>%
      filter(web_ID < 13) %>%
      mutate(deviation = inferred_m / actual_m) %>%
      pull(deviation)}) %>%
  rename(deviation = `{`) %>%
  write_rds("data/empirical_omnivory_deviation.rds", compress = "xz")