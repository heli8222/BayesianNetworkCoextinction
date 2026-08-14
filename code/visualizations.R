suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(patchwork))
suppressPackageStartupMessages(library(figpatch))



# Functions ---------------------------------------------------------------------------

overlay_lines_EE_e <- function(m, color, type, lwd, a, b) {
  if (a == 1 & b == 1) {
    return(geom_function(fun = function(x) ((x*(m+1) + (1-x)^(m+1) - 1)/(m*x^2)),
                         col = color, alpha = 0.5, linetype = type, linewidth = lwd))
  }
  
  if (a == 2 & b == 1) {
    return(geom_function(fun = function(x) ((m*x + log(((1-x)^2 - x^2*exp((1-2*x)/2))/((1-x)^2 - x^2*exp((m-0.5)*(2*x-1)))))/(m*x*(1-x))),
                         col = color, alpha = 0.5, linetype = type, linewidth = lwd))
  }
  
  if (a == 1 & b == 2) {
    return(geom_function(fun = function(x) ((log((x*(2-x)*exp(m-0.5) + (1-x)^2)/(x*(2-x)*exp(-0.5) + (1-x)^2)) - m*x)/(m*x*(1-x))),
                         col = color, alpha = 0.5, linetype = type, linewidth = lwd))
  }
}


overlay_lines_F_e <- function(m, color, a, b) {
  if (a == 1 & b == 1) {
    return(geom_function(fun = function(x) ((x/(1-x)) * (((x*(m+1) + (1-x)^(m+1) - 1)/(m*x^2)) - 1)),
                         col = color, alpha = 0.5, linewidth  = 0.5))
  }
  
  if (a == 2 & b == 1) {
    return(geom_function(fun = function(x) ((x/(1-x)) * (((m*x + log(((1-x)^2 - x^2*exp((1-2*x)/2))/((1-x)^2 - x^2*exp((m-0.5)*(2*x-1)))))/(m*x*(1-x))) - 1)),
                         col = color, alpha = 0.5, linewidth = 0.5))
  }
  
  if (a == 1 & b == 2) {
    return(geom_function(fun = function(x) ((x/(1-x)) * ((log((x*(2-x)*exp(m-0.5) + (1-x)^2)/(x*(2-x)*exp(-0.5) + (1-x)^2)) - m*x)/(m*x*(1-x)) - 1)),
                         col = color, alpha = 0.5, linewidth = 0.5))
  }
}


overlay_lines_EE_TL <- function(e, color, a, b) {
  if (a == 1 & b == 1) {
    return(geom_function(fun = function(x) ((e*(x+1) + (1-e)^(x+1) - 1)/(x*e^2)),
                         col = color, alpha = 0.5, linewidth = 0.5))
  }
  
  if (a == 2 & b == 1) {
    return(geom_function(fun = function(x) ((x*e + log(((1-e)^2 - e^2*exp((1-2*e)/2))/((1-e)^2 - e^2*exp((x-0.5)*(2*e-1)))))/(x*e*(1-e))),
                         col = color, alpha = 0.5, linewidth = 0.5))
  }
  
  if (a == 1 & b == 2) {
    return(geom_function(fun = function(x) ((log((e*(2-e)*exp(x-0.5) + (1-e)^2)/(e*(2-e)*exp(-0.5) + (1-e)^2)) - x*e)/(x*e*(1-e))),
                         col = color, alpha = 0.5, linewidth = 0.5))
  }
}


overlay_lines_F_TL <- function(e, color, a, b) {
  if (a == 1 & b == 1) {
    return(geom_function(fun = function(x) ((e/(1-e)) * (((e*(x+1) + (1-e)^(x+1) - 1)/(x*e^2)) - 1)),
                         col = color, alpha = 0.5, linewidth  = 0.5))
  }
  
  if (a == 2 & b == 1) {
    return(geom_function(fun = function(x) ((e/(1-e)) * (((x*e + log(((1-e)^2 - e^2*exp((1-2*e)/2))/((1-e)^2 - e^2*exp((x-0.5)*(2*e-1)))))/(x*e*(1-e))) - 1)),
                         col = color, alpha = 0.5, linewidth = 0.5))
  }
  
  if (a == 1 & b == 2) {
    return(geom_function(fun = function(x) ((e/(1-e)) * ((log((e*(2-e)*exp(x-0.5) + (1-e)^2)/(e*(2-e)*exp(-0.5) + (1-e)^2)) - x*e)/(x*e*(1-e)) - 1)),
                         col = color, alpha = 0.5, linewidth = 0.5))
  }
}


EE_strat_e <- function(data, threatType, a, b) {
  data %>%
    filter(threat_type == threatType & alpha == a & beta == b) %>%
    {
      ggplot(., aes(x = threat_magnitude, y = excess_extinctions, col = trophic_levels)) +
        geom_point(size = 1, alpha = 0.5, position = position_jitter(width = 0.01, seed = 63)) +
        {
          if (threatType != "basal-only") {
            unlist(pmap(list(2:10, viridis::plasma(9, end = 0.9), "solid", 0.5, a, b),
                        overlay_lines_EE_e))
          }
        } +
        coord_cartesian(xlim = c(0, 1), ylim = c(1, 1.1*max(.$excess_extinctions))) +
        scale_x_continuous(name = "Direct Extinction Risk (\u03b5)", breaks = seq(0, 1, 0.2)) +
        scale_y_continuous(name = "Excess Extinction (E)") +
        scale_colour_viridis_d(name = "Trophic\nlevels", option = "C", end = 0.9,
                               guide = guide_legend(override.aes = list(alpha = 1),
                                                    reverse = TRUE,
                                                    ncol = 1)) +
        theme_bw()
    }
}


F_strat_e <- function(data, threatType, a, b) {
  data %>%
    filter(threat_type == threatType & alpha == a & beta == b) %>%
    {
      ggplot(., aes(x = threat_magnitude, y = fragility, col = trophic_levels)) +
        geom_point(size = 1, alpha = 0.5, position = position_jitter(width = 0.01, seed = 63)) +
        {
          if (threatType != "basal-only") {
            unlist(pmap(list(2:10, viridis::plasma(9, end = 0.9), a, b),
                        overlay_lines_F_e))
          }
        } +
        coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) +
        scale_x_continuous(name = "Direct Extinction Risk (\u03b5)", breaks = seq(0, 1, 0.2)) +
        scale_y_continuous(name = "Fragility (F)", breaks = seq(0, 1, 0.25)) +
        scale_colour_viridis_d(name = "Trophic\nlevels", option = "C", end = 0.9,
                               guide = guide_legend(override.aes = list(alpha = 1),
                                                    reverse = TRUE,
                                                    ncol = 1)) +
        theme_bw()
    }
}


EE_strat_TL <- function(data, threatType, a, b) {
  data %>%
    mutate(trophic_levels = as.numeric(as.character(trophic_levels))) %>%
    filter(threat_type == threatType & alpha == a & beta == b) %>%
    {
      ggplot(., aes(x = trophic_levels, y = excess_extinctions, col = threat_magnitude)) +
        geom_point(size = 1, alpha = 0.5, position = position_jitter(width = 0.01, seed = 63)) +
        {
          if (threatType != "basal-only") {
            unlist(pmap(list(c(1:4, 4.999, 6:9)/10, viridis::plasma(9, end = 0.9), a, b),
                        overlay_lines_EE_TL))
          }
        } +
        coord_cartesian(xlim = c(2, 10), ylim = c(1, 1.1*max(.$excess_extinctions))) +
        scale_x_continuous(name = "Number of Trophic Levels", breaks = seq(2, 10, 2)) +
        scale_y_continuous(name = "Excess Extinction (E)") +
        scale_colour_viridis_c(name = "Direct Extinction\nRisk (\u03b5)", option = "C", end = 0.9,
                               guide = guide_colorbar(alpha = 1)) +
        theme_bw()
    }
}


F_strat_TL <- function(data, threatType, a, b) {
  data %>%
    mutate(trophic_levels = as.numeric(as.character(trophic_levels))) %>%
    filter(threat_type == threatType & alpha == a & beta == b) %>%
    {
      ggplot(., aes(x = trophic_levels, y = fragility, col = threat_magnitude)) +
        geom_point(size = 1, alpha = 0.5, position = position_jitter(width = 0.01, seed = 63)) +
        {
          if (threatType != "basal-only") {
            unlist(pmap(list(c(1:4, 4.999, 6:9)/10, viridis::plasma(9, end = 0.9), a, b),
                        overlay_lines_F_TL))
          }
        } +
        coord_cartesian(xlim = c(2, 10), ylim = c(0, 1)) +
        scale_x_continuous(name = "Number of Trophic Levels", breaks = seq(2, 10, 2)) +
        scale_y_continuous(name = "Fragility (F)", breaks = seq(0, 1, 0.25)) +
        scale_colour_viridis_c(name = "Direct Extinction\nRisk (\u03b5)", option = "C", end = 0.9,
                               guide = guide_colorbar(alpha = 1)) +
        theme_bw()
    }
}


EE_niche_S <- function(data, threatType, a, b) {
  data %>%
    filter(threat_type == threatType & alpha == a & beta == b) %>%
    ggplot(aes(x = S, y = excess_extinctions, col = threat_magnitude)) +
    geom_point(size = 1, alpha = 0.5, position = position_jitter(width = 0.1, seed = 63)) +
    scale_x_discrete(name = "Food Web Size") +
    scale_y_continuous(name = "Excess Extinction (E)") +
    scale_colour_viridis_c(name = "Direct Extinction\nRisk (\u03b5)", option = "C", end = 0.9,
                           guide = guide_colorbar(alpha = 1)) +
    theme_bw()
}


F_niche_S <- function(data, threatType, a, b) {
  data %>%
    filter(threat_type == threatType & alpha == a & beta == b) %>%
    ggplot(aes(x = S, y = fragility, col = threat_magnitude)) +
    geom_point(size = 1, alpha = 0.5, position = position_jitter(width = 0.1, seed = 63)) +
    scale_x_discrete(name = "Food Web Size") +
    scale_y_continuous(name = "Fragility (F)",
                       limits = c(0, 1), breaks = seq(0, 1, 0.25)) +
    scale_colour_viridis_c(name = "Direct Extinction\nRisk (\u03b5)", option = "C", end = 0.9,
                           guide = guide_colorbar(alpha = 1)) +
    theme_bw()
}


EE_niche_C <- function(data, threatType, a, b) {
  data %>%
    filter(threat_type == threatType & alpha == a & beta == b) %>%
    ggplot(aes(x = C, y = excess_extinctions, col = threat_magnitude)) +
    geom_point(size = 1, alpha = 0.5, position = position_jitter(width = 0.1, seed = 63)) +
    scale_x_discrete(name = "Food Web Connectance") +
    scale_y_continuous(name = "Excess Extinction (E)") +
    scale_colour_viridis_c(name = "Direct Extinction\nRisk (\u03b5)", option = "C", end = 0.9,
                           guide = guide_colorbar(alpha = 1)) +
    theme_bw()
}


F_niche_C <- function(data, threatType, a, b) {
  data %>%
    filter(threat_type == threatType & alpha == a & beta == b) %>%
    ggplot(aes(x = C, y = fragility, col = threat_magnitude)) +
    geom_point(size = 1, alpha = 0.5, position = position_jitter(width = 0.1, seed = 63)) +
    scale_x_discrete(name = "Food Web Connectance") +
    scale_y_continuous(name = "Fragility (F)",
                       limits = c(0, 1), breaks = seq(0, 1, 0.25)) +
    scale_colour_viridis_c(name = "Direct Extinction\nRisk (\u03b5)", option = "C", end = 0.9,
                           guide = guide_colorbar(alpha = 1)) +
    theme_bw()
}


EE_niche_TL <- function(data, params, threatType, a, b) {
  data %>%
    filter(threat_type == threatType & alpha == a & beta == b) %>%
    left_join(select(params, web_ID, actual_m), by = join_by(web_ID)) %>%
    ggplot(aes(x = actual_m, y = excess_extinctions, col = threat_magnitude)) +
    geom_point(size = 1, alpha = 0.5, position = position_jitter(width = 0.01, seed = 63)) +
    scale_x_continuous(name = "Number of Trophic Levels",
                       limits = c(2, 10), breaks = seq(2, 10, 2)) +
    scale_y_continuous(name = "Excess Extinction (E)") +
    scale_colour_viridis_c(name = "Direct Extinction\nRisk (\u03b5)", option = "C", end = 0.9,
                           guide = guide_colorbar(alpha = 1)) +
    theme_bw()
}


F_niche_TL <- function(data, params, threatType, a, b) {
  data %>%
    filter(threat_type == threatType & alpha == a & beta == b) %>%
    left_join(select(params, web_ID, actual_m), by = join_by(web_ID)) %>%
    ggplot(aes(x = actual_m, y = fragility, col = threat_magnitude)) +
    geom_point(size = 1, alpha = 0.5, position = position_jitter(width = 0.01, seed = 63)) +
    scale_x_continuous(name = "Number of Trophic Levels",
                       limits = c(2, 10), breaks = seq(2, 10, 2)) +
    scale_y_continuous(name = "Fragility (F)",
                       limits = c(0, 1), breaks = seq(0, 1, 0.25)) +
    scale_colour_viridis_c(name = "Direct Extinction\nRisk (\u03b5)", option = "C", end = 0.9,
                           guide = guide_colorbar(alpha = 1)) +
    theme_bw()
}


EE_niche_boxplot <- function(data, params, threatType, a, b, conn) {
  data %>%
    filter(threat_type == threatType & alpha == a & beta == b & C == conn) %>%
    ggplot(aes(x = threat_magnitude, y = excess_extinctions, fill = S,
               group = interaction(threat_magnitude, S))) +
    geom_boxplot(linewidth = 0.5, outlier.size = 1) +
    unlist(pmap(list(filter(params, C == conn)$mean_actual_m,
                     viridis::plasma(length(unique(data$S)), end = 0.9),
                     "solid", 0.5, a, b),
                overlay_lines_EE_e)) +
    unlist(pmap(list(filter(params, C == conn)$mean_corrected_m,
                     viridis::plasma(length(unique(data$S)), end = 0.9),
                     "dashed", 0.5, a, b),
                overlay_lines_EE_e)) +
    scale_x_continuous(name = "Direct Extinction Risk (\u03b5)",
                       limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
    scale_y_continuous(name = "Excess Extinction (E)") +
    scale_fill_viridis_d(name = "Food Web\nSize", option = "C", end = 0.9,
                         guide = guide_legend(reverse = TRUE)) +
    labs(title = str_c("C = ", conn)) +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5))
}


EE_empirical_S <- function(data, threatType, a, b) {
  data %>%
    filter(threat_type == threatType & alpha == a & beta == b & replicate == 0) %>%
    ggplot(aes(x = S, y = excess_extinctions, col = threat_magnitude)) +
    geom_point(size = 1, alpha = 0.5, position = position_jitter(width = 2, seed = 63)) +
    scale_x_continuous(name = "Food Web Size") +
    scale_y_continuous(name = "Excess Extinction (E)") +
    scale_colour_viridis_c(name = "Direct Extinction\nRisk (\u03b5)", option = "C", end = 0.9,
                           guide = guide_colorbar(alpha = 1)) +
    theme_bw(base_size = 8)
}


F_empirical_S <- function(data, threatType, a, b) {
  data %>%
    filter(threat_type == threatType & alpha == a & beta == b & replicate == 0) %>%
    ggplot(aes(x = S, y = fragility, col = threat_magnitude)) +
    geom_point(size = 1, alpha = 0.5, position = position_jitter(width = 2, seed = 63)) +
    scale_x_continuous(name = "Food Web Size") +
    scale_y_continuous(name = "Fragility (F)",
                       limits = c(0, 1), breaks = seq(0, 1, 0.25)) +
    scale_colour_viridis_c(name = "Direct Extinction\nRisk (\u03b5)", option = "C", end = 0.9,
                           guide = guide_colorbar(alpha = 1)) +
    theme_bw(base_size = 8)
}


EE_empirical_C <- function(data, threatType, a, b) {
  data %>%
    filter(threat_type == threatType & alpha == a & beta == b & replicate == 0) %>%
    ggplot(aes(x = C, y = excess_extinctions, col = threat_magnitude)) +
    geom_point(size = 1, alpha = 0.5, position = position_jitter(width = 0.0025, seed = 63)) +
    scale_x_continuous(name = "Food Web Connectance") +
    scale_y_continuous(name = "Excess Extinction (E)") +
    scale_colour_viridis_c(name = "Direct Extinction\nRisk (\u03b5)", option = "C", end = 0.9,
                           guide = guide_colorbar(alpha = 1)) +
    theme_bw(base_size = 8)
}


F_empirical_C <- function(data, threatType, a, b) {
  data %>%
    filter(threat_type == threatType & alpha == a & beta == b & replicate == 0) %>%
    ggplot(aes(x = C, y = fragility, col = threat_magnitude)) +
    geom_point(size = 1, alpha = 0.5, position = position_jitter(width = 0.0025, seed = 63)) +
    scale_x_continuous(name = "Food Web Connectance") +
    scale_y_continuous(name = "Fragility (F)",
                       limits = c(0, 1), breaks = seq(0, 1, 0.25)) +
    scale_colour_viridis_c(name = "Direct Extinction\nRisk (\u03b5)", option = "C", end = 0.9,
                           guide = guide_colorbar(alpha = 1)) +
    theme_bw(base_size = 8)
}


EE_empirical_TL <- function(data, params, threatType, a, b) {
  data %>%
    filter(threat_type == threatType & alpha == a & beta == b & replicate == 0) %>%
    left_join(select(params, web_ID, actual_m), by = join_by(web_ID)) %>%
    ggplot(aes(x = actual_m, y = excess_extinctions, col = threat_magnitude)) +
    geom_point(size = 1, alpha = 0.5, position = position_jitter(width = 0.01, seed = 63)) +
    scale_x_continuous(name = "Number of Trophic Levels") +
    scale_y_continuous(name = "Excess Extinction (E)") +
    scale_colour_viridis_c(name = "Direct Extinction\nRisk (\u03b5)", option = "C", end = 0.9,
                           guide = guide_colorbar(alpha = 1)) +
    theme_bw()
}


F_empirical_TL <- function(data, params, threatType, a, b) {
  data %>%
    filter(threat_type == threatType & alpha == a & beta == b & replicate == 0) %>%
    left_join(select(params, web_ID, actual_m), by = join_by(web_ID)) %>%
    ggplot(aes(x = actual_m, y = fragility, col = threat_magnitude)) +
    geom_point(size = 1, alpha = 0.5, position = position_jitter(width = 0.01, seed = 63)) +
    scale_x_continuous(name = "Number of Trophic Levels") +
    scale_y_continuous(name = "Fragility (F)",
                       limits = c(0, 1), breaks = seq(0, 1, 0.25)) +
    scale_colour_viridis_c(name = "Direct Extinction\nRisk (\u03b5)", option = "C", end = 0.9,
                           guide = guide_colorbar(alpha = 1)) +
    theme_bw()
}


correlation_plot <- function(data) {
  lm <- data %>%
    filter(type == "Niche") %>%
    lm(inferred_m ~ actual_m, .)
  data %>%
    ggplot(aes(x = actual_m, y = inferred_m, col = type)) +
    geom_point(data = ~subset(., type == "Niche"), size = 1, alpha = 0.5) +
    geom_point(data = ~subset(., type == "Empirical"), size = 1, alpha = 1) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", linewidth = 0.5, alpha = 0.75) +
    geom_abline(intercept = lm$coefficients[[1]], slope = lm$coefficients[[2]],
                linetype = "solid", linewidth = 0.5, alpha = 0.75) +
    scale_x_continuous(name = "True Max Trophic Level",
                       limits = c(min(c(data$actual_m, data$inferred_m)),
                                  max(c(data$actual_m, data$inferred_m)))) +
    scale_y_continuous(name = "Effective Max Trophic Level",
                       limits = c(min(c(data$actual_m, data$inferred_m)),
                                  max(c(data$actual_m, data$inferred_m)))) +
    scale_colour_viridis_d(name = "Food Web Type", option = "C", begin = 0.2, end = 0.8,
                           guide = guide_legend(override.aes = list(alpha = 1),
                                                reverse = TRUE, nrow = 1)) +
    theme_bw(base_size = 8) +
    theme(aspect.ratio = 1, legend.position = "bottom")
}


EE_empirical_comparison <- function(data, params, threatType, a, b, web) {
  data %>%
    filter(threat_type == threatType & alpha == a & beta == b & web_name == web) %>%
    ggplot(aes(x = threat_magnitude, y = excess_extinctions, group = threat_magnitude)) +
    geom_boxplot(data = ~subset(., replicate != 0), col = "#fca636", outlier.size = 0.5) +
    {
      overlay_lines_EE_e(filter(params, web_name == web & replicate == 0)$actual_m,
                         "black", "solid", 0.25, a, b)
    } +
    geom_point(data = ~subset(., replicate == 0), size = 0.5, col = "#6a00a8") +
    geom_text(x = Inf, y = Inf , hjust = 1.1, vjust = 1.5, 
              label = str_c("S = ",
                            filter(params, web_name == web)$S[1],
                            ", C = ",
                            round(filter(params, web_name == web)$C[1], 2)),
              size = 2.5, check_overlap = TRUE) +
    scale_x_continuous(name = "Direct Extinction Risk (\u03b5)",
                       limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
    scale_y_continuous(name = "Excess Extinction (E)") +
    labs(title = str_c(web)) +
    theme_bw(base_size = 8) +
    theme(plot.title = element_text(hjust = 0.5))
}


EE_empirical_comparison_full <- function(data, params, threatType, a, b, web) {
  data %>%
    filter(threat_type == threatType & alpha == a & beta == b & web_name == web) %>%
    ggplot(aes(x = threat_magnitude, y = excess_extinctions, group = threat_magnitude)) +
    geom_boxplot(data = ~subset(., replicate != 0), col = "#fca636", outlier.size = 0.5) +
    {
      overlay_lines_EE_e(filter(params, web_name == web & replicate == 0)$actual_m,
                         "black", "solid", 0.25, a, b)
    } +
    geom_point(data = ~subset(., replicate == 0), size = 0.5, col = "#6a00a8") +
    geom_text(x = Inf, y = Inf , hjust = 1.1, vjust = 1.5, 
              label = str_c("S = ",
                            filter(params, web_name == web)$S[1],
                            ", C = ",
                            round(filter(params, web_name == web)$C[1], 2)),
              size = 2.5, check_overlap = TRUE) +
    scale_x_continuous(name = "Direct Extinction Risk (\u03b5)",
                       limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
    scale_y_continuous(name = "Excess Extinction (E)") +
    labs(title = str_c(web)) +
    theme_bw(base_size = 8) +
    theme(plot.title = element_text(hjust = 0.5))
}


excess_loss_ES_by_web <- function(data, threatType, a, b) {
  data %>%
    filter(threat_type == threatType & alpha == a & beta == b) %>%
    ggplot(aes(x = threat_magnitude, y = excess_extinctions, color = ES_name)) +
    geom_line() +
    facet_wrap(~web_name) +
    scale_x_continuous(name = "Direct Extinction Risk (\u03b5)",
                       limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
    scale_y_continuous(name = "Excess Loss of ES") +
    scale_color_manual(values = ES_web_color_key, name = "Ecosystem\nService",
                       guide = guide_legend(reverse = TRUE), drop = F) +
    theme_bw(base_size = 8)
}


niche_S_TL <- function(data) {
  data %>%
    mutate(S = as.numeric(as.character(S))) %>%
    ggplot(aes(x = S, y = actual_m)) +
    geom_point(size = 1, alpha = 0.5, position = position_jitter(width = 2, seed = 63)) +
    stat_smooth(method = "lm") +
    scale_x_continuous(name = "Food Web Size") +
    scale_y_continuous(name = "Maximum TL of Food Web") +
    theme_bw()
}


niche_C_TL <- function(data) {
  data %>%
    mutate(C = as.numeric(as.character(C))) %>%
    ggplot(aes(x = C, y = actual_m)) +
    geom_point(size = 1, alpha = 0.5, position = position_jitter(width = 0.0025, seed = 63)) +
    stat_smooth(method = "lm") +
    scale_x_continuous(name = "Food Web Connectance") +
    scale_y_continuous(name = "Maximum TL of Food Web") +
    theme_bw()
}


omnivory_plot <- function(data) {
  data %>%
    ggplot(aes(x = omnivory, y = deviation)) +
    geom_point(size = 1, col = "#6a00a8") +
    annotate(geom = "text", x = -Inf, y = Inf, hjust = -0.05, vjust = 2,
             label = str_c("Correlation coefficient = ",
                           round(cor(data$omnivory, data$deviation), digits = 4))) +
    scale_x_continuous(name = "Omnivory Index",
                       limits = c(min(data$omnivory), max(data$omnivory))) +
    scale_y_continuous(name = "Deviation from Stratified Web Prediction",
                       limits = c(min(data$deviation), max(data$deviation))) +
    theme_bw(base_size = 8)
}


TL_empirical_comparison <- function(web, params) {
  params %>%
    filter(web_name == web) %>%
    ggplot(aes(x = actual_m)) +
    geom_histogram(data = ~subset(., replicate != 0), bins = 15) +
    geom_vline(xintercept = filter(params, replicate == 0 & web_name == web)$actual_m,
               col = "red", linetype = "dashed") +
    scale_x_continuous(name = "Maximum Trophic Level") +
    scale_y_continuous(name = "Count") +
    labs(title = str_c(web)) +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5))
}



# Import data -------------------------------------------------------------------------

strat_outcome <- read_rds("data/stratified_outcome_small.rds") %>%
  mutate(fragility = (indirect_extinctions - direct_extinctions) / (trophic_levels * species_per_level - direct_extinctions)) %>%
  mutate(across(c(web_ID, threat_type, trophic_levels), as_factor))
strat_param <- read_rds("data/stratified_param_fit.rds") %>%
  select(!threat_type & !alpha & !beta) %>%
  mutate(across(c(web_ID, trophic_levels), as_factor))


niche_outcome <- read_rds("data/niche_outcome_small.rds")  %>%
  mutate(fragility = (indirect_extinctions - direct_extinctions) / (S - direct_extinctions)) %>%
  mutate(across(c(web_ID, threat_type, S, C), as_factor))
niche_param <- read_rds("data/niche_param_fit.rds") %>%
  select(!threat_type & !alpha & !beta) %>%
  mutate(across(c(web_ID, S, C), as_factor))
niche_param_adjusted <- read_rds("data/niche_param_adjusted.rds") %>%
  mutate(across(c(S, C), as_factor))


empirical_outcome <- read_rds("data/empirical_outcome_small.rds") %>%
  mutate(web_name = case_match(
    web_ID %% 12,
    1 ~ "BSQ", 2 ~ "CHB", 3 ~ "CSM", 4 ~ "EMB", 5 ~ "EPB", 6 ~ "EWB",
    7 ~ "LRL", 8 ~ "MPC", 9 ~ "MRM", 10 ~ "PCR", 11 ~ "PRV", 0 ~ "STM")) %>%
  mutate(fragility = (indirect_extinctions - direct_extinctions) / (S - direct_extinctions)) %>%
  mutate(across(c(web_ID, threat_type, web_name), as_factor))
empirical_param <- read_rds("data/empirical_param_fit.rds") %>%
  select(!threat_type & !alpha & !beta) %>%
  mutate(web_name = case_match(web_ID %% 12,
                               1 ~ "BSQ", 2 ~ "CHB", 3 ~ "CSM", 4 ~ "EMB", 5 ~ "EPB",
                               6 ~ "EWB", 7 ~ "LRL", 8 ~ "MPC", 9 ~ "MRM", 10 ~ "PCR",
                               11 ~ "PRV", 0 ~ "STM")) %>%
  mutate(replicate = rep(0:50, each = 12)) %>%
  mutate(across(c(web_ID, web_name), as_factor))
empirical_web_color_key <- viridis::plasma(length(unique(empirical_outcome$web_name)),
                                           end = 0.9)
names(empirical_web_color_key) <- sort(unique(empirical_outcome$web_name))


empirical_omnivory <- read_rds("data/empirical_omnivory_deviation.rds") %>%
  mutate(web_name = case_match(web_ID %% 12,
                               1 ~ "BSQ", 2 ~ "CHB", 3 ~ "CSM", 4 ~ "EMB", 5 ~ "EPB",
                               6 ~ "EWB", 7 ~ "LRL", 8 ~ "MPC", 9 ~ "MRM", 10 ~ "PCR",
                               11 ~ "PRV", 0 ~ "STM")) %>%
  mutate(across(c(web_ID, web_name), as_factor))


ES_outcome <- read_rds("data/empirical_ES_outcome_small.rds") %>%
  mutate(web_name = case_match(web_ID %% 12,
                               1 ~ "BSQ", 2 ~ "CHB", 3 ~ "CSM", 4 ~ "EMB", 5 ~ "EPB",
                               6 ~ "EWB", 7 ~ "LRL", 8 ~ "MPC", 9 ~ "MRM", 10 ~ "PCR",
                               11 ~ "PRV", 0 ~ "STM")) %>%
  mutate(ES_name = paste0(toupper(substr(ES_name, 1, 1)),
                          substring(ES_name, 2))) %>%
  mutate(across(c(web_ID, threat_type, web_name, ES_name), as_factor))
ES_web_color_key <- viridis::plasma(length(unique(ES_outcome$ES_name)),
                                    end = 0.9)
names(ES_web_color_key) <- sort(unique(ES_outcome$ES_name))



# Create figures ----------------------------------------------------------------------

# figure 2
strat_list_EE <- tibble(threatType = c("deterministic")) %>%
  crossing(tibble(a = c(1, 1, 2), b = c(1, 2, 1))) %>%
  mutate(plot = pmap(., EE_strat_e, data = strat_outcome)) %>%
  pull(plot)
strat_list_F <- tibble(threatType = c("deterministic")) %>%
  crossing(tibble(a = c(1, 1, 2), b = c(1, 2, 1))) %>%
  mutate(plot = pmap(., F_strat_e, data = strat_outcome)) %>%
  pull(plot)
fig_2_list <- c(strat_list_EE, strat_list_F)
fig_2_bottom <- wrap_plots(fig_2_list) +
  plot_layout(axes = "collect", guides = "collect") +
  plot_annotation(tag_levels = 'a')
ggsave("figures_and_tables/figure_2_bottom.svg", fig_2_bottom, device = "svg",
       width = 173, height = 100, units = "mm")


# figure 3
fig_3_list <- tibble(threatType = c("deterministic")) %>%
  crossing(tibble(a = c(1), b = c(1))) %>%
  crossing(tibble(conn = c(0.05, 0.10, 0.15, 0.20))) %>%
  arrange(threatType, a, b, conn) %>%
  mutate(plot = pmap(., EE_niche_boxplot,
                     data = niche_outcome, params = niche_param_adjusted)) %>%
  pull(plot)
fig_3 <- wrap_plots(fig_3_list) +
  plot_layout(ncol = 2, axes = "collect", guides = "collect") +
  plot_annotation(tag_levels = 'a')
ggsave("figures_and_tables/figure_3.svg", fig_3, device = "svg",
       width = 173, height = 125, units = "mm")


# figure 4
niche_points <- niche_param %>%
  select(actual_m, inferred_m) %>%
  mutate(type = "Niche")
empirical_points <- empirical_param %>%
  filter(replicate == 0) %>%
  select(actual_m, inferred_m) %>%
  mutate(type = "Empirical")
empirical_comparison_list <- tibble(threatType = c("deterministic")) %>%
  crossing(tibble(a = c(1), b = c(1))) %>%
  crossing(web = c("BSQ", "CHB", "EPB", "EWB", "MRM", "PRV")) %>%
  arrange(web, threatType, a, b) %>%
  mutate(plot = pmap(., EE_empirical_comparison,
                     data = empirical_outcome, params = empirical_param)) %>%
  pull(plot)
fig_4_left <- rbind(empirical_points, niche_points) %>%
  correlation_plot()
fig_4_right <- wrap_plots(empirical_comparison_list) +
  plot_layout(ncol = 3, axes = "collect", guides = "collect") &
  ylim(1, 4.5) & ylab("Excess Extinction (E)")
ggsave("figures_and_tables/figure_4_left.svg", fig_4_left, device = "svg",
       width = 73, height = 73, units = "mm")
ggsave("figures_and_tables/figure_4_right.svg", fig_4_right, device = "svg",
       width = 100, height = 73, units = "mm")


# figure 5
fig_5 <- excess_loss_ES_by_web(ES_outcome, "deterministic", 1, 1)
ggsave("figures_and_tables/figure_5.svg", fig_5, device = "svg",
       width = 173, height = 100, units = "mm")


# figure S2
niche_list_EE_S_all <- tibble(threatType = c("deterministic")) %>%
  crossing(tibble(a = c(1, 1, 2), b = c(1, 2, 1))) %>%
  mutate(plot = pmap(., EE_niche_S, data = niche_outcome)) %>%
  pull(plot)
niche_list_F_S_all <- tibble(threatType = c("deterministic")) %>%
  crossing(tibble(a = c(1, 1, 2), b = c(1, 2, 1))) %>%
  mutate(plot = pmap(., F_niche_S, data = niche_outcome)) %>%
  pull(plot)
fig_S2_list <- c(niche_list_EE_S_all, niche_list_F_S_all)
fig_S2 <- wrap_plots(fig_S2_list) +
  plot_layout(axes = "collect", guides = "collect") +
  plot_annotation(tag_levels = 'a')
ggsave("figures_and_tables/figure_S2.svg", fig_S2, device = "svg",
       width = 173, height = 100, units = "mm")


# figure S3
niche_list_EE_C_all <- tibble(threatType = c("deterministic")) %>%
  crossing(tibble(a = c(1, 1, 2), b = c(1, 2, 1))) %>%
  mutate(plot = pmap(., EE_niche_C, data = niche_outcome)) %>%
  pull(plot)
niche_list_F_C_all <- tibble(threatType = c("deterministic")) %>%
  crossing(tibble(a = c(1, 1, 2), b = c(1, 2, 1))) %>%
  mutate(plot = pmap(., F_niche_C, data = niche_outcome)) %>%
  pull(plot)
fig_S3_list <- c(niche_list_EE_C_all, niche_list_F_C_all)
fig_S3 <- wrap_plots(fig_S3_list) +
  plot_layout(axes = "collect", guides = "collect") +
  plot_annotation(tag_levels = 'a')
ggsave("figures_and_tables/figure_S3.svg", fig_S3, device = "svg",
       width = 173, height = 100, units = "mm")


# figure S4
empirical_list_EE_S <- tibble(threatType = c("deterministic")) %>%
  crossing(tibble(a = c(1, 2), b = c(1, 1))) %>%
  mutate(plot = pmap(., EE_empirical_S, data = empirical_outcome)) %>%
  pull(plot)
empirical_list_F_S <- tibble(threatType = c("deterministic")) %>%
  crossing(tibble(a = c(1, 2), b = c(1, 1))) %>%
  mutate(plot = pmap(., F_empirical_S, data = empirical_outcome)) %>%
  pull(plot)
empirical_list_EE_C <- tibble(threatType = c("deterministic")) %>%
  crossing(tibble(a = c(1, 2), b = c(1, 1))) %>%
  mutate(plot = pmap(., EE_empirical_C, data = empirical_outcome)) %>%
  pull(plot)
empirical_list_F_C <- tibble(threatType = c("deterministic")) %>%
  crossing(tibble(a = c(1, 2), b = c(1, 1))) %>%
  mutate(plot = pmap(., F_empirical_C, data = empirical_outcome)) %>%
  pull(plot)
fig_S4_list <- c(empirical_list_EE_S, empirical_list_EE_C,
                 empirical_list_F_S, empirical_list_F_C)  
fig_S4 <- wrap_plots(fig_S4_list) +
  plot_layout(ncol = 4, axes = "collect", guides = "collect") +
  plot_annotation(tag_levels = 'a')
ggsave("figures_and_tables/figure_S4.svg", fig_S4, device = "svg",
       width = 173, height = 100, units = "mm")


# figure S5
empirical_list_EE_S_all <- tibble(threatType = c("deterministic")) %>%
  crossing(tibble(a = c(1, 1, 2), b = c(1, 2, 1))) %>%
  mutate(plot = pmap(., EE_empirical_S, data = empirical_outcome)) %>%
  pull(plot)
empirical_list_F_S_all <- tibble(threatType = c("deterministic")) %>%
  crossing(tibble(a = c(1, 1, 2), b = c(1, 2, 1))) %>%
  mutate(plot = pmap(., F_empirical_S, data = empirical_outcome)) %>%
  pull(plot)
fig_S5_list <- c(empirical_list_EE_S_all, empirical_list_F_S_all)
fig_S5 <- wrap_plots(fig_S5_list) +
  plot_layout(axes = "collect", guides = "collect") +
  plot_annotation(tag_levels = 'a')
ggsave("figures_and_tables/figure_S5.svg", fig_S5, device = "svg",
       width = 173, height = 100, units = "mm")


# figure S6
empirical_list_EE_C_all <- tibble(threatType = c("deterministic")) %>%
  crossing(tibble(a = c(1, 1, 2), b = c(1, 2, 1))) %>%
  mutate(plot = pmap(., EE_empirical_C, data = empirical_outcome)) %>%
  pull(plot)
empirical_list_F_C_all <- tibble(threatType = c("deterministic")) %>%
  crossing(tibble(a = c(1, 1, 2), b = c(1, 2, 1))) %>%
  mutate(plot = pmap(., F_empirical_C, data = empirical_outcome)) %>%
  pull(plot)
fig_S6_list <- c(empirical_list_EE_C_all, empirical_list_F_C_all)
fig_S6 <- wrap_plots(fig_S6_list) +
  plot_layout(axes = "collect", guides = "collect") +
  plot_annotation(tag_levels = 'a')
ggsave("figures_and_tables/figure_S6.svg", fig_S6, device = "svg",
       width = 173, height = 100, units = "mm")


# figure S7
fig_S7_left <- niche_S_TL(niche_param)
fig_S7_right <- niche_C_TL(niche_param)
fig_S7 <- fig_S7_left + fig_S7_right + plot_layout(axes = "collect")
ggsave("figures_and_tables/figure_S7.svg", fig_S7, device = "svg",
       width = 110, height = 66, units = "mm")


# figure S8
fig_S8 <- omnivory_plot(empirical_omnivory)
ggsave("figures_and_tables/figure_S8.svg", fig_S8, device = "svg",
       width = 82, height = 62, units = "mm")


# figure S9
fig_S9 <- tibble(web = c("BSQ", "CHB", "CSM", "EMB", "EPB", "EWB",
                         "LRL", "MPC", "MRM", "PCR", "PRV", "STM")) %>%
  mutate(plot = pmap(., TL_empirical_comparison, params = empirical_param)) %>%
  pull(plot) %>%
  wrap_plots() + plot_layout(axes = "collect") &
  ylim(0, 13) & ylab("Count") & xlim(3, 9) & xlab("Maximum Trophic Level")
ggsave("figures_and_tables/figure_S9.svg", fig_S9, device = "svg",
       width = 173, height = 100, units = "mm")


# figure S10
strat_list_EE_det_all <- tibble(threatType = c("deterministic")) %>%
  crossing(tibble(a = c(1, 1, 2), b = c(1, 2, 1))) %>%
  mutate(plot = pmap(., EE_strat_TL, data = strat_outcome)) %>%
  pull(plot)
strat_list_EE_sto_all <- tibble(threatType = c("stochastic")) %>%
  crossing(tibble(a = c(1, 1, 2), b = c(1, 2, 1))) %>%
  mutate(plot = pmap(., EE_strat_TL, data = strat_outcome)) %>%
  pull(plot)
strat_list_EE_basal_all <- tibble(threatType = c("basal-only")) %>%
  crossing(tibble(a = c(1, 1, 2), b = c(1, 2, 1))) %>%
  mutate(plot = pmap(., EE_strat_TL, data = strat_outcome)) %>%
  pull(plot)
fig_S10_list <- c(strat_list_EE_det_all, strat_list_EE_sto_all, strat_list_EE_basal_all)
fig_S10 <- wrap_plots(fig_S10_list) +
  plot_layout(axes = "collect", guides = "collect") +
  plot_annotation(tag_levels = 'a')
ggsave("figures_and_tables/figure_S10.svg", fig_S10, device = "svg",
       width = 173, height = 100, units = "mm")


# figure S11
strat_list_F_det_all <- tibble(threatType = c("deterministic")) %>%
  crossing(tibble(a = c(1, 1, 2), b = c(1, 2, 1))) %>%
  mutate(plot = pmap(., F_strat_TL, data = strat_outcome)) %>%
  pull(plot)
strat_list_F_sto_all <- tibble(threatType = c("stochastic")) %>%
  crossing(tibble(a = c(1, 1, 2), b = c(1, 2, 1))) %>%
  mutate(plot = pmap(., F_strat_TL, data = strat_outcome)) %>%
  pull(plot)
strat_list_F_basal_all <- tibble(threatType = c("basal-only")) %>%
  crossing(tibble(a = c(1, 1, 2), b = c(1, 2, 1))) %>%
  mutate(plot = pmap(., F_strat_TL, data = strat_outcome)) %>%
  pull(plot)
fig_S11_list <- c(strat_list_F_det_all, strat_list_F_sto_all, strat_list_F_basal_all)
fig_S11 <- wrap_plots(fig_S11_list) +
  plot_layout(axes = "collect", guides = "collect") +
  plot_annotation(tag_levels = 'a')
ggsave("figures_and_tables/figure_S11.svg", fig_S11, device = "svg",
       width = 173, height = 100, units = "mm")


# figure S12
niche_list_EE_det_all <- tibble(threatType = c("deterministic")) %>%
  crossing(tibble(a = c(1, 1, 2), b = c(1, 2, 1))) %>%
  mutate(plot = pmap(., EE_niche_TL, data = niche_outcome, params = niche_param)) %>%
  pull(plot)
niche_list_EE_sto_all <- tibble(threatType = c("stochastic")) %>%
  crossing(tibble(a = c(1, 1, 2), b = c(1, 2, 1))) %>%
  mutate(plot = pmap(., EE_niche_TL, data = niche_outcome, params = niche_param)) %>%
  pull(plot)
niche_list_EE_basal_all <- tibble(threatType = c("basal-only")) %>%
  crossing(tibble(a = c(1, 1, 2), b = c(1, 2, 1))) %>%
  mutate(plot = pmap(., EE_niche_TL, data = niche_outcome, params = niche_param)) %>%
  pull(plot)
fig_S12_list <- c(niche_list_EE_det_all, niche_list_EE_sto_all, niche_list_EE_basal_all)
fig_S12 <- wrap_plots(fig_S12_list) +
  plot_layout(axes = "collect", guides = "collect") +
  plot_annotation(tag_levels = 'a')
ggsave("figures_and_tables/figure_S12.svg", fig_S12, device = "svg",
       width = 173, height = 100, units = "mm")


# figure S13
niche_list_F_det_all <- tibble(threatType = c("deterministic")) %>%
  crossing(tibble(a = c(1, 1, 2), b = c(1, 2, 1))) %>%
  mutate(plot = pmap(., F_niche_TL, data = niche_outcome, params = niche_param)) %>%
  pull(plot)
niche_list_F_sto_all <- tibble(threatType = c("stochastic")) %>%
  crossing(tibble(a = c(1, 1, 2), b = c(1, 2, 1))) %>%
  mutate(plot = pmap(., F_niche_TL, data = niche_outcome, params = niche_param)) %>%
  pull(plot)
niche_list_F_basal_all <- tibble(threatType = c("basal-only")) %>%
  crossing(tibble(a = c(1, 1, 2), b = c(1, 2, 1))) %>%
  mutate(plot = pmap(., F_niche_TL, data = niche_outcome, params = niche_param)) %>%
  pull(plot)
fig_S13_list <- c(niche_list_F_det_all, niche_list_F_sto_all, niche_list_F_basal_all)
fig_S13 <- wrap_plots(fig_S13_list) +
  plot_layout(axes = "collect", guides = "collect") +
  plot_annotation(tag_levels = 'a')
ggsave("figures_and_tables/figure_S13.svg", fig_S13, device = "svg",
       width = 173, height = 100, units = "mm")


# figure S14
empirical_list_EE_det_all <- tibble(threatType = c("deterministic")) %>%
  crossing(tibble(a = c(1, 1, 2), b = c(1, 2, 1))) %>%
  mutate(plot = pmap(., EE_empirical_TL,
                     data = empirical_outcome, params = empirical_param)) %>%
  pull(plot)
empirical_list_EE_sto_all <- tibble(threatType = c("stochastic")) %>%
  crossing(tibble(a = c(1, 1, 2), b = c(1, 2, 1))) %>%
  mutate(plot = pmap(., EE_empirical_TL,
                     data = empirical_outcome, params = empirical_param)) %>%
  pull(plot)
empirical_list_EE_basal_all <- tibble(threatType = c("basal-only")) %>%
  crossing(tibble(a = c(1, 1, 2), b = c(1, 2, 1))) %>%
  mutate(plot = pmap(., EE_empirical_TL,
                     data = empirical_outcome, params = empirical_param)) %>%
  pull(plot)
fig_S14_list <- c(empirical_list_EE_det_all,
                  empirical_list_EE_sto_all,
                  empirical_list_EE_basal_all)
fig_S14 <- wrap_plots(fig_S14_list) +
  plot_layout(axes = "collect", guides = "collect") +
  plot_annotation(tag_levels = 'a')
ggsave("figures_and_tables/figure_S14.svg", fig_S14, device = "svg",
       width = 173, height = 100, units = "mm")


# figure S15
empirical_list_F_det_all <- tibble(threatType = c("deterministic")) %>%
  crossing(tibble(a = c(1, 1, 2), b = c(1, 2, 1))) %>%
  mutate(plot = pmap(., F_empirical_TL,
                     data = empirical_outcome, params = empirical_param)) %>%
  pull(plot)
empirical_list_F_sto_all <- tibble(threatType = c("stochastic")) %>%
  crossing(tibble(a = c(1, 1, 2), b = c(1, 2, 1))) %>%
  mutate(plot = pmap(., F_empirical_TL,
                     data = empirical_outcome, params = empirical_param)) %>%
  pull(plot)
empirical_list_F_basal_all <- tibble(threatType = c("basal-only")) %>%
  crossing(tibble(a = c(1, 1, 2), b = c(1, 2, 1))) %>%
  mutate(plot = pmap(., F_empirical_TL,
                     data = empirical_outcome, params = empirical_param)) %>%
  pull(plot)
fig_S15_list <- c(empirical_list_F_det_all,
                  empirical_list_F_sto_all,
                  empirical_list_F_basal_all)
fig_S15 <- wrap_plots(fig_S15_list) +
  plot_layout(axes = "collect", guides = "collect") +
  plot_annotation(tag_levels = 'a')
ggsave("figures_and_tables/figure_S15.svg", fig_S15, device = "svg",
       width = 173, height = 100, units = "mm")


# figure S16
full_empirical_comparison_list <- tibble(threatType = c("deterministic")) %>%
  crossing(tibble(a = c(1), b = c(1))) %>%
  crossing(web = c("BSQ", "CHB", "CSM", "EMB", "EPB", "EWB",
                   "LRL", "MPC", "MRM", "PCR", "PRV", "STM")) %>%
  arrange(web, threatType, a, b) %>%
  mutate(plot = pmap(., EE_empirical_comparison_full,
                     data = empirical_outcome, params = empirical_param)) %>%
  pull(plot)
fig_S16 <- wrap_plots(full_empirical_comparison_list) +
  plot_layout(ncol = 4, axes = "collect", guides = "collect") &
  ylim(1, 4.5) & ylab("Excess Extinction (E)")
ggsave("figures_and_tables/figure_S16.svg", fig_S16, device = "svg",
       width = 173, height = 100, units = "mm")