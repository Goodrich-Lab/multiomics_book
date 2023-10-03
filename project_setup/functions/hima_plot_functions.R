## ---- HIMA_early_plot ----
#' Plot of High Dimensional Mediation Analysis with Early Integration
#' 
#' Given a tidy dataframe summarizing the results of HIMA analysis 
#' and return a plot of the results
#'
#' @param result_hima_early  a tidy dataframe summarizing the results of HIMA analysis
#'
#' @return a figure of the result of hima early intergration 
#'
#' @import dplyr
#' @importFrom ggplot2 ggplot
#' 
hima_early_integration_plot <- function(result_hima_early) {
  
  # Pivot longer for figure
  early_long <- result_hima_early %>% 
    pivot_longer(cols = c(`TE (%)`, alpha, beta), 
                 names_to = "name") %>%
    mutate(name = factor(name, levels = c("alpha", "beta", "TE (%)")))
  
  # Plot features
  fig <- ggplot(early_long, 
         aes(x = ftr_name, 
             y = value,
             fill = omic_layer)) + 
    geom_bar(stat = "identity") +
    facet_grid(name ~ omic_layer, 
               scales = "free",
               space = "free_x") +
    scale_fill_brewer(type = "qual", palette = 2) +
    geom_hline(yintercept = 0, linetype = 1, color = "grey50") + 
    ylab(NULL) + xlab(NULL) +
    theme(
      strip.background = element_blank(),
      strip.text.x = element_blank(),
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5),
      legend.title = element_blank(), 
      legend.position = "bottom", 
      legend.justification = c(1, 0))
  
  return(fig)
}


## ---- HIMA_intermediate_plot ----

#' Plot of High Dimensional Mediation Analysis with Intermediate Integration
#' 
#' Given a tidy dataframe summarizing the results of HIMA analysis
#'
#' @param result_hima_intermediate  a tidy dataframe summarizing the results of HIMA analysis
#'
#' @return a figure of the result of hima intermediate intergration 
#'
#' @import dplyr
#' @importFrom ggplot2 ggplot
#' 
hima_intermediate_integration_plot <- function(result_hima_intermediate) {
  # Pivot longer for figure
  intermediate_long <- result_hima_intermediate %>% 
    pivot_longer(cols = c(alpha, beta,`TE (%)`), 
                 names_to = "name") %>%
    mutate(name = factor(name, levels = c("alpha", "beta", "TE (%)")))
  
  # Plot features
  p <- ggplot(intermediate_long, 
         aes(x = fct_inorder(ftr_name), 
             y = value,
             fill = omic_layer)) + 
    geom_bar(stat = "identity") +
    facet_grid(name ~ omic_layer, 
               scales = "free",
               space = "free_x") +
    scale_fill_brewer(type = "qual", palette = 2) +
    geom_hline(yintercept = 0, linetype = 1, color = "grey50") + 
    ylab(NULL) + xlab(NULL) +
    theme(
      strip.background = element_blank(),
      strip.text.x = element_blank(),
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5),
      legend.title = element_blank(), 
      legend.position = "bottom", # Place the legend at the bottom
      legend.justification = c(1, 0))
  
  return(p)
}


## ---- HIMA_late_plot ----
#' Plot of High Dimensional Mediation Analysis with Intermediate Integration
#' 
#' Given a tidy dataframe summarizing the results of HIMA analysis
#'
#' @param result_hima_late  a tidy dataframe summarizing the results of HIMA analysis
#'
#' @return a figure of the result of hima intermediate intergration 
#'
#' @import dplyr
#' @importFrom ggplot2 ggplot
#' 
hima_late_integration_plot <- function(result_hima_late) {
  
  # Pivot longer for figure
  late_long <- result_hima_late %>% 
    pivot_longer(cols = c(`TE (%)`, beta, alpha), 
                 names_to = "name") %>%
    mutate(name = factor(name, levels = c("alpha", "beta", "TE (%)")))
  
  # Plot features
  p <- ggplot(late_long, 
         aes(x = ftr_name, 
             y = value,
             fill = omic_layer)) + 
    geom_bar(stat = "identity") +
    facet_grid(name ~ omic_layer, 
               scales = "free",
               space = "free_x") +
    scale_fill_brewer(type = "qual", palette = 2) +
    geom_hline(yintercept = 0, linetype = 1, color = "grey50") + 
    ylab(NULL) + xlab(NULL) +
    theme(
      strip.background = element_blank(),
      strip.text.x = element_blank(),
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5),
      legend.title = element_blank(), 
      legend.position = "bottom", 
      legend.justification = c(1, 0))
  
  return(p)
}