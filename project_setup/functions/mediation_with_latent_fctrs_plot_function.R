## ---- Mediation with latent fctrs early ----
#' Plot result of mediation with latent factors in early integration
#' 
#' Given a list including two tidy dataframes summarizing the results of HIMA analysis,
#' function runs plotting and returns a plot
#'
#' @param med_early_list A list including two tidy dataframes summarizing the results of HIMA analysis

#'
#' @return a figure of the result mediation with latent factors in intermediate integration
#' 
#' @import dplyr
#' @importFrom dplyr left_join
#' @importFrom ggplot2 ggplot
#' @importFrom base apply
#' 

med_with_latent_fctrs_early_integration_plot <- function(med_early_list) {
  
  # Panel A Bargraph of mediation effects of PC's --------------
  pca_early_long <- med_early_list$result_hima_pca_early_sig%>%
    pivot_longer(cols = c(Alpha, Beta, `TE (%)`))
  
  # Plot 
  panel_a_early <- ggplot(pca_early_long, aes(x = pcs_ordered, y = value)) +
    geom_bar(stat = "identity", fill = "grey50") + 
    geom_hline(yintercept = 0) +
    facet_grid(name ~ ., scales = "free", space = "free_x", switch = "y") + 
    ggh4x::facetted_pos_scales(
      y = list(name == "Alpha"  ~ scale_y_continuous(
        limits = c(-.45,.45), breaks = c(-0.4, 0, .4)),
        name == "Beta"   ~ scale_y_continuous(
          limits = c(-.45,.45), breaks = c(-0.4, 0, .4)),
        name == "TE (%)" ~ scale_y_continuous(
          limits = c(-1,55), n.breaks = 4))) + 
    theme(axis.title = element_blank(), 
          strip.placement = "outside",
          axis.text.x = element_blank(),
          strip.background = element_blank(),
          text = element_text(size = 15))
  
  # Panel B: heatmap of correlation of features vs PC's --------------
  
  # Select features for plots -------
  ## Select rows with values in the top 10 of their respective columns 
  top <- apply(med_early_list$result_ftr_cor_sig_pcs_early %>% select(-omic_layer, -omic_num) %>% column_to_rownames("feature"), 2, function(x) x %in% tail(sort(abs(x)), 10)) 
  
  ## Filter selected rows
  ftr_cor_sig_pcs_top <- med_early_list$result_ftr_cor_sig_pcs_early[which(rowSums(top) > 0), ] 
  
  # Pivot longer
  ftr_early_pcs_long_top_ft <-
    ftr_cor_sig_pcs_top %>%
    pivot_longer(cols = all_of(med_early_list$result_hima_pca_early_sig$pc_num),
                 names_to = "pc_num",
                 values_to = "Correlation")
  
  # Join with PCA long to get the ordered pcs
  panel_b_dat_early_top_ft <- left_join(ftr_early_pcs_long_top_ft,
                                        pca_early_long %>%
                                          filter(name == "Alpha") %>%
                                          dplyr::select(pc_num, pcs_ordered))
  
  # Plot
  panel_b_early <- ggplot(data = panel_b_dat_early_top_ft,
                          aes(y = feature,
                              x = fct_rev(pcs_ordered), 
                              fill = Correlation)) +
    geom_tile(color = "white") +
    facet_grid(omic_layer ~ ., scales = "free", space = "free") +
    scale_fill_gradient2(low  = "blue",
                         mid  = "white",
                         high = "red",
                         midpoint = 0,
                         limits = c(-1, 1),
                         breaks = c(-1, 0, 1),
                         na.value = "grey20") +
    theme(
      axis.text.x = element_text(size = 15,angle = 90, hjust = 1, vjust = .5),
      strip.text = element_blank(),
      axis.title = element_blank(), 
      axis.ticks.x = element_blank(),
      legend.position = "none",
      text = element_text(size = 20)) 
  
  # Combine Figures 
  p <- cowplot::plot_grid(
    NULL, panel_a_early,  NULL, panel_b_early, 
    ncol = 1, align = "v", axis = "lr",
    rel_heights  = c(.05, .75, .1, 1.75),
    labels = c("a)","", "b) "))

  return(p)
}

## ---- Mediation with latent fctrs intermediate plot ----
#' Plot result of mediation with latent factors in intermediate integration
#' 
#' Given a list including two tidy dataframes summarizing the results of HIMA analysis,
#' function runs plotting and returns a plot
#'
#' @param med_int_list A list including two tidy dataframes summarizing the results of HIMA analysis

#'
#' @return a figure of the result mediation with latent factors in early integration
#' 
#' @import dplyr
#' @importFrom dplyr left_join
#' @importFrom ggplot2 ggplot
#' @importFrom base apply
#'
med_with_latent_fctrs_intermediate_integration_plot <- function(med_int_list) {
  # Panel A: Bargraph of mediation effects of PC's --------------
  result_hima_jive_long <- med_int_list$result_hima_jive_sig %>%
    pivot_longer(cols = c(Alpha, Beta, `TE (%)`))
  
  # Plot  
  panel_a_jive <- ggplot(result_hima_jive_long, 
                         aes(x = components_ord, 
                             y = value)) +
    geom_bar(stat = "identity", fill = "grey50") + 
    geom_hline(yintercept = 0) +
    facet_grid(name ~ .,
               scales = "free",
               space = "free_x", switch = "y") + 
    # labs(y = NULL) +
    facetted_pos_scales(
      y = list(
        name == "Alpha"  ~ scale_y_continuous(
          limits = c(-.45,.45), breaks = c(-0.4, 0, .4)),
        name == "Beta"   ~ scale_y_continuous(
          limits = c(-.45,.45), breaks = c(-0.4, 0, .4)),
        name == "TE (%)" ~ scale_y_continuous(
          limits = c(-1,55), n.breaks = 4))) + 
    theme(axis.title = element_blank(), 
          title = element_text(hjust = 0),
          strip.placement = "outside",
          axis.text.y = element_text(size = 15),
          axis.text.x = element_blank(),
          strip.background = element_blank(),
          text = element_text(size = 15))
  
  # Panel B: heatmap of correlation of features vs PC's --------------
  
  # Select features for plots -------
  # Select rows with values in the top 10 of their respective columns
  top10 <- apply(med_int_list$result_ftr_cor_sig_pcs_jive%>% select(-omic_layer, -omic_num) %>% column_to_rownames("feature"), 2, function(x) x %in% tail(sort(abs(x)), 10))
  
  ftr_cor_sig_pcs_top <- med_int_list$result_ftr_cor_sig_pcs_jive[which(rowSums(top10) > 0), ]
  
  # Pivot longer
  ftr_jive_pcs_long_top_ft <-
    ftr_cor_sig_pcs_top %>%
    pivot_longer(cols = all_of(med_int_list$result_hima_jive_sig$component),
                 names_to = "component",
                 values_to = "Correlation")
  
  # Change features which are in wrong JIVE individual component to zero
  panel_b_dat_jive_fin <- ftr_jive_pcs_long_top_ft %>%
    mutate(
      in_ind_omic = str_sub(component, 1, 5) == str_sub(omic_layer, 1, 5),
      Correlation = ifelse(in_ind_omic | str_detect(component, "Joint"),
                           Correlation, NA)) %>%
    left_join(result_hima_jive_long %>% 
                filter(name == "Alpha") %>%
                dplyr::select(component,components_ord, ind_joint_num),
              by= "component",
              relationship = "many-to-many") %>%
    mutate(component = str_replace_all(component, "_", " ") %>% toTitleCase())
  
  # Plot
  panel_b_jive <- ggplot(data = panel_b_dat_jive_fin,
                         aes(y = feature,
                             x = components_ord, 
                             fill = Correlation)) +
    geom_tile(color = "white") +
    facet_grid(omic_layer ~ ., scales = "free", space = "free") +
    scale_fill_gradient2(low  = "blue",
                         mid  = "white",
                         high = "red",
                         midpoint = 0,
                         limits = c(-1, 1),
                         breaks = c(-1, 0, 1),
                         na.value = "grey") +
    theme(
      axis.text.x = element_text(size = 15,angle = 90, hjust = 1, vjust = .5),
      strip.text = element_blank(),
      axis.title = element_blank(), 
      axis.ticks.x = element_blank(),
      legend.position = "none",
      text = element_text(size = 20)) 
  
  # Combine Figures 
  p <- cowplot::plot_grid(
    NULL, panel_a_jive,  NULL, panel_b_jive, 
    ncol = 1, align = "v", axis = "lr",
    rel_heights  = c(.05, .75, .1, 1.75),
    labels = c("a)","", "b) ")) 
  
  return(p)
}


## ---- Mediation with latent fctrs late plot ----
#' Plot result of mediation with latent factors in late integration
#' 
#' Given a list including two tidy dataframes summarizing the results of HIMA analysis,
#' function runs plotting and returns a plot
#'
#' @param med_late_list A list including two tidy dataframes summarizing the results of HIMA analysis

#'
#' @return a figure of the result mediation with latent factors in late integration
#' 
#' @import dplyr
#' @importFrom dplyr left_join
#' @importFrom ggplot2 ggplot
#' @importFrom base apply
#'
med_with_latent_fctrs_late_integration_plot <- function(med_late_list) {
  # i. Data: Panel A Bargraph of mediation effects of PC's --------------
  pca_late_long <- med_late_list$result_hima_late_sig %>%
    pivot_longer(cols = c(Alpha, Beta, `TE (%)`))
  
  # Plot 
  panel_a_late <- ggplot(pca_late_long, aes(x = pcs_ordered, y = value)) +
    geom_bar(stat = "identity", fill = "grey50") + 
    geom_hline(yintercept = 0) +
    facet_grid(name ~ omic_num, scales = "free", space = "free_x", switch = "y") + 
    ggh4x::facetted_pos_scales(
      y = list(name == "Alpha"  ~ scale_y_continuous(
        limits = c(-.45,.45), breaks = c(-0.4, 0, .4)),
        name == "Beta"   ~ scale_y_continuous(
          limits = c(-.45,.45), breaks = c(-0.4, 0, .4)),
        name == "TE (%)" ~ scale_y_continuous(
          limits = c(-1,55), n.breaks = 4))) + 
    theme(axis.title = element_blank(), 
          strip.placement = "outside",
          strip.text = element_blank(),
          axis.text.x = element_blank(),
          strip.background = element_blank(),
          text = element_text(size = 15))
  
  # ii. Data: Panel B: heatmap of correlation of features vs PC's-----
  # Select features for plots -------
  # Select rows with values in the top 10 of their respective columns
  top_late_int <- apply(med_late_list$result_ftr_cor_sig_pcs_late %>% select(-feature,-omic_layer, -omic_num) %>% janitor::remove_empty(which = "rows"), 2,
                        function(x) x %in% tail(sort(abs(x)), 10))
  selected_rows <- which(rowSums(top_late_int, na.rm = TRUE) > 0)
  
  ftr_cor_sig_pcs_late <- med_late_list$result_ftr_cor_sig_pcs_late[selected_rows, ]
  
  # Pivot longer
  ftr_late_pcs_long_top_ft  <-  pivot_longer(ftr_cor_sig_pcs_late, 
                                             cols = all_of(med_late_list$result_hima_late_sig$pc_num), 
                                             names_to = "pc_num", 
                                             values_to = "Correlation") 
  
  # Join with PCA long to get the ordered pcs 
  panel_b_dat_late_top <- left_join(ftr_late_pcs_long_top_ft, 
                                    pca_late_long %>%
                                      filter(name == "Alpha") %>%
                                      dplyr::select(pc_num, pcs_ordered), 
                                    by = "pc_num") %>%
    mutate(omic_num2 = case_when(str_detect(pc_num, "meth") ~ 1, 
                                 str_detect(pc_num, "transc") ~ 2, 
                                 str_detect(pc_num, "miR") ~ 3,
                                 str_detect(pc_num, "pro") ~ 4, 
                                 str_detect(pc_num, "met") ~ 5))
  
  # Plot
  panel_b_late <- ggplot(data = panel_b_dat_late_top,
                         aes(y = feature,
                             x = pcs_ordered, 
                             fill = Correlation)) + 
    geom_tile(color = "white") +
    facet_grid(omic_num ~ omic_num2, scales = "free", space = "free") +
    scale_fill_gradient2(low  = "blue",
                         mid  = "white",
                         high = "red",
                         midpoint = 0,
                         limits = c(-1, 1),
                         breaks = c(-1, 0, 1),
                         na.value = "grey") +
    theme(
      axis.text.x = element_text(size = 15,angle = 90, hjust = 1, vjust = .5),
      strip.text = element_blank(),
      axis.title = element_blank(), 
      axis.ticks.x = element_blank(),
      legend.position = "none",
      text = element_text(size = 20)) 
  
  # Combine Figures 
  p <- cowplot::plot_grid(
    NULL, panel_a_late,  NULL, panel_b_late, 
    ncol = 1, align = "v", axis = "lr",
    rel_heights  = c(.05, .75, .1, 1.75),
    labels = c("a)","", "b) "))
  
  return(p)
  
}
