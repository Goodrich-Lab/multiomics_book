## ---- plot_Omics_Profiles ----
#' Plot of Omics profiles for each cluster using LUCID
#' 
#' Given an object of class from LUCID
#'
#' @param fit an object of class from LUCID
#' @param integration_type type of integreation,, "Early" or "Intermediate"
#'
#' @return a figure of Omics profiles for each cluster using LUCID
#'
#' @import dplyr
#' @importFrom ggplot2 ggplot


plot_omics_profiles <- function(fit, integration_type) {
  if(integration_type == "Early"){
    M_mean = as.data.frame(fit$pars$mu)
    M_mean$cluster = as.factor(1:2)
    # Reshape the data
    M_mean_melt <- M_mean %>% 
      pivot_longer(cols = -cluster, names_to = "variable", values_to = "value")
    
    M_mean_melt <- M_mean_melt %>% 
      mutate(cluster = ifelse(cluster == 2, "High Risk", "Low Risk"))
    # add color label for omics layer
    M_mean_melt = M_mean_melt %>%
      mutate(color_label = case_when(str_detect(variable,  "cg") ~ "1", 
                                     str_detect(variable, "TC") ~ "2", 
                                     TRUE ~ "3"))
    
    fig <- ggplot(M_mean_melt, 
                  aes(fill = color_label, y = value, x = variable)) +
      geom_bar(position="dodge", stat="identity") +
      ggtitle("Omics profiles for 2 latent clusters") +
      # facet_wrap(~(cluster)) +
      facet_grid(rows = vars(cluster), scales = "free_y") +
      theme(legend.position="none") +
      geom_hline(yintercept = 0) +
      xlab("") +
      theme(text = element_text(size=10),
            axis.text.x = element_text(angle = 90, vjust = 1,
                                       hjust = 1),
            plot.margin = margin(10, 10, 10, 80),
            panel.background = element_rect(fill="white"), 
            strip.background = element_rect(fill = "white"),
            axis.line.x = element_line(color = "black"),
            axis.line.y = element_line(color = "black"),) +
      scale_fill_manual(values = c("#2fa4da", "#A77E69", "#e7b6c1"))
  }else if(integration_type == "Intermediate"){
    M_mean = as_tibble(fit$res_Mu_Sigma$Mu[[1]], rownames = "variable") %>%
      bind_rows(as_tibble(fit$res_Mu_Sigma$Mu[[2]], rownames = "variable")) %>%
      bind_rows(as_tibble(fit$res_Mu_Sigma$Mu[[3]], rownames = "variable"))
    
    # Reorder results because mirna order is reversed
    M_mean1 <- M_mean %>% 
      left_join(meta_df, by = c("variable" = "ftr_name")) %>%
      mutate(`Low Risk`  =  if_else(omic_layer == "miRna", V2, V1), 
             `High Risk` =  if_else(omic_layer == "miRna", V1, V2)) %>%
      dplyr::select(-c("V1", "V2"))
    
    # Pivot longer for figure 
    M_mean_l <- M_mean1 %>% 
      pivot_longer(cols = c(`Low Risk`, `High Risk`),
                   names_to = "cluster",
                   values_to = "value")
    
    # add color label for omics layer
    M_mean2 = M_mean_l %>%
      mutate(color_label = case_when(omic_layer == "methylome" ~ "1", 
                                     omic_layer == "transcriptome" ~ "2", 
                                     omic_layer == "miRna" ~ "3"), 
             low_high = if_else(str_detect(cluster, "Low"), 0,1),
             omic = if_else(omic_layer == "miRna", 
                            "miR",
                            str_sub(omic_layer, end = 1) %>% toupper()),
             omic_cluster = str_c(omic, low_high))
    
    # Filter only the top ## differential expressed features 
    M_mean2_top <- M_mean2 %>% 
      group_by(variable) %>% 
      filter(abs(value) == max(abs(value))) %>% 
      ungroup() %>% 
      arrange(max(abs(value))) %>% 
      group_by(omic_layer) %>% 
      slice_head(n=12) %>%
      ungroup()
    
    # Plots top 12 features
    fig <- ggplot(M_mean2  %>% filter(variable %in% M_mean2_top$variable),
           aes(fill = color_label, y = value, x = variable)) +
      geom_bar(position="dodge", stat="identity") +
      ggtitle("Omics profiles for 2 latent clusters - Lucid in Parallel") +
      facet_grid(rows = vars(cluster),
                 cols = vars(omic_layer), scales = "free_x", space = "free") +
      theme(legend.position="none") +
      geom_hline(yintercept = 0) +
      xlab("") +
      theme(text = element_text(size=10),
            axis.text.x = element_text(angle = 90, vjust = 1,
                                       hjust = 1),
            plot.margin = margin(10, 10, 10, 80),
            panel.background = element_rect(fill="white"),
            strip.background = element_rect(fill = "white"),
            axis.line.x = element_line(color = "black"),
            axis.line.y = element_line(color = "black"),) +
      scale_fill_manual(values = c("#2fa4da", "#A77E69", "#e7b6c1"))
  }
  
  return(fig)
}