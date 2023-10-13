## ---- plot_LUCID_Early ----
#' Plot Sankey Diagram for LUCID in Early integration
#' 
#' Given an object of class from LUCID
#'
#' @param lucid_fit1  an object of class from LUCID
#' @param text_size  size of the text in sankey diagram
#'
#' @return a Sankey Diagram for LUCID in Early integration
#'
#' @import dplyr
#' @importFrom ggplot2 ggplot



sankey_early_integration <- function(lucid_fit1, text_size = 15) {
  # Get sankey dataframe ----
  get_sankey_df <- function(x,
                            G_color = "dimgray", 
                            X_color = "#eb8c30",
                            Z_color = "#2fa4da", 
                            Y_color = "#afa58e", 
                            pos_link_color = "#67928b", 
                            neg_link_color = "#d1e5eb", 
                            fontsize = 10) {
    K <- x$K
    var.names <- x$var.names
    pars <- x$pars
    dimG <- length(var.names$Gnames)
    dimZ <- length(var.names$Znames)
    valueGtoX <- as.vector(t(x$pars$beta[, -1]))
    valueXtoZ <- as.vector(t(x$pars$mu))
    valueXtoY <- as.vector(x$pars$gamma$beta)[1:K]
    
    # GtoX
    GtoX <- data.frame(
      source = rep(x$var.names$Gnames, K), 
      target = paste0("Latent Cluster", 
                      as.vector(sapply(1:K, function(x) rep(x, dimG)))), 
      value = abs(valueGtoX), 
      group = as.factor(valueGtoX > 0))
    
    # XtoZ
    XtoZ <- data.frame(
      source = paste0("Latent Cluster", 
                      as.vector(sapply(1:K, 
                                       function(x) rep(x, dimZ)))), 
      target = rep(var.names$Znames, 
                   K), value = abs(valueXtoZ),
      group = as.factor(valueXtoZ > 
                          0))
    
    # subset top 25% of each omics layer
    top25<- XtoZ %>%
      filter(source == "Latent Cluster1") %>%
      mutate(omics = case_when(grepl("cg", target) ~ "Methylation",
                               grepl("tc", target) ~ "Transcriptome",
                               grepl("miR", target) ~ "miRNA")) %>%
      group_by(omics) %>%
      arrange(desc(value)) %>%
      slice(1:7) %>%
      ungroup()
    
    XtoZ_sub<- XtoZ %>%
      filter(target %in% top25$target)
    
    
    # XtoY
    XtoY <- data.frame(source = paste0("Latent Cluster", 1:K), 
                       target = rep(var.names$Ynames, K), value = abs(valueXtoY), 
                       group = as.factor(valueXtoY > 0))
    links <- rbind(GtoX, XtoZ_sub, XtoY)
    # links <- rbind(GtoX, XtoZ, XtoY)
    
    nodes <- data.frame(
      name = unique(c(as.character(links$source), 
                      as.character(links$target))), 
      group = as.factor(c(rep("exposure",
                              dimG), rep("lc", K), rep("biomarker", nrow(XtoZ_sub)/2), "outcome")))
    # group = as.factor(c(rep("exposure", 
    # dimG), rep("lc", K), rep("biomarker", dimZ), "outcome")))
    ## the following two lines were used to exclude covars from the plot
    links <- links %>% filter(!grepl("cohort", source) & 
                                !grepl("age", source) & 
                                !grepl("fish", source) &
                                !grepl("sex", source))
    nodes <- nodes %>% filter(!grepl("cohort", name) &
                                !grepl("age", name) & 
                                !grepl("fish", name) &
                                !grepl("sex", name)) 
    
    links$IDsource <- match(links$source, nodes$name) - 1
    links$IDtarget <- match(links$target, nodes$name) - 1
    
    color_scale <- data.frame(
      domain = c("exposure", "lc", "biomarker", 
                 "outcome", "TRUE", "FALSE"), 
      range = c(G_color, X_color, 
                Z_color, Y_color, pos_link_color, neg_link_color))
    
    sankey_df = list(links = links, 
                     nodes = nodes)
    return(sankey_df)
  }
  # 1. Get sankey dataframes ----
  sankey_dat <- get_sankey_df(lucid_fit1)
  n_omics <- length(lucid_fit1$var.names$Znames)
  # link data
  links <- sankey_dat[["links"]] 
  # node data
  nodes <- sankey_dat[["nodes"]] 
  
  nodes1 <- nodes %>% 
    mutate(group = case_when(str_detect(name,"Cluster") ~ "lc",
                             str_detect(name, "cg") ~ "CpG",
                             str_detect(name, "outcome") ~ "outcome",
                             str_detect(name, "pro") ~ "Prot",
                             str_detect(name, "met") ~ "Met",
                             str_detect(name, "tc") ~ "TC",
                             str_detect(name, "miR") ~ "miRNA",
                             str_detect(name, "G1") ~ "exposure"),
           name = ifelse(name == "G1", "Hg",name))
  links1 <- links %>%
    mutate(source = ifelse(source == "G1", "Hg",source))
  # 6. Plotly Version ----
  
  ## 6.1 Set Node Color Scheme: ----
  color_pal_sankey <- matrix(
    c("exposure", sankey_colors$range[sankey_colors$domain == "exposure"],
      "lc",       "#b3d8ff",
      "CpG",     sankey_colors$range[sankey_colors$domain == "layer1"],
      "TC",      sankey_colors$range[sankey_colors$domain == "layer2"],
      "miRNA", sankey_colors$range[sankey_colors$domain == "layer3"],
      "outcome",  sankey_colors$range[sankey_colors$domain == "Outcome"]), 
    ncol = 2, byrow = TRUE) %>%
    as_tibble(.name_repair = "unique") %>% 
    janitor::clean_names() %>%
    dplyr::rename(group = x1, color = x2)
  
  # Add color scheme to nodes
  nodes_new_plotly <- nodes1 %>% 
    left_join(color_pal_sankey) %>%
    mutate(
      x = case_when(
        group == "exposure" ~ 0,
        str_detect(name, "Cluster") ~ 1/3,
        str_detect(name, "cg")|
          str_detect(name, "tc")|
          str_detect(name, "miR")|
          str_detect(name, "outcome")~ 2/3
      ))
  
  nodes_new_plotly1 <- nodes_new_plotly %>%
    # Modify names of features for plotting
   dplyr::select(group, color, x, name)%>% 
    mutate(name = case_when(name == "value" ~ "<b>Hg</b>",
                            name == "Latent Cluster1" ~ "<b>Joint Omics\nProfile 0</b>",
                            name == "Latent Cluster2" ~ "<b>Joint Omics\nProfile 1</b>",
                            TRUE ~ name))
    
  
  ## 6.2 Get links for Plotly, set color ----
  links_new <- links1  %>%
    mutate(
      link_color = case_when(
        # Ref link color
        value == 0 ~     "#f3f6f4",
        # # Cluster 
        # str_detect(source, "Cluster1") &  group == TRUE  ~  "#706C6C",
        # str_detect(source, "Cluster1") &  group == FALSE ~  "#D3D3D3",
        # str_detect(source, "Cluster2") &  group == TRUE  ~  "#706C6C",
        # str_detect(source, "Cluster2") &  group == FALSE ~  "#D3D3D3",
        ##############
        # Exposure
        str_detect(source, "Hg") &  group == TRUE  ~  "red",
            # Outcome
        str_detect(target, "outcome") &  group == TRUE  ~  "red",
        # Methylation 
        str_detect(target, "tc") &  group == TRUE  ~  "#bf9000",
        str_detect(target, "tc") &  group == FALSE ~  "#ffd966",
        # Transcriptome
        str_detect(target, "cg") &  group == TRUE  ~  "#38761d",
        str_detect(target, "cg") &  group == FALSE ~  "#b6d7a8",
        # proteome
        str_detect(target, "miR") &  group == TRUE  ~  "#a64d79",
        str_detect(target, "miR") &  group == FALSE ~  "#ead1dc",
        ##
        group == FALSE ~ "#D3D3D3", # Negative association
        group == TRUE ~  "#706C6C")) # Positive association
  
  links_new1<- links_new %>%
   dplyr::select(colnames(links_new), target)
    
  plotly_link <- list(
    source = links_new1$IDsource,
    target = links_new1$IDtarget,
    value = links_new1$value+.00000000000000000000001, 
    color = links_new1$link_color)  
  
  # Get list of nodes for Plotly
  plotly_node <- list(
    label = nodes_new_plotly1$name, 
    color = nodes_new_plotly1$color,
    pad = 15,
    thickness = 20,
    line = list(color = "black",width = 0.5),
    x = nodes_new_plotly1$x, 
    # y = c(0.01, 
    #       0.3, 0.7, # clusters
    #       seq(from = .01, to = 1, by = 0.04)[1:(dimZ * 0.25)], # biomaker
    #       .95
    y = c(0.01,
          0.1, 0.5, # clusters
          seq(from = .05, to = 1, by = 0.04)[1:21],
          # seq(from = (.01+0.06*7), to = 1, by = 0.08)[1:5],
          # 0.9,
          # biomaker
          0.98
  ))
  
  
  ## 6.3 Plot Figure ----
  (fig <- plot_ly(
    type = "sankey",
    domain = list(
      x =  c(0,1),
      y =  c(0,1)),
    orientation = "h",
    node = plotly_node,
    link = plotly_link))
  
  (fig <- fig %>% layout(
    # title = "Basic Sankey Diagram",
    font = list(
      size = text_size
    ))
  )
  return(fig)
}
