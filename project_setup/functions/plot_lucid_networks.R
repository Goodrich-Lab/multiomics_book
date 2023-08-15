# Functions for plotting netowrks for LUCIDusM


# Function that calculates the centroid of a cluster:
calculate_centroid <- function(nodes) {
  x <- mean(nodes[,1])
  y <- mean(nodes[,2])
  return(c(x, y))
}

# function that shifts the nodes in a cluster
shift_nodes <- function(nodes, shift, scale_shift = c(1,1)) {
  nodes[,1] <- nodes[,1] + scale_shift[1]*shift[1]
  nodes[,2] <- nodes[,2] + scale_shift[2]*shift[2]
  return(nodes)
}



ggraph_fxn <- function(g, x, y, label_groups = FALSE, 
                       concavity = 2, 
                       label.fontsize = 8) {
  if(label_groups){
    # ggraph(g,layout="manual",x=bb$xy[,1],y=bb$xy[,2]) +
    ggraph(g,layout="manual",x=x,y=y) +
      geom_mark_hull(aes(x, y, group = cluster_assignment,
                         color = cluster_assignment,
                         fill = cluster_assignment,
                         label=group_name
      ),
      radius = unit(1, "mm"),
      label.margin = margin(1,1,1,1,"mm"),
      label.fontsize = label.fontsize, 
      label.fontface = c("plain"),
      con.cap = unit(1, "mm"), 
      concavity = concavity,
      expand = unit(1.1, "mm"),
      alpha = 0.4) +
      geom_edge_link0(aes(edge_color = col, edge_alpha = col), width = 0.2) +
      geom_node_point(aes(col = cluster_assignment), alpha = .6) +
      scale_color_brewer(palette = "Set1") +
      scale_fill_brewer(palette = "Set1") +
      scale_edge_alpha_manual(values = c(.1, .4)) + #low risk, high risk
      scale_edge_color_manual(values = c("grey50", "black")) +
      theme_graph()+
      theme(legend.position = "none")
  } else {
    ggraph(g,layout="manual",x=x,y=y) +
      geom_mark_hull(aes(x, y, group = cluster_assignment,
                         color = cluster_assignment,
                         fill = cluster_assignment),
                     radius = unit(1, "mm"),
                     con.cap = unit(1, "mm"), 
                     concavity = concavity,
                     expand = unit(1.1, "mm"),
                     alpha = 0.4) +
      geom_edge_link0(aes(edge_color = col, edge_alpha = col), width = 0.2) +
      geom_node_point(aes(col = cluster_assignment), alpha = .6) +
      scale_color_brewer(palette = "Set1") +
      scale_fill_brewer(palette = "Set1") +
      scale_edge_alpha_manual(values = c(.1, .4)) + #low risk, high risk
      scale_edge_color_manual(values = c("grey50", "black")) +
      theme_graph()+
      theme(legend.position = "none")
  }
}
