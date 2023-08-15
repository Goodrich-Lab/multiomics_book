### Function --plot_mediator
# Set data matrix
plot_mediation <- function(n_mediators, mediator_names, 
                           position,
                           alpha, beta, 
                           percent_total_effect, 
                           color, 
                           main_title){
  
  dir_effect <- "'Total Effect: 0.097'"
  n_nodes <- n_mediators+2
  
  M <- matrix(nrow=n_nodes, ncol=n_nodes, data=0)
  rownames(M) <- c(mediator_names, "AS", "Liver Enzyme")
  colnames(M) <- c(mediator_names, "AS", "Liver Enzyme")
  
  M[,n_mediators+1] <- c(alpha, 0,0) 
  M[n_mediators+2,] <- c(beta, 0,0)
  M[n_mediators+2, n_mediators+1] <- dir_effect
  
  # Add color for exposure and outcome 
  color = c(color, "white", "white")
  
  # Get position
  pos <- cbind(c(rep(0.5, n_mediators), 0.1, 0.9), 
               c(position, .1, .1))
  
  pp <- plotmat(M, 
                pos = pos, 
                name = rownames(M),
                lwd = 1, 
                box.lwd = 1, 
                
                # Arrow lines
                arr.lcol = "grey50", # Color of arrow line 
                arr.tcol = "black", # Color of arrow text
                # arr.lwd = (abs(M)*)^2,# Arrow size
                arr.col = "grey50",
                arr.pos = .5,
                arr.length = 0.1,
                arr.width = 0.1,
                dtext = -0.5, # Arrow text relative to arrow position
                cex.txt = 1, # Text size for arrows
                
                shadow.col = NULL,
                box.size = c(percent_total_effect/4, 0.1, 0.1),
                box.col = color,
                # box.cex = c((abs(percent_total_effect)*4), 1, 1),
                box.cex = c(rep(1, n_mediators), 1, 1),
                # box.lcol = 
                box.type = c(rep("circle", n_mediators), "rect", "rect"),
                box.prop = 0.25, 
                curve = 0, 
                main = main_title)
  
  
}


plot_mediation_methylome <- function(n_mediators, mediator_names, 
                                     position,
                                     alpha, beta, 
                                     percent_total_effect, 
                                     color, 
                                     main_title){
  
  dir_effect <- "'Total Effect: 0.097'"
  n_nodes <- n_mediators+2
  
  M <- matrix(nrow=n_nodes, ncol=n_nodes, data=0)
  rownames(M) <- c(mediator_names, "AS", "Liver Enzyme")
  colnames(M) <- c(mediator_names, "AS", "Liver Enzyme")
  
  M[,n_mediators+1] <- c(alpha, 0,0) 
  M[n_mediators+2,] <- c(beta, 0,0)
  M[n_mediators+2, n_mediators+1] <- dir_effect
  
  # Add color for exposure and outcome 
  color = c(color, "white", "white")
  
  # Get position
  pos <- cbind(c(rep(0.5, n_mediators), 0.1, 0.9), 
               c(position, .1, .1))
  
  pp <- plotmat(M, 
                pos = pos, 
                name = rownames(M),
                lwd = 1, 
                box.lwd = 1, 
                
                # Arrow lines
                arr.lcol = "grey50", # Color of arrow line 
                arr.tcol = "black", # Color of arrow text
                # arr.lwd = (abs(M)*)^2,# Arrow size
                arr.col = "grey50",
                arr.pos = .5,
                arr.length = 0.1,
                arr.width = 0.1,
                dtext = -0.5, # Arrow text relative to arrow position
                cex.txt = 1, # Text size for arrows
                
                shadow.col = NULL,
                box.size = c(percent_total_effect/8, 0.1, 0.1),
                box.col = color,
                # box.cex = c((abs(percent_total_effect)*4), 1, 1),
                box.cex = c(rep(1, n_mediators), 1, 1),
                # box.lcol = 
                box.type = c(rep("circle", n_mediators), "rect", "rect"),
                box.prop = 0.25, 
                curve = 0, 
                main = main_title)
}