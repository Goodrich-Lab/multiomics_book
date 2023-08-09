
result <- result_jive
# Show var explained ----
# Define a function with two inputs: a result object and a vector of colors
var_explained_JIVE <- function(result, col = c("grey20", "grey43", "grey65")) {
  
  # Determine the length of the list of data frames contained in the result object
  l <- length(result$data)
  
  # Calculate the proportion of variation explained by the joint and individual variation
  # for each data frame using the Frobenius norm
  # Join variation
  VarJoint = rep(0, l)
  for (i in 1:l) VarJoint[i] = norm(result$joint[[i]], type = "F")^2/norm(result$data[[i]], type = "F")^2
  # Individual variation
  VarIndiv = rep(0, l)
  for (i in 1:l) VarIndiv[i] = norm(result$individual[[i]],type = "F")^2/norm(result$data[[i]], type = "F")^2
  
  # Calculate the proportion of residual variation for each data frame
  VarResid = 1 - VarJoint - VarIndiv
  
  # Set the plot margins and layout for the bar plot
  par(mar = c(5.1, 4.1, 4.1, 0))
  layout(matrix(c(1, 2), 1, 2), heights = c(5, 5), widths = c(5, 2))
  
  # Create a bar plot of the variation explained, with the joint, individual, and residual variation shown as stacked bars
  # The names of the data frames are used as the labels for the x-axis
  barplot(
    rbind(VarJoint, VarIndiv, VarResid), col = col, main = "Variation Explained",
          names.arg = names(result$data))
}

