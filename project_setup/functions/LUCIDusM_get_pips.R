# This function extracts the PIPs from cluster 2 for each omics layer. 
# Curretnly only works for lucid in parallel with 3 omics layers

get_pips <- function(lucid_fit, exposure){
  z <- lucid_fit$z
  N = lucid_fit$N
  
  ## PIPs for omic 1
  z_margin_1 <- t(sapply(1:N, function(j) {marginSums(z[, , , j], margin = 1)})) 
  
  z_margin_1 <- data.frame(z_margin_1) %>% 
    rename(pip_l1_c1 = X1, pip_l1_c2 = X2) %>% 
    mutate(l1_pred = ifelse(pip_l1_c1 > 0.5, 0, 1)) 
  
  ## PIPs for omic layer 2
  z_margin_2 <- t(sapply(1:N, function(j) {marginSums(z[, ,, j], margin = 2)}))
  
  z_margin_2 <- data.frame(z_margin_2) %>% 
    rename(pip_l2_c1 = X1, pip_l2_c2 = X2) %>% 
    mutate(l2_pred = ifelse(pip_l2_c1 > 0.5, 0, 1))
  
  
  ## PIPs for for omic layer 3
  z_margin_3 <- t(sapply(1:N, function(j) {marginSums(z[, , , j], margin = 3)}))
  
  z_margin_3 <- data.frame(z_margin_3) %>% 
    rename(pip_l3_c1 = X1, pip_l3_c2 = X2) %>% 
    mutate(l3_pred = ifelse(pip_l3_c1 > 0.5, 0, 1))
  
  lucid_pips <- bind_cols(z_margin_1, 
                          z_margin_2,
                          z_margin_3)
  
  # Remove cluster 1 pips (they are redundant)
  lucid_pips <- lucid_pips %>% 
    dplyr::select(-contains("c1"))
  
}
