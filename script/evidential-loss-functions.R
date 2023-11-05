# This script uses loss functions that estimate the RUV, and then a log-likelihood loss based on a log-normal distribution.

evidential_loss_2cpt_bolus <- function(y_true, y_pred) {
  
  concentrations <- y_true[,1,]
  times <- y_true[,2,]
  mask <- y_true[,3,]
  dose <- y_true[,4,1:1]
  
  y_pred_exp <- tf$exp(y_pred) # Assume the params are logged from the NN
  
  CL <-  y_pred_exp[,1:1]
  V <- y_pred_exp[,2:2]
  Q <- y_pred_exp[,3:3]
  V2 <- y_pred_exp[,4:4]
  sigma <- y_pred_exp[,5:5]
  
  conc_pred <- two_cmpt_bolus_param_to_predictions(dose, times, CL, V, Q, V2)
  
  flat_mask <- tf$reshape(mask, list(-1L))
  
  flat_preds <- tf$boolean_mask(tf$reshape(conc_pred, list(-1L)), flat_mask)
  
  flat_true <- tf$boolean_mask(tf$reshape(concentrations, list(-1L)), flat_mask)
  
  mse <- keras::k_sqrt(keras::k_mean(keras::k_square(k_log(flat_true + 1E-6) - k_log(flat_preds + 1E-6))))
  
  pred_dist <- tfp$distributions$Normal(keras::k_log(flat_preds + 1E-6), sigma)
  
  ll <- keras::k_sum(pred_dist$log_prob(keras::k_log(flat_true)))
  
  -2.0 * ll
}
