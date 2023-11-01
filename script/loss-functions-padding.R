source(here::here("script", "loss-functions.R"))
source(here::here("script", "loss-functions-two-cmpt.R"))


logitudinal_loss_1cpt_bolus_padded <- function(y_true, y_pred) {
  
  concentrations <- y_true[,1,]
  
  #ones_shape <- tf$ones_like(concentrations)
  
  times <- y_true[,2,]
  
  mask <- y_true[,3,]
  
  dose <- y_true[,4,1:1]
  
  #times <- ones_shape*tf$constant(list(times_vec))
  # Minimal example variables for development
  # y_true <- tf$constant(list(list(1, 2, 3), list(2, 3, 4)))
  # y_pred <- tf$constant(list(list(1, 1, 1, 1), list(1.1, 1.2, 1.3, 1.4)))
  #times = tf$constant(list(list(0.5, 1, 2), list(0.5, 1, 2)))
  
  y_pred_exp <- tf$exp(y_pred) # Assume the params are logged from the NN
  
  CL <-  y_pred_exp[,1:1]
  V <- y_pred_exp[,2:2]
  
  #browser()
  
  
  #But note that the Tensor.ndim and Tensor.shape attributes don't return Tensor objects. If you need a Tensor use the tf.rank or tf.shape function. This difference is subtle, but it can be important when building graphs (later).
  
  conc_pred <- one_cmpt_bolus_param_to_predictions(dose, times, CL, V)
  
  # Would be good to log this, but the predictions can be numerically zero, which is a problem
  # This relies on the reshape operation always behaving the same here
  # See also tf$boolean_mask
  flat_mask <- tf$reshape(mask, list(-1L))
  
  flat_preds <- tf$multiply(tf$reshape(log(conc_pred + 1E-12), list(-1L)), flat_mask)
  
  flat_true <- tf$multiply(tf$reshape(log(concentrations + 1E-12), list(-1L)), flat_mask) 
  
  mse <- keras::k_sum(keras::k_square(flat_true - flat_preds)) / keras::k_sum(mask)
  
  mse
}


#' Calculate loss function 2 compartment
#' 
#' @param y_true 3D array of [subjects, dosing info, times]
#' the second dimension has rows of (observed concentration, observation time, flag for observed 1=osbserved 0=masked, dose)
#' The third dimension may be padded at the end with zeros
#' @param y_pred 2D array of [subjects, PK parameters]
#' There are 4 PK parameters (CL, V1, Q, V2)
#' @return Mean squared error loss of the logged concentrations
logitudinal_loss_2cpt_bolus_padded <- function(y_true, 
                                               y_pred) {
  
  concentrations <- y_true[,1,]
  times <- y_true[,2,]
  mask <- y_true[,3,]
  dose <- y_true[,4,1:1]
  
  y_pred_exp <- tf$exp(y_pred) # Assume the params are logged from the NN
  
  CL <-  y_pred_exp[,1:1]
  V <- y_pred_exp[,2:2]
  Q <- y_pred_exp[,3:3]
  V2 <- y_pred_exp[,4:4]
  
  conc_pred <- two_cmpt_bolus_param_to_predictions(dose, times, CL, V, Q, V2)
  
  flat_mask <- tf$reshape(mask, list(-1L))
  
  flat_preds <- tf$multiply(tf$reshape(log(conc_pred + 1E-12), list(-1L)), flat_mask)
  
  flat_true <- tf$multiply(tf$reshape(log(concentrations + 1E-12), list(-1L)), flat_mask) 
  
  mse <- keras::k_sum(keras::k_square(flat_true - flat_preds)) / keras::k_sum(mask)
  
  mse
}


logitudinal_loss_2cpt_bolus_padded_ode <- function(y_true, y_pred) {
  
  concentrations <- y_true[,1,]
  times <- y_true[,2,]
  mask <- y_true[,3,]
  dose <- y_true[,4,1:1]
  
  y_pred_exp <- tf$exp(y_pred) # Assume the params are logged from the NN
  
  CL <-  y_pred_exp[,1:1]
  V <- y_pred_exp[,2:2]
  Q <- y_pred_exp[,3:3]
  V2 <- y_pred_exp[,4:4]
  
  conc_pred <- two_cmpt_bolus_param_to_predictions_ode(dose, times, CL, V, Q, V2)
  flat_mask <- tf$reshape(mask, list(-1L))
  
  flat_preds <- tf$multiply(tf$reshape(log(conc_pred + 1E-12), list(-1L)), flat_mask)
  
  flat_true <- tf$multiply(tf$reshape(log(concentrations + 1E-12), list(-1L)), flat_mask) 
  
  mse <- keras::k_sum(keras::k_square(flat_true - flat_preds)) / keras::k_sum(mask)
  
  mse
}