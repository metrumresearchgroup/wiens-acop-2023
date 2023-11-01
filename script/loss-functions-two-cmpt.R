## load tf probability for use in loss function
import("tensorflow_probability")
tfp <- tf_probability()

two_cmpt_bolus_loss <- function(y_true, y_pred) {
  times_vec <- list(0.5, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 24, 48, 72, 96)
  two_cmpt_bolus_loss_time(y_true, y_pred)
}

two_cmpt_bolus_loss_time <- function(y_true, y_pred, times_vec) {
  ones_shape <- tf$ones_like(y_true)
  
  times <- ones_shape*tf$constant(times_vec)
  
  y_pred_exp <- tf$exp(y_pred) # Assume the params are logged from the NN
  
  .CL = y_pred_exp[,1:1]
  .V1 = y_pred_exp[,2:2]
  .Q = y_pred_exp[,3:3]
  .V2 = y_pred_exp[,4:4]
  .VAR = y_pred_exp[,5:5]
  .dose_times = times
  ## TODO: pull dose values from y_true instead of hard-coding
  .dose_levels <- tf$ones_like(.CL)*5
  preds <- two_cmpt_bolus_param_to_predictions(.dose_levels, .dose_times, .CL, .V1, .Q, .V2)
  
  ## calculate log-likelihood and return
  ## this part fits the VAR param
  llh <- tfp$distributions$Normal$log_prob(tfp$distributions$Normal(preds[[1]], .VAR), y_true)
  sumllh <- tf$reduce_sum(llh)
  return(sumllh)
}

tensor_reshape <- function(.param, .dose_times){
  tf$`repeat`(.param, repeats = .dose_times$shape[2], axis = 1L )
}

#' Calculate concentrations given 4 PK parameters
#' @param dose_levels Rank 1 tensor of the dose for each patient
#' @param obs_times Rank 2 tensor of the observation times for each patient
#' @param CL Rank 2 tensor of dimensions (n, 1) tensor of the CL for each patient 
#' @param V1 Rank 2 tensor of dimensions (n, 1) of the central volume for each patient
#' @param CL Rank 2 tensor of dimensions (n, 1) tensor of the Q for each patient
#' @param V2 Rank 2 tensor of dimensions (n, 1) tensor of the peripheral volume for each patient
two_cmpt_bolus_param_to_predictions <- function(dose_levels, 
                                                obs_times,
                                                CL, 
                                                V1,
                                                Q, 
                                                V2, 
                                                .scale = 1000.0) {

  CL  = tensor_reshape(CL, obs_times)
  V1  = tensor_reshape(V1, obs_times)
  Q   = tensor_reshape(Q, obs_times)
  V2  = tensor_reshape(V2, obs_times)
  
  kelim =  tf$realdiv(CL,V1)
  k12 =  tf$realdiv(Q,V1)
  k21 = k12 * tf$realdiv(V1, V2)
  beta = 0.5 * (k12 + k21 + kelim- tf$pow( tf$pow(k12 + k21 + kelim,2) - (4 * k21 * kelim) , 0.5))
  alpha = k21 * kelim / beta
  A1 = (alpha - k21) / (V1 * (alpha - beta))
  B1 = (beta - k21) / (V1 * (beta - alpha))
  
  dose_matrix <- tf$`repeat`(dose_levels, repeats = obs_times$shape[2], axis = 1L )
  
  conc_amt <- dose_matrix * ( A1 * tf$exp(-1 * alpha * obs_times) + B1 * tf$exp(-1 * beta * obs_times))
  
  conc_pred <- .scale * conc_amt
  
  return(conc_pred)
}

two_cmpt_bolus_param_to_predictions_ode <- 
  function(
    dose_levels,
    dose_times,
    CL, 
    V1,
    Q,
    V2,
    .scale = 1000.0
  ) { 
#    t_init = 0
    ode_fn = function(t, y, A) {
      return(tf$linalg$matvec(A,y)) # multiplies matrix A by vector y (initial conditions)
    }
    
    # ode system:  
    # dC1/dt = -(CL/V1 + Q/V1)*C1 + Q/V2*C2 + R_in
    # dC2/dt = Q/V1*C1 - Q/V2*C2
    
    # rate constants
    # k10 = CL/V1
    # k12 = Q/V1
    # k21 = Q/V2
    
    # "A" matrix components for system
    # A(1,1) = -(k10 +k12)
    # A(1,2) = k21
    # A(2,2) = -k21
    # A(2,1) = k12
    # where rows represent 'compartments' and columns represent 'concentrations'
    
    dose_levels <- tf$realdiv(dose_levels, V1)
    dose_levels1 <- tf$concat(list(dose_levels, tf$zeros_like(dose_levels)), axis = 1L)
    
        
    CL = tf$squeeze(CL) # using these 'squeezes to change shape from (n,1) to (n) to fit solver
    V1 = tf$squeeze(V1)
    Q = tf$squeeze(Q)
    V2 = tf$squeeze(V2)
    
    # A11 = 0
    # A12 = 0
    # A21 = 0
    # A22 = 0
# 
#     .A = tf$constant(list(list(A11,A12), list(A21,A22)), shape = list(2L, 2L))
#     
    res = tf$map_fn(
      function(.x) {
        t_init = tf$constant(0.0)
        k10 <- tf$realdiv(.x[[1]], .x[[2]])
        k12 <- tf$realdiv(.x[[3]], .x[[2]])
        k21 <- tf$realdiv(.x[[3]], .x[[4]])
        
        A11 = -(k10 +k12)
        A12 = k21
        A22 = -k21
        A21 = k12
        
        #.A = tf$constant(list(list(A11,A12), list(A21,A22)), shape = list(2L, 2L))
        
        
        #.A$assign(list(list(A11, A12), list(A21, A22)))
         A1 = tf$stack(list(A11, A12), axis = 0L)
         A2 = tf$stack(list(A21, A22), axis = 0L)
         .A = tf$stack(list(A1,A2), axis = 0L)
        
        # browser()
   #     .A = tf$constant(.A, shape(2L, 2L))
   #     .y_init = tf$Variable(.x[[5]], shape = NULL)
        .y_init = .x[[5]]
        .solution_t = .x[[6]]
        # print("before")
        # print(.y_init)
         #print(t_init)
        #print(.solution_t)
        .y_init$set_shape(2L)
        .A$set_shape(list(2L, 2L))
        #print(.A)
        res = tfp$math$ode$DormandPrince()$solve(ode_fn, t_init, .y_init,
                                                 solution_times = .solution_t,
                                                 constants=list("A" = .A))$states[,1:1]
       # res
     #  print(CL)
      }, list(CL, V1, Q, V2, dose_levels1, dose_times),
   
      fn_output_signature = tf$float32
    )
    res *.scale
    #res[,1] * .scale
   # val = res[,1] * .scale
   # val
  }
