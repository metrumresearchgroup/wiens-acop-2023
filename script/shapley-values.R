mean_abs <- function(x) {mean(abs(x))}


#' Create a tibble of shapley values
#' 
#' @param .model fit keras model
#' @param .shapley_input_df data frame of one row per ID, 
#' with columns of ID, covariates to explain
#' e.g.:
#' bake(baked, new_data = pkdat) %>% 
#'  select(ID, all_of(covs))  %>%
#'  distinct()
#' @param .covariate_data unscaled covariates, if needed. Defaults to .shapley_input_df
#' @param .used_dummy_input Flag if the model used a column of 1s as input to a functional API model
#' @param nsim number of simulations to approximate shapley values
#' @return tibble of 5 columns, ID, covariate,
#'  shapley value, covariate value, covariate quantile
compute_shapley_analysis <- function(.model,
                                     .shapley_input_df,
                                     .covariate_data = NA,
                                     .used_dummy_input = F,
                                     .used_wt_input = F,
                                     nsim = 1) {
  
  
  if(!.used_dummy_input & !.used_wt_input) {
    shap_CL_pred_wrapper <- function(object, newdata) {
      
      model_preds <- predict(object, newdata)
      model_preds[,1]
      
    }} else if(.used_wt_input) {
      shap_CL_pred_wrapper <- function(object, newdata) {
        
        model_preds <- predict(object, list(X = newdata, WT = as.matrix(newdata[,1], ncol = 1)))
        model_preds[,1]
        
      }} else {
      shap_CL_pred_wrapper <- function(object, newdata) {
        
        model_preds <- predict(object, list(X = newdata, Constant = as.matrix(rep(1, nrow(newdata)), ncol = 1)))
        #model_preds <- object(list(X = newdata, const = as.matrix(rep(1, nrow(newdata)), ncol = 1)))
        model_preds[,1]
        
      }
    }
  
  if(is.na(.covariate_data)) {cov_dat <- .shapley_input_df} else
    {cov_dat <- .covariate_data}
  
  covariate_data_long <- 
    cov_dat %>%
    select(ID, all_of(covs)) %>%
    distinct() %>%
    pivot_longer(cols = -ID, names_to = "covariate", values_to = "cov_value") %>%
    group_by(covariate) %>%
    mutate(cov_quantile = percent_rank(cov_value))
  
  
  shapley_values <- 
    fastshap::explain(
      object = .model,
      pred_wrapper = shap_CL_pred_wrapper,
      X = .shapley_input_df %>% 
        select(-ID) %>%
        as.matrix,
      newdata = .shapley_input_df %>% 
        select(-ID) %>%
        as.matrix,
      nsim = nsim)
  
  
  shapley_df <- bind_cols(.shapley_input_df %>% select(ID), 
                          as.tibble(shapley_values)) %>%
    pivot_longer(cols = -ID, names_to = "covariate", values_to = "shapley_value") %>%
    inner_join(covariate_data_long, by = c("ID", "covariate")) %>%
    mutate(covariate = factor(covariate)) %>% 
    mutate(covariate = 
             forcats::fct_reorder(covariate, shapley_value, .fun = mean_abs))
  
  shapley_df
  
}

#' Create a beeswarm plot, analogous to a variable importance plot
#' 
#' @param shapley_df output from compute_shapley_analysis()
#' @return ggplot
beeswarm_shapley_plot <- function(shapley_df) {
  shapley_df %>%
    ggplot(aes(x = covariate, y = shapley_value,  color = cov_quantile)) +
    geom_jitter(alpha = 0.5, width = 0.1) +
    scale_colour_viridis_c(limits = c(0, 1), option = "plasma", direction = -1) +
    geom_abline(slope = 0, intercept = 0, colour = "darkgrey") +
    coord_flip() +
    labs(x = "Shapley Value", "Covariate", color = "Covariate Quantile")
  
}


#' Create a scatterplot of the estimated shapley values
#' 
#' @param shapley_df output from compute_shapley_analysis()
#' @param .shapley_input_df data frame of one row per ID, 
#' with columns of ID, covariates to explain
#' @param categorical_covariate string of a column in shapley_input_df to facet by;
#' Does not need to be a variable used to explain the model (e.g. a label)
#' @param label plot label for the categorical covariate
#' @param n number of individuals to plot (default all rows)
#' @return ggplot
shapley_scatter_plots <- function(shapley_df,
                                  shapley_input_df, 
                                  categorical_covariate,
                                  label = NA_character_,
                                  n = NULL) {
  
  n_indiv = coalesce(n, nrow(shapley_input_df))
  
  shapley_df %>%
    inner_join(slice_sample(shapley_input_df, n = n_indiv), by = "ID") %>%
    filter(covariate != categorical_covariate) %>%
    ggplot(aes(x = cov_value, 
               y = shapley_value,
               color = as.factor(.data[[categorical_covariate]] ))) +
    geom_point() +
    facet_wrap(~covariate, scales = "free_x") +
    labs(x = "Covariate Value",
         y = "Shapley Value", 
         color = coalesce(label, categorical_covariate))}



#' Create a smoothed plot of the estimated shapley values
#' 
#' @param shapley_df output from compute_shapley_analysis()
#' @param .shapley_input_df data frame of one row per ID, 
#' with columns of ID, covariates to explain
#' @param categorical_covariate string of a column in shapley_input_df to facet by;
#' Does not need to be a variable used to explain the model (e.g. a label)
#' @param label plot label for the categorical covariate
#' @return ggplot
shapley_smoother_plots <- function(shapley_df,
                                   shapley_input_df, 
                                   categorical_covariate,
                                   label = NA_character_) {
  shapley_df %>%
    inner_join(shapley_input_df, by = "ID") %>%
    filter(covariate != categorical_covariate) %>%
    ggplot(aes(x = cov_value, 
               y = shapley_value,
               color = as.factor(.data[[categorical_covariate]] ))) +
    geom_smooth() +
    facet_wrap(~covariate, scales = "free_x") +
    labs(x = "Covariate Value",
         y = "Shapley Value", 
         color = coalesce(label, categorical_covariate))
}



