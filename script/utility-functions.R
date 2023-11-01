#' Create a y dataset of concentrations with additional information
#' need to compute concentration curves
#' 
#' @output 3-D array, of dim (N subjects, 4, max timepoints)
#' The second dimension is of order concentration, times, mask, dose
#' @param pk_df pk data in a long format with columns of ID, TIME, DV, DOSE. Can be a superset of IDS of ids_df
#' @param ids_df Tibble of ids to pull pk data for. This order will be retained
create_masking <- function(pk_df, ids_df) {
  
  assertthat::see_if(all(c("TIME", "DV", "DOSE", "ID") %in% colnames(pk_df)), 
                     msg = "pk_df missing required column")
  
  assertthat::see_if("ID" %in% colnames(ids_df), 
                     msg = "ids_df must have ID column")
  assertthat::see_if(length(ids_df$ID) == length(unique(ids_df$ID)), 
                     msg = "IDs in ids_df must be unique")
  
  
  y_train1 = ids_df %>%
    select(ID) %>%
    inner_join(pk_df, by = "ID") %>%
    group_by(ID) %>%
    nest(TIME = TIME, DV = DV, DOSE = DOSE) %>%
    mutate(TIME = map(TIME, ~.x[[1]]),
           DV = map(DV, ~.x[[1]]),
           DOSE = map(DOSE, ~.x[[1]]))
  
  concs <-  tf$keras$preprocessing$sequence$pad_sequences(y_train1$DV, padding = "post", dtype = "float")
  
  times <- tf$keras$preprocessing$sequence$pad_sequences(y_train1$TIME, padding = "post", dtype = "float", value = max(pk_df$TIME) + 1)
  
  # 1/TRUE means to keep value, 0/FALSE to discard
  mask <- 1L*(tf$keras$preprocessing$sequence$pad_sequences(y_train1$DV, padding = "post", dtype = "bool"))
  
  dose <-  tf$keras$preprocessing$sequence$pad_sequences(y_train1$DOSE, padding = "post", dtype = "float")
  
  aperm(simplify2array(list(concs, times, mask, dose)), c(1, 3, 2))
}