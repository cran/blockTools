extract_conditions <- function(assg.obj, data, id.var){
  
  condition <- rep(NA, nrow(data))
  class(condition) <- "integer"
  
  flat_assgs <- bind_rows(assg.obj$assg)
  
  for(col_idx in 1:(ncol(assg.obj$assg[[1]]) - 1)){
    
    wh_this_condition <- data[[id.var]] %in% flat_assgs[, col_idx]
    
    condition[wh_this_condition] <- col_idx
    
  }
  
  return(condition)
}
