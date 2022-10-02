#' Difference between survival function and a given percentile
#'
#' @description Given a percentile perc, returns the difference between this percentile and the value of survival function
#' Useful for uniroot function
#' 
#' @export 
#' @keywords internal 
#'
#'


Sstar_func_perc <- function(z,perc,...) Sstar(x=z,...) - perc