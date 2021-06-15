#' Schoendfeld's formula
#'
#' @description Schoendfeld's formula to calculate the number of events required
#' @export 
#' @keywords internal 
#'
#'
#'
schoendfeld_formula <- function(alpha,power,HR) E <- 4*(qnorm(1-alpha/2) +  qnorm(power))^2 / (log(HR))^2
