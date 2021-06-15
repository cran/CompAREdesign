#' Freedman's formula
#'
#' @description Freedman's formula to calculate the number of events required
#' @export 
#' @keywords internal 
#'
#'
#'
freedman_formula  <- function(alpha,power,HR) E <- (HR+1)^2 * (qnorm(1-alpha/2) + qnorm(power))^2/(HR-1)^2