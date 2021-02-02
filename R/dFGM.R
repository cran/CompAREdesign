#' FGM copula
#'
#' @description returns the density function of the copulaFGM
#' @export 
#' @keywords internal 
#'
#'

dFGM <- function(u,v,theta) 1 + theta * (1- 2*u) * (1 - 2*v)