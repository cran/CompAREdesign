#' Clayton's copula
#'
#' @description returns the density function of Clayton's copula 
#' @export 
#' @keywords internal 
#'
#'

dClayton <- function(u,v,theta){ (u*v)^(-theta-1) * (theta+1) * (u^(-theta) + v^(-theta) - 1)^(-2 - 1/theta)}
