#' Frank's copula
#'
#' @description returns the density function of Frank's copula  
#' @export 
#' @keywords internal 
#'
#'

dFrank <- function(u,v,theta){(theta*(1-exp(-theta))*exp(-theta*(u+v)))/ (exp(-theta) +  exp(-theta*(u+v))- exp(-theta*u)-exp(-theta*v))^2}