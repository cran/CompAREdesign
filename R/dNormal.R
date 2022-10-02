#' Normal copula
#'
#' @description returns the density function of Normal copula  
#' @export 
#' @keywords internal 
#'
#'

dNormal <- function(u,v,theta){
  x <- qnorm(u)
  y <- qnorm(v)
  
  (1 - theta^2)^(-1/2) * exp(-(x^2 + y^2 - 2*theta*x*y)/(2*(1-theta^2))) * exp((x^2 + y^2)/2)
}