#' Gumbel's copula
#'
#' @description returns the density function of Gumbel's copula 
#' @export 
#' @keywords internal 
#'
#'

dGumbel <- function(u,v,theta){
  u1 <- -log(u)
  u2 <- -log(v)
  num1 <- exp(-(u1^theta + u2^theta)^(1/theta))
  num2 <- (u*v)^(-1)
  num3 <- (u1*u2)^(theta-1)
  num4 <- (u1^theta + u2^theta)^(1/theta) + theta - 1
  num <- num1 * num2 * num3 * num4
  den <- (u1^theta + u2^theta)^(2-1/theta)
  num/den
}