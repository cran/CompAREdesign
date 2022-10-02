#' Survival function for composite endpoint
#' 
#' @description  Returns the value of the survival function of S* at point x given the marginal distributions 
#' and the bivariate distributions via copula
#'
#' @param x	        Point in which to be evaluated
#' @param dist1     Distribution function of the marginal T1 (pweibull) 
#' @param dist2     Distribution function of the marginal T2 (pweibull) 
#' @param param1    Parameters of the marginal distribution function T1 (pweibull) 
#' @param param2    Parameters of the marginal distribution function T2 (pweibull) 
#' @param dist_biv  Distribution function of the bivariate distribution via copula
#' @export 
#' @keywords internal 
#'
#'
#'

Sstar <- function(x,dist1,dist2,param1,param2,dist_biv) { 
  y <- if(length(x) == 1) c(x,x) else cbind(x,x)
  
  
  return(
    1
    - do.call(dist1,c(list(q=x),param1))
    - do.call(dist2,c(list(q=x),param2))
    + (pMvdc(y, dist_biv))
  )
}