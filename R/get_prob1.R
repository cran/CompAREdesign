#' Probability of observing the event
#' 
#' @description   Returns the probability of observing a single endpoint
#' 
#' @param beta1	   Shape parameter for a Weibull law for the relevant event
#' @param beta2    Shape parameter for a Weibull law for the additional event 
#' @param b11      Scale parameter of the Weibull distribution in treated arm for the relevant event
#' @param b21      Scale parameter of the Weibull distribution in treated arm for the additional event
#' @param case     Censoring case: 1, 2, 3 or 4 
#' @param rho      Spearman's coefficient between the 2 marginal distributions
#' @param copula   Copula to use
#' @param endpoint Endpoint from which it is wanted to obtain the probability of being observed
#'
#' @export 
#' @keywords internal 
#'
#'
get_prob1 <- function(beta1,beta2,b11,b21,case,rho,copula='Frank',endpoint, seed = 12345){
  
  ##-- Build copula
  copula0 <- CopulaSelection(copula=copula,rho=rho,rho_type='Spearman')
  which.copula <- copula0[[1]]
  theta <- copula0[[2]]
  if(copula=="Frank")   which.copula <-  archmCopula(family = "frank", dim = 2, param = theta)
  if(copula=="Gumbel")  which.copula <-  archmCopula(family = "gumbel", dim = 2, param = theta)
  if(copula=="Clayton") which.copula <-  archmCopula(family = "clayton", dim = 2, param = theta)
  distribution1 <- mvdc(copula = which.copula, 
                        margins = c("weibull", "weibull"), 
                        paramMargins = list(list(shape = beta1, scale = b11),
                                            list(shape = beta2, scale = b21)))
  
  ##-- Function to calculate probabilities
  prop_p <- function(d,p) sum(d[,p]<d[,3-p])/nrow(d)
  
  ##-- Calculate probabilities
  if(!is.na(seed)) set.seed(seed)
  d  <- rMvdc(10000,distribution1)
  p1 <- prop_p(d=d,p=endpoint)
  
  return(p1)
}

