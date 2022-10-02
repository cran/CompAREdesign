#' Probability of composite binary endpoints
#'
#'
#' @description This function calculates the probability of the composite binary endpoint formed by a combination of two events (E1 and E2).
#' This probability is calculated by means of the probabilities of the composite components (E1 and E2) and the correlation between them in terms of Pearson's correlation coefficient.
#'
#'
#' @param p_e1 numeric parameter, probability of the event E1
#' @param p_e2 numeric parameter, probability of the event E2
#' @param rho numeric parameter, Pearson's correlation between E1 and E2
#'
#' @export
#'
#' @return Returns the probability of the composite endpoint (E1 or E2).
#'
#' @details The input parameters stand for the probability of the composite components and Pearson's correlation between the two components.
#' Note that Pearson's correlation takes values between two bounds that depend on the probabilities p0_e1 and p0_e2.
#' To calculate the correlation bounds you can use the R functions lower_corr and upper_corr, available in this package.
#' 
#' @examples 
#' CompAREdesign::prob_cbe(p_e1=0.1, p_e2=0.2, rho=0)
#' 
#' @references Bofill Roig, M., & Gomez Melis, G. (2019). A new approach for sizing trials with composite binary endpoints using anticipated marginal values and accounting for the correlation between components. Statistics in Medicine, 38(11), 1935-1956. https://doi.org/10.1002/sim.8092
#'
#'
prob_cbe <- function(p_e1, p_e2, rho){
  if(p_e1 < 0 || p_e1 > 1){
    stop("The probability of observing the event E1 (p_e1) must be number between 0 and 1")
  }else if(p_e2 < 0 || p_e2 > 1){
    stop("The probability of observing the event E2 (p_e2) must be number between 0 and 1")
  }else if(rho < lower_corr(p_e1,p_e2)  ||  rho > upper_corr(p_e1,p_e2)){
    stop("The correlation must be in the correct interval") 
  }else{  
    prob_ce <- 1- (1-p_e1)*(1-p_e2)*( 1+ rho*sqrt(p_e1*p_e2/((1-p_e1)*(1-p_e2))))
    return(prob_ce)
  }
}
