#' Lower bound for Pearson's Correlation
#'
#'
#' @description Pearson's correlation between two binary outcomes takes values between two bounds defined according to the probabilities of the binary outcomes.
#' This function calculates the lower bound of the correlation based on the probabilities of two binary outcomes.
#'
#'
#' @param p_e1 numeric parameter, probability of the event E1
#' @param p_e2 numeric parameter, probability of the event E2
#'
#' @export
#' 
#' @examples 
#' CompAREdesign::lower_corr(p_e1=0.1, p_e2=0.6)
#'
#' @return Returns the minimum value that the correlation between the two outcomes can take.
#' @details lower_corr returns a numeric value between -1 and 0.
#'
lower_corr <- function(p_e1,p_e2){
  if(p_e1 < 0 || p_e1 > 1){
    stop("The probability of observing the event E1 (p_e1) must be number between 0 and 1")
  }else if(p_e2 < 0 || p_e2 > 1){
    stop("The probability of observing the event E2 (p_e2) must be number between 0 and 1")
  }else{
    lower_corr <- max(  -sqrt(p_e1*p_e2/((1-p_e1)*(1-p_e2))), -sqrt(((1-p_e1)*(1-p_e2))/(p_e1*p_e2)  ) )
    return(lower_corr)
  }
}
