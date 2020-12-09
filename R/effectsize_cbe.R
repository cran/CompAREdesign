#' Effect for composite binary endpoints
#'
#' @description This function calculates different effect measures for binary composite outcomes. 
#' The composite endpoint is assumed to be a binary endpoint formed by a combination of two events (E1 and E2).
#' The effect size is calculated on the basis of anticipated information on the composite components and the correlation between them. The function allows to compute the effect size in terms of the risk difference, risk ratio and odds ratio. 
#'
#' @param p0_e1 numeric parameter, probability of occurrence E1 in the control group
#' @param p0_e2 numeric parameter, probability of occurrence E2 in the control group
#' @param effm_e1 Effect measure used for the event E1  (effm_e1 = "diff" for difference of proportions, effm_e1 = "rr" for risk ratio, effm_e1 = "or" for odds ratio)
#' @param eff_e1 numeric parameter, anticipated effect for the composite component E1
#' @param effm_e2 Effect measure used for the event E2  (effm_e2 = "diff" for difference of proportions, effm_e2 = "rr" for risk ratio, effm_e2 = "or" for odds ratio)
#' @param eff_e2 numeric parameter, anticipated effect for the composite component E2
#' @param effm_ce Effect measure used for the composite endpoint (effm_ce = "diff" for difference of proportions, effm_ce = "rr" for risk ratio, effm_ce = "or" for odds ratio)
#' @param rho numeric parameter, Pearson's correlation between the two events E1 and E2
#'
#' @export
#'
#' @return Returns the effect for the composite binary endpoint and the effects for the composite components.
#'
#' @details The input parameters stand for the probability of the composite components and Pearson's correlation between the two components.
#' Note that Pearson's correlation takes values between two bounds that depend on the probabilities p0_e1 and p0_e2.
#' To calculate the correlation bounds you can use the R functions lower_corr and upper_corr, available in this package.
#'
#'
#' @references Bofill Roig, M., & Gómez Melis, G. (2019). A new approach for sizing trials with composite binary endpoints using anticipated marginal values and accounting for the correlation between components. Statistics in Medicine, 38(11), 1935-1956. https://doi.org/10.1002/sim.8092
#'
#'
effectsize_cbe <- function(p0_e1, p0_e2, eff_e1, effm_e1, eff_e2, effm_e2, effm_ce="diff", rho){ 

  if(p0_e1 < 0 || p0_e1 > 1){
    stop("The probability of observing the event E1 (p_e1) must be number between 0 and 1")
  }else if(p0_e2 < 0 || p0_e2 > 1){
    stop("The probability of observing the event E2 (p_e2) must be number between 0 and 1")
  }else if(effm_e1 != "diff" && effm_e1 != "rr" && effm_e1 != "or"){
    stop("You have to choose between odds ratio, relative risk or difference in proportions")
  }else if((effm_e1 == "diff" && eff_e1 > 0) || (effm_e1 == "or" && (eff_e1 < 0 || eff_e1 > 1)) || (effm_e1 == "rr" && (eff_e1 < 0 || eff_e1 > 1))){
    stop("The effect of the event E1 is not right")
  }else if(effm_e2 != "diff" && effm_e2 != "rr" && effm_e2 != "or"){
    stop("You have to choose between odds ratio, relative risk or difference in proportions")
  }else if((effm_e2 == "diff" && eff_e2 > 0) || (effm_e2 == "or" && (eff_e2 < 0 || eff_e2 > 1)) || (effm_e2 == "rr" && (eff_e2 < 0 || eff_e2 > 1))){
    stop("The effect of the event E2 is not right")
  }else if(effm_ce != "diff" && effm_ce != "rr" && effm_ce != "or"){
    stop("You have to choose between odds ratio, relative risk or difference in proportions")
  }
  # else if(rho <= lower_corr(p0_e1,p0_e2)  ||  rho >= upper_corr(p0_e1,p0_e2)){
  #   stop("The correlation must be in the correct interval")
  # }


  if(effm_e1 == "or"){
    p1_e1= (eff_e1*p0_e1/(1-p0_e1))/(1+(eff_e1*p0_e1/(1-p0_e1)))
  }else if(effm_e1 == "rr"){
    p1_e1 = eff_e1 * p0_e1
  }else if(effm_e1 == "diff"){
    p1_e1 = eff_e1 + p0_e1
  }

  if(effm_e2 == "or"){
    p1_e2 = (eff_e2*p0_e2/(1-p0_e2))/(1+(eff_e2*p0_e2/(1-p0_e2)))
  }else if(effm_e2 == "rr"){
    p1_e2 = eff_e2 * p0_e2
  }else if(effm_e2 == "diff"){
    p1_e2 = eff_e2 + p0_e2
  }
  
  if(rho < max(c(lower_corr(p0_e1,p0_e2),lower_corr(p1_e1,p1_e2)))  ||  rho > min(c(upper_corr(p0_e1,p0_e2),upper_corr(p1_e1,p1_e2)))){
    stop("The correlation must be in the correct interval")
  }
  
  p0_CBE = prob_cbe(p_e1=p0_e1, p_e2=p0_e2, rho=rho)
  p1_CBE = prob_cbe(p_e1=p1_e1, p_e2=p1_e2, rho=rho)

  if(effm_ce == "diff"){
    diff_e1 = p1_e1 - p0_e1
    diff_e2 = p1_e2 - p0_e2
    effect = p1_CBE - p1_CBE
    effect_out <- data.frame(diff_e1,diff_e2,effect)
  }else if(effm_ce == "rr"){
    rr_e1 = p1_e1 / p0_e1
    rr_e2 = p1_e2 / p0_e2
    effect = p1_CBE/p0_CBE
    effect_out = data.frame(rr_e1,rr_e2,effect)
  }else if(effm_ce == "or"){
    O10= p0_e1/(1-p0_e1)
    O20= p0_e2/(1-p0_e2)
    or_e1 = (p1_e1/(1-p1_e1))/(p0_e1/(1-p0_e1))
    or_e2 = (p1_e2/(1-p1_e2))/(p0_e2/(1-p0_e2))
    effect = ((O10*or_e1+1)*(O20*or_e2+1)-1-rho*sqrt(or_e1*or_e2*O10*O20))*(1+rho*sqrt(O10*O20))/
      (((1+O10)*(1+O20)-1-rho*sqrt(O10*O20))*(1+rho*sqrt(or_e1*or_e2*O10*O20)))
    effect_out = data.frame(or_e1,or_e2,effect)
  }

  colnames(effect_out) <- c("Effect E1","Effect E2","Effect CE")
  return(effect_out)

}
