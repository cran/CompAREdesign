#' Simulation of binary composite endpoints
#'
#' @description This simulates two-arm randomised controlled trials with binary composite endpoints. The composite endpoint is assumed to be an endpoint formed by a combination of two events (E1 and E2).
#'
#' @param p0_e1 numeric parameter, probability of occurrence E1 in the control group
#' @param p0_e2 numeric parameter, probability of occurrence E2 in the control group
#' @param effm_e1 Effect measure used for the event E1  (effm_e1 = "diff" for difference of proportions, effm_e1 = "rr" for risk ratio, effm_e1 = "or" for odds ratio)
#' @param eff_e1 numeric parameter, anticipated effect for the composite component E1
#' @param effm_e2 Effect measure used for the event E2  (effm_e2 = "diff" for difference of proportions, effm_e2 = "rr" for risk ratio, effm_e2 = "or" for odds ratio)
#' @param eff_e2 numeric parameter, anticipated effect for the composite component E2
#' @param rho numeric parameter, Pearson's correlation between the two events E1 and E2
#' @param samplesize sample size per arm
#'
#' @export
#' @import stats
#' @return Simulated data
#'
#' @details The input parameters stand for the probability of the composite components and Pearson's correlation between the two components.
#' Note that Pearson's correlation takes values between two bounds that depend on the probabilities p0_e1 and p0_e2.
#' To calculate the correlation bounds you can use the R functions lower_corr and upper_corr, available in this package.
#'
#'
#'

simula_cbe <- function(p0_e1, p0_e2, eff_e1, effm_e1, eff_e2, effm_e2, rho, samplesize){

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
  }
  # effect sizes e1
  if(effm_e1 == "or"){
    p1_e1= (eff_e1*p0_e1/(1-p0_e1))/(1+(eff_e1*p0_e1/(1-p0_e1)))
  }else if(effm_e1 == "rr"){
    p1_e1 = eff_e1 * p0_e1
  }else if(effm_e1 == "diff"){
    p1_e1 = eff_e1 + p0_e1
  }
  # effect sizes e2
  if(effm_e2 == "or"){
    p1_e2 = (eff_e2*p0_e2/(1-p0_e2))/(1+(eff_e2*p0_e2/(1-p0_e2)))
  }else if(effm_e2 == "rr"){
    p1_e2 = eff_e2 * p0_e2
  }else if(effm_e2 == "diff"){
    p1_e2 = eff_e2 + p0_e2
  }
  # correlation
  if(rho < max(c(lower_corr(p0_e1,p0_e2),lower_corr(p1_e1,p1_e2)))  ||  rho > min(c(upper_corr(p0_e1,p0_e2),upper_corr(p1_e1,p1_e2)))){
    stop("The correlation must be in the correct interval")
  }

  p0_CBE = prob_cbe(p_e1=p0_e1, p_e2=p0_e2, rho=rho)
  p1_CBE = prob_cbe(p_e1=p1_e1, p_e2=p1_e2, rho=rho)

  # Group0
  data0=f_sim(samplesize=samplesize,p_e1=p0_e1,p_e2=p0_e2,p_ce=p0_CBE)
  data0$treatment=rep(0, samplesize)

  # Group1
  data1=f_sim(samplesize=samplesize,p_e1=p1_e1,p_e2=p1_e2,p_ce=p1_CBE)
  data1$treatment=rep(1, samplesize)

  # output
  data=rbind(data0,data1)

  return(data)

}


#######################





