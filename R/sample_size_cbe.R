#' Sample size for composite binary endpoints
#'
#' @description This function calculates the required sample size for trials with a composite binary endpoint as primary endpoint.
#' The primary endpoint is assumed to be a composite binary endpoint formed by a combination of two events (E1 and E2).
#' The sample size is computed to evaluate differences between two groups  in terms of the risk difference, risk ratio or odds ratio.
#' The sample size is calculated on the basis of anticipated information on the composite components and the correlation between them.
#'
#' @param p0_e1 numeric parameter, probability of occurrence E1 in the control group
#' @param p0_e2 numeric parameter, probability of occurrence E2 in the control group
#' @param effm_e1 Effect measure used for the event E1  (effm_e1 = "diff" for difference of proportions, effm_e1 = "rr" for risk ratio, effm_e1 = "or" for odds ratio)
#' @param eff_e1 numeric parameter, anticipated effect for the composite component E1
#' @param effm_e2 Effect measure used for the event E2  (effm_e2 = "diff" for difference of proportions, effm_e2 = "rr" for risk ratio, effm_e2 = "or" for odds ratio)
#' @param eff_e2 numeric parameter, anticipated effect for the composite component E2
#' @param effm_ce Effect measure used for the composite endpoint (effm_ce = "diff" for difference of proportions, effm_ce = "rr" for risk ratio, effm_ce = "or" for odds ratio)
#' @param rho numeric parameter, Pearson's correlation between the two events E1 and E2
#' @param alpha Type I error
#' @param beta Type II error
#' @param unpooled Variance estimate used for the sample size calculation ("TRUE" for unpooled variance estimate, and "FALSE" for pooled variance estimate).
#'
#' @export
#' @import stats
#' @return Return the total sample size for composite binary endpoints based on the anticipated values of the composite components
#' and the association between them in terms of Pearson's correlation.
#'
#' @details The input parameters stand for the probability of the composite components and Pearson's correlation between the two components.
#' Note that Pearson's correlation takes values between two bounds that depend on the probabilities p0_e1 and p0_e2.
#' To calculate the correlation bounds you can use the R functions lower_corr and upper_corr, available in this package.
#' 
#' @references Bofill Roig, M., & Gomez Melis, G. (2019). A new approach for sizing trials with composite binary endpoints using anticipated marginal values and accounting for the correlation between components. Statistics in Medicine, 38(11), 1935-1956. https://doi.org/10.1002/sim.8092
#'
#'
samplesize_cbe <- function(p0_e1, p0_e2, eff_e1, effm_e1, eff_e2, effm_e2, effm_ce="diff", rho, alpha = 0.05, beta = 0.2, unpooled = TRUE){

  requireNamespace("stats")
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
  # else if(rho <= max(c(lower_corr(p0_e1,p0_e2),lower_corr(p1_e1,p1_e2)))  ||  rho >= max(c(upper_corr(p0_e1,p0_e2),upper_corr(p1_e1,p1_e2)))){
  #   stop("The correlation must be in the correct interval")
  # }
  else if( 0 > alpha || alpha > 1){
    stop("Alpha value must be number between 0 and 1")
  }else if( 0 > beta || beta > 1){
    stop("Beta value must be number between 0 and 1")
  }else if(unpooled != TRUE && unpooled != FALSE){
    stop("You must choose between pooled and unpooled variance")
  }


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

  # p0_CBE = 1- (1-p0_e1)*(1-p0_e2)*( 1+ rho*sqrt(p0_e1*p0_e2/((1-p0_e1)*(1-p0_e2)) ))
  # p1_CBE = 1- (1-p1_e1)*(1-p1_e2)*( 1+ rho*sqrt(p1_e1*p1_e2/((1-p1_e1)*(1-p1_e2)) ))


  if(effm_ce == "rr"){
    rr_CBE <- effectsize_cbe(p0_e1, p0_e2, eff_e1=p1_e1-p0_e1, eff_e2=p1_e2-p0_e2, effm_e1 = "diff", effm_e2 = "diff", rho, effm_ce = "rr")[1,3]
    if(unpooled == TRUE){
      samp = 2*(((qnorm(1-alpha,0,1)+qnorm(1-beta,0,1))/(log(rr_CBE)))^2*( (1-rr_CBE*p0_CBE)/(rr_CBE*p0_CBE) + (1-p0_CBE)/p0_CBE))
    }else{
      p = (rr_CBE*p0_CBE + p0_CBE)/2 
      samp = 2*(( (qnorm(1-alpha,0,1)* sqrt(2*(1-p)/p) + qnorm(1-beta,0,1)* sqrt( (1-rr_CBE*p0_CBE)/(rr_CBE*p0_CBE) + (1-p0_CBE)/p0_CBE) )/(log(rr_CBE)) )^2  )
    }
  }else if(effm_ce == "or"){
    or_CBE = effectsize_cbe(p0_e1, p0_e2, eff_e1=p1_e1-p0_e1, eff_e2=p1_e2-p0_e2, effm_e1 = "diff", effm_e2 = "diff", rho, effm_ce = "or")[1,3]
    p1 = (or_CBE*p0_CBE/(1-p0_CBE))/(1+(or_CBE*p0_CBE/(1-p0_CBE)))
    if(unpooled == TRUE){
      samp = 2*((qnorm(1-alpha,0,1)+qnorm(1-beta,0,1))/(log(or_CBE)))^2*( 1/(p0_CBE*(1-p0_CBE)) + 1/(p1*(1-p1)) )
    }else{
      samp = 2*((qnorm(1-alpha,0,1)* sqrt(2/(((p1 + p0_CBE)/2)*(1-((p1 + p0_CBE)/2)))) + qnorm(1-beta,0,1)* sqrt(1/(p0_CBE*(1-p0_CBE)) + 1/(p1*(1-p1))))/(log(or_CBE)))^2
    }
  }else{
    diff_CBE <- effectsize_cbe(p0_e1, p0_e2, eff_e1=p1_e1-p0_e1, eff_e2=p1_e2-p0_e2, effm_e1 = "diff", effm_e2 = "diff", rho, effm_ce = "diff")[1,3]
    if(unpooled== TRUE){
      samp = 2*((qnorm(1-alpha,0,1) +qnorm(1-beta,0,1))/diff_CBE)^2*( p0_CBE*(1-p0_CBE) + (diff_CBE+p0_CBE)*(1-p0_CBE-diff_CBE))
    }else{
      p = (diff_CBE + 2 * p0_CBE)/2 
      samp = 2*((qnorm(1-alpha,0,1)* sqrt(2*p*(1-p)) +  qnorm(1-beta,0,1)* sqrt( p0_CBE*(1-p0_CBE) + (diff_CBE+p0_CBE)*(1-p0_CBE-diff_CBE)))/diff_CBE)^2
    }

  }
  return(samp)

}
