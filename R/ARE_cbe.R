#' ARE method for composite binary endpoints
#'
#' @description The composite endpoint is assumed to be a binary endpoint formed by a combination of two events (E1 and E2). We assume that the endpoint 1 is more relevant for the clinical question than endpoint 2. 
#' This function calculates the ARE method for binary endpoints. The method quantifies the differences in efficiency of using the composite or the relevant as primary endpoint to lead the trial and, moreover, provides a decision rule to choose the primary endpoint. If the ARE is larger than 1, the composite endpoint may be considered the best option as primary endpoint. Otherwise, the relevant endpoint is preferred. 
#' 
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
#' @return Returns the ARE value. If the ARE value is larger than 1 then the composite endpoint is preferred over the relevant endpoint. Otherwise, the endpoint 1 is preferred as the primary endpoint of the study.
#'
#' @details The input parameters stand for the probability of the composite components and Pearson's correlation between the two components.
#' Note that Pearson's correlation takes values between two bounds that depend on the probabilities p0_e1 and p0_e2.
#' To calculate the correlation bounds you can use the R functions lower_corr and upper_corr, available in this package.
#'
#'
#' @references Bofill Roig, M., & Gomez Melis, G. (2018). Selection of composite binary endpoints in clinical trials. Biometrical Journal, 60(2), 246-261. https://doi.org/10.1002/bimj.201600229
#'
#'
ARE_cbe <- function(p0_e1, p0_e2, eff_e1, effm_e1= "or", eff_e2, effm_e2= "or", effm_ce="or", rho){
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
  }else if(effm_ce != "diff" && effm_ce != "or"){
    stop("You have to choose between odds ratio, or difference in proportions")
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

    if(effm_ce == "diff"){
      diff_e1 = p1_e1 - p0_e1
      diff_e2 = p1_e2 - p0_e2
      effect = prob_cbe(p1_e1,p1_e2,rho) - prob_cbe(p0_e1,p0_e2,rho)
      effect_out <- data.frame(diff_e1,diff_e2,effect)

      are <- ( (effect)^2*(p0_e1*(1-p0_e1)) )/( (diff_e1)^2*(p0_CBE*(1-p0_CBE))  )


    }else if(effm_ce == "or"){
      O10= p0_e1/(1-p0_e1)
      O20= p0_e2/(1-p0_e2)
      or_e1 = (p1_e1/(1-p1_e1))/(p0_e1/(1-p0_e1))
      or_e2 = (p1_e2/(1-p1_e2))/(p0_e2/(1-p0_e2))
      effect = ((O10*or_e1+1)*(O20*or_e2+1)-1-rho*sqrt(or_e1*or_e2*O10*O20))*(1+rho*sqrt(O10*O20))/
        (((1+O10)*(1+O20)-1-rho*sqrt(O10*O20))*(1+rho*sqrt(or_e1*or_e2*O10*O20)))
      effect_out = data.frame(or_e1,or_e2,effect)

      are <- log(effect)^2*p0_CBE*(1- p0_CBE)/(log(or_e1)^2*p0_e1*(1-p0_e1))
    }


    return(are)

}
