#' Plot graphics related to the composite endpoint.
#'
#' @description Plot the survival function and the HR for composite endpoint over time and the ARE (Assymptotic Relative Efficiency) and sample size
#' size according to the correlation.The composite endpoint is assumed to be a time to event endpoint formed by a combination of two events (E1 and E2). We assume that the endpoint 1 is more relevant for the clinical question than endpoint 2. 
#' #' 
#'
#' @param p0_e1 numeric parameter between 0 and 1, expected proportion of observed events for the endpoint E1
#' @param p0_e2 numeric parameter between 0 and 1, expected proportion of observed events for the endpoint E2
#' @param HR_e1 numeric parameter between 0 and 1, expected cause specific hazard Ratio the endpoint E1
#' @param HR_e2 numeric parameter between 0 and 1, expected cause specific hazard Ratio the endpoint E2
#' @param beta_e1 numeric positive parameter, shape parameter (\eqn{\beta_1}) for a Weibull distribution for the endpoint E1 in the control group. See details for more info.
#' @param beta_e2 numeric positive parameter, shape parameter (\eqn{\beta_2}) for a Weibull distribution for the endpoint E2 in the control group. See details for more info.
#' @param case integer parameter in \{1,2,3,4\}: (1) none of the endpoints is death; (2) endpoint 2 is death; (3) endpoint 1 is death; (4) both endpoints are death by different causes. 
#' @param copula character indicating the copula to be used: "Frank" (default), "Gumbel" or "Clayton". See details for more info.
#' @param rho numeric parameter between -1 and 1, Spearman's correlation coefficient o Kendall Tau between the marginal distribution of the times to the two events E1 and E2. See details for more info.
#' @param rho_type character indicating the type of correlation to be used: "Spearman" (default) or "Tau". See details for more info.
#' @param followup_time numeric parameter indicating the maximum follow up time (in any unit). Default is 1.
#' @param alpha numeric parameter. The probability of type I error. By default \eqn{\alpha=0.05}
#' @param power numeric parameter. The power to detect the treatment effect. By default \eqn{1-\beta=0.80}
#' @param ss_formula character indicating the formula to be used for the sample size calculation on the single components: 'schoenfeld' (default) or 'freedman' 
#' 
#' @import ggpubr
#' @export 
#'
#' @return Four plots related to composite endpoint are returned:
#' \describe{
#' \item{S}{Survival curve for the composite endpoint over time}
#' \item{HR}{Hazard Ratio for the composite endpoint over time}
#' \item{ARE}{ARE according to correlation (\eqn{\rho})}
#' \item{SS}{Sample size for the composite endpoint according to correlation (\eqn{\rho})} 
#' }
#' 
#' @details Some parameters might be difficult to anticipate, especially the shape parameters of Weibull distributions and those referred to the relationship between the marginal distributions. 
#' For the shape parameters (beta_e1, beta_e2) of the Weibull distribution, we recommend to use \eqn{\beta_j=0.5}, \eqn{\beta_j=1} or \eqn{\beta_j=2} if a decreasing, constant or increasing rates over time are expected, respectively.
#' For the correlation (rho) between both endpoints, generally a positive value is expected as it has no sense to design an study with two endpoints negatively correlated. We recommend to use \eqn{\rho=0.1}, \eqn{\rho=0.3} or \eqn{\rho=0.5} for weak, mild and moderate correlations, respectively.
#' For the type of correlation (rho_type), although two different type of correlations are implemented, we recommend the use of the Spearman's correlation.
#' In any case, if no information is available on these parameters, we recommend to use the default values provided by the function.
#'
#'
#'


plot_tte <- function(p0_e1, p0_e2, HR_e1, HR_e2, beta_e1=1, beta_e2=1, case, 
                     copula = 'Frank', rho=0.3, rho_type='Spearman',
                     followup_time=1,
                     alpha=0.05, power=0.80 ,ss_formula='schoenfeld'){ 
  
  
  requireNamespace("stats")
  if(p0_e1 < 0 || p0_e1 > 1){
    stop("The probability of observing the event E1 (p_e1) must be a number between 0 and 1")
  }else if(p0_e2 < 0 || p0_e2 > 1){
    stop("The probability of observing the event E2 (p_e2) must be a number between 0 and 1")
  }else if(HR_e1 < 0 || HR_e1 > 1){
    stop("The hazard ratio for the relevant endpoint E1 (HR_e1) must be a number between 0 and 1")
  }else if(HR_e2 < 0 || HR_e2 > 1){
    stop("The hazard ratio for the secondary endpoint E2 (HR_e2) must be a number between 0 and 1")
  }else if(beta_e1 <= 0){
    stop("The shape parameter for the marginal weibull distribution of the relevant endpoint E1 (beta_e1) must be a positive number")
  }else if(beta_e2 <= 0){
    stop("The shape parameter for the marginal weibull distribution of the secondary endpoint E2 (beta_e2) must be a positive number")
  }else if(!case %in% 1:4){
    stop("The case (case) must be a number in {1,2,3,4}. See ?ARE_tte")
  }else if(!copula %in% c('Frank','Gumbel','Clayton')){
    stop("The copula (copula) must be one of 'Frank','Gumbel' or 'Clayton'")
  }else if(rho < -1 || rho > 1){
    stop("The correlation (rho) must be a number between -1 and 1")
  }else if(!rho_type %in% c('Spearman','Kendall')){
    stop("The correlation type (rho_type) must be one of 'Spearman' or 'Kendall'")
  }else if(case==4 && p0_e1 + p0_e2 > 1){
    stop("The sum of the proportions of observed events in both endpoints in case 4 must be lower than 1")
  }else if(!(is.numeric(followup_time) && followup_time>0)){
    stop("The followup_time must be a positive numeric value")
  }else if(alpha<=0 || alpha>=1){
    stop("The probability of type I error (alpha) must be a numeric value between 0 and 1")
  }else if(power<=0 || power>=1){
    stop("The power must be a numeric value between 0 and 1")
  }else if(!ss_formula %in% c('schoenfeld','freedman')){
    stop("The selected formula (ss_formula) must be one of 'schoenfeld' (default) or 'freedman'")
  }
  

  invisible(capture.output(plot_surv   <- surv_tte(p0_e1, p0_e2, HR_e1, HR_e2, beta_e1=beta_e1, beta_e2=beta_e2, case=case, 
                                copula = copula, rho=rho, rho_type=rho_type,
                                followup_time=followup_time,
                                plot_res=FALSE,plot_store=TRUE)$gg_object))
  
  invisible(capture.output(plot_effect <- effectsize_tte(p0_e1, p0_e2, HR_e1, HR_e2, beta_e1=beta_e1, beta_e2=beta_e2, case=case,  
                                copula = copula, rho=rho, rho_type=rho_type,
                                followup_time=followup_time,
                                plot_res=FALSE,plot_store=TRUE)$gg_object))
  
  invisible(capture.output(plot_ARE    <- ARE_tte(p0_e1, p0_e2, HR_e1, HR_e2, beta_e1=beta_e1, beta_e2=beta_e2, case=case,  
                                copula = copula, rho=rho, rho_type=rho_type, 
                                plot_res=FALSE,plot_store=TRUE)$gg_object))
  
  invisible(capture.output(plot_ss     <- samplesize_tte(p0_e1, p0_e2, HR_e1, HR_e2, beta_e1=1, beta_e2=1, case, 
                                copula = copula, rho=rho, rho_type=rho_type, 
                                plot_res=FALSE,plot_store=TRUE, 
                                alpha=0.05, power=0.80 ,ss_formula='schoenfeld')$gg_object))
  
    
  print(ggarrange(plot_surv,
                  plot_effect,
                  plot_ARE,
                  plot_ss,ncol=2,nrow=2))
  
  
}


