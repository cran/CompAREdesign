#' Sample size for composite time to event endpoints
#'
#' @description This function calculates the required sample size for trials with a composite time to event endpoint as primary endpoint.
#' The primary endpoint is assumed to be a composite time to event endpoint formed by a combination of two events (E1 and E2).
#' The sample size is computed to evaluate differences between two groups based on the log rank test.
#' The sample size is calculated on the basis of anticipated information on the composite components and the correlation between them.
#'
#' @param alpha numeric parameter. The probability of type I error. By default \eqn{\alpha=0.05}
#' @param power numeric parameter. The power to detect the treatment effect. By default \eqn{1-\beta=0.80}
#' @param ss_formula character indicating the formula to be used for the sample size calculation on the single components: 'schoendfeld' (default) or 'freedman' 
#' @inheritParams ARE_tte
#' 
#' @import copula
#' 
#' @export 
#'
#' @return A list containing the following components:
#'
#' \describe{
#'   \item{\code{ss_E1}}{Total sample size (both groups) for a trial using endpoint 1 as primary endpoint}
#'   \item{\code{ss_E2}}{Total sample size (both groups) for a trial using endpoint 2 as primary endpoint}
#'   \item{\code{ss_Ec}}{Total sample size (both groups) for a trial using composite endpoint as primary endpoint}
#' } 
#'
#' @details Some parameters might be difficult to anticipate, especially the shape parameters of Weibull distributions and those referred to the relationship between the marginal distributions. 
#' For the shape parameters (beta_e1, beta_e2) of the Weibull distribution, we recommend to use \eqn{\beta_j=0.5}, \eqn{\beta_j=1} or \eqn{\beta_j=2} if a decreasing, constant or increasing rates over time are expected, respectively.
#' For the correlation (rho) between both endpoints, generally a positive value is expected as it has no sense to design an study with two endpoints negatively correlated. We recommend to use \eqn{\rho=0.1}, \eqn{\rho=0.3} or \eqn{\rho=0.5} for weak, mild and moderate correlations, respectively.
#' For the type of correlation (rho_type), although two different type of correlations are implemented, we recommend the use of the Spearman's correlation.
#' In any case, if no information is available on these parameters, we recommend to use the default values provided by the function.
#' 
#' The user can choose between the two most common formulae (Schoendfeld and Freedman) for the sample size calculation for the single components. 
#' Schoendfeld formula always be used for the composite endpoint.
#'
#'
#' @references 
#' Friedman L.M., Furberg C.D., DeMets D.L. Fundamentals of Clinical Trials. 3rd ed. New York: Springer; 1998.
#' Cortés Martínez, J., Geskus, R.B., Kim, K. et al. Using the geometric average hazard ratio in sample size calculation for time-to-event data with composite endpoints. BMC Med Res Methodol 21, 99 (2021). https://doi.org/10.1186/s12874-021-01286-x
#'
#'
samplesize_tte <- function(p0_e1, p0_e2, HR_e1, HR_e2, beta_e1=1, beta_e2=1, case, copula = 'Frank', rho=0.3, rho_type='Spearman', alpha=0.05, power=0.80 ,ss_formula='schoendfeld'){
  
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
    stop("The case (case) must be a number in {1,2,3,4}. See ?effectsize_tte")
  }else if(!copula %in% c('Frank','Gumbel','Clayton')){
    stop("The copula (copula) must be one of 'Frank','Gumbel','Clayton'")
  }else if(rho < -1 || rho > 1){
    stop("The correlation (rho) must be a number between -1 and 1")
  }else if(!rho_type %in% c('Spearman','Kendall')){
    stop("The correlation type (rho_type) must be one of 'Spearman' or 'Kendall'")
  }else if(alpha<=0 || alpha>=1){
    stop("The probability of type I error (alpha) must be a numeric value between 0 and 1")
  }else if(power<=0 || power>=1){
    stop("The power must be a numeric value between 0 and 1")
  }else if(!ss_formula %in% c('schoendfeld','freedman')){
    stop("The selected formula (ss_formula) must be one of 'schoendfeld' (default) or 'freedman'")
  }
  
  eff_size <- effectsize_tte(p0_e1, p0_e2, HR_e1, HR_e2, beta_e1, beta_e2, case, copula, rho, rho_type, subdivisions=1000,plot_HR = FALSE)
  gAHR <- eff_size$effect_size$gAHR
  
  ##-- Events
  events_1 <- ifelse(ss_formula=='schoendfeld',
                     schoendfeld_formula(alpha,power,HR_e1),
                     freedman_formula(alpha,power,HR_e1))
  events_2 <- ifelse(ss_formula=='schoendfeld',
                     schoendfeld_formula(alpha,power,HR_e2),
                     freedman_formula(alpha,power,HR_e2))
  events_c <- schoendfeld_formula(alpha,power,gAHR)
  
  ##-- Sample size --> Falta p1_e1 y p1_e2
  p1_e1 <- eff_size$measures_by_group$p_e1[2]
  p1_e2 <- eff_size$measures_by_group$p_e2[2]
  p0_star <- eff_size$measures_by_group$pstar[1]
  p1_star <- eff_size$measures_by_group$pstar[2]
  
  ss_1 <- as.numeric(2*ceiling(events_1/(p0_e1 + p1_e1)))
  ss_2 <- as.numeric(2*ceiling(events_2/(p0_e2 + p1_e2)))
  ss_c <- as.numeric(2*ceiling(events_c/(p0_star + p1_star)))
  
  ##-- Output text
  cat('The total sample size required to conduct a trial with the first component is',format(ss_1,digits = 0,big.mark = ','),'\n',
      'The total sample size required to conduct a trial with the second component is',format(ss_2,digits = 0,big.mark = ','),'\n',
      'The total sample size required to conduct a trial with the composite endpoint is',format(ss_c,digits = 0,big.mark = ','),'\n')

  return(list('ss_E1' = ss_1,
              'ss_E2' = ss_2,
              'ss_Ec' = ss_c))
}
