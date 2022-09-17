#' Sample size for composite time to event endpoints
#'
#' @description This function calculates the required sample size for trials with a composite time to event endpoint as primary endpoint.
#' The primary endpoint is assumed to be a composite time to event endpoint formed by a combination of two events (E1 and E2).
#' The sample size is computed to evaluate differences between two groups based on the log rank test.
#' The sample size is calculated on the basis of anticipated information on the composite components and the correlation between them.
#'
#' @param alpha numeric parameter. The probability of type I error. By default \eqn{\alpha=0.05}
#' @param power numeric parameter. The power to detect the treatment effect. By default \eqn{1-\beta=0.80}
#' @param ss_formula character indicating the formula to be used for the sample size calculation on the single components: 'schoenfeld' (default) or 'freedman' 
#' @param subdivisions integer parameter greater than or equal to 10. Number of points used to plot the sample size according to correlation. The default is 50. Ignored if plot_res=FALSE and plot_store=FALSE.
#' @param plot_res logical indicating if the sample size according to the correlation should be displayed. The default is FALSE
#' @param plot_store logical indicating if the plot of sample size according to the correlation is stored for future customization. The default is FALSE
#' @inheritParams ARE_tte
#' 
#' @rawNamespace import(copula, except = c(profile,coef,logLik,confint))
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
#' In addition, if \code{plot_store=TRUE} an object of class \code{ggplot} with
#' the sample size for composite endpoint according to correlation is stored 
#' in the list.
#' 
#' @details Some parameters might be difficult to anticipate, especially the shape parameters of Weibull distributions and those referred to the relationship between the marginal distributions. 
#' For the shape parameters (beta_e1, beta_e2) of the Weibull distribution, we recommend to use \eqn{\beta_j=0.5}, \eqn{\beta_j=1} or \eqn{\beta_j=2} if a decreasing, constant or increasing rates over time are expected, respectively.
#' For the correlation (rho) between both endpoints, generally a positive value is expected as it has no sense to design an study with two endpoints negatively correlated. We recommend to use \eqn{\rho=0.1}, \eqn{\rho=0.3} or \eqn{\rho=0.5} for weak, mild and moderate correlations, respectively.
#' For the type of correlation (rho_type), although two different type of correlations are implemented, we recommend the use of the Spearman's correlation.
#' In any case, if no information is available on these parameters, we recommend to use the default values provided by the function.
#' 
#' The user can choose between the two most common formulae (Schoenfeld and Freedman) for the sample size calculation for the single components. 
#' Schoenfeld formula always be used for the composite endpoint.
#'
#'
#' @references 
#' Friedman L.M., Furberg C.D., DeMets D.L. Fundamentals of Clinical Trials. 3rd ed. New York: Springer; 1998.
#' Cortés Martínez, J., Geskus, R.B., Kim, K. et al. Using the geometric average hazard ratio in sample size calculation for time-to-event data with composite endpoints. BMC Med Res Methodol 21, 99 (2021). https://doi.org/10.1186/s12874-021-01286-x
#'
#'
samplesize_tte <- function(p0_e1, p0_e2, HR_e1, HR_e2, beta_e1=1, beta_e2=1, 
                           case, copula = 'Frank', rho=0.3, rho_type='Spearman', 
                           alpha=0.05, power=0.80 ,ss_formula='schoenfeld', 
                           subdivisions=50, plot_res=FALSE, plot_store=FALSE){
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
  }else if(!ss_formula %in% c('schoenfeld','freedman')){
    stop("The selected formula (ss_formula) must be one of 'schoenfeld' (default) or 'freedman'")
  }else if(!is.logical(plot_res)){
    stop("The parameter plot_res must be logical")
  }else if(!is.logical(plot_store)){
    stop("The parameter plot_store must be logical")
  }
  
  # Values of rho where to calculate Sample size
  rho_sel <- rho
  if(plot_res | plot_store){
    rho_seq <- unique(c(rho,seq(0.01,0.98,length=subdivisions)))
  }else{
    rho_seq <- rho
  }
  
  # Storage
  SS_array_1 <- SS_array_2 <- SS_array_c <- c()
  
  # Calculate Sample size for each rho
  pb = txtProgressBar(min = 0, max = length(rho_seq), initial = 0)
  for(rho in rho_seq){
    setTxtProgressBar(pb,which(rho_seq==rho))
  
    ##-- Effect size
    invisible(capture.output(eff_size <- effectsize_tte(p0_e1, p0_e2, HR_e1, HR_e2, beta_e1, beta_e2, case, copula, rho, rho_type, subdivisions=1000,plot_res = FALSE)))
    gAHR <- eff_size$effect_size$gAHR
    
    ##-- Events
    events_1 <- ifelse(ss_formula=='schoenfeld',
                       schoenfeld_formula(alpha,power,HR_e1),
                       freedman_formula(alpha,power,HR_e1))
    events_2 <- ifelse(ss_formula=='schoenfeld',
                       schoenfeld_formula(alpha,power,HR_e2),
                       freedman_formula(alpha,power,HR_e2))
    events_c <- schoenfeld_formula(alpha,power,gAHR)
    
    ##-- Probabilities of observing the event
    p1_e1 <- eff_size$measures_by_group$p_e1[2]
    p1_e2 <- eff_size$measures_by_group$p_e2[2]
    p0_star <- eff_size$measures_by_group$pstar[1]
    p1_star <- eff_size$measures_by_group$pstar[2]
    
    ##-- Sample size
    ss_1 <- as.numeric(2*ceiling(events_1/(p0_e1 + p1_e1)))
    ss_2 <- as.numeric(2*ceiling(events_2/(p0_e2 + p1_e2)))
    ss_c <- as.numeric(2*ceiling(events_c/(p0_star + p1_star)))
    
    SS_array_1 <- c(SS_array_1,ss_1)
    SS_array_2 <- c(SS_array_2,ss_2)
    SS_array_c <- c(SS_array_c,ss_c)
    
  }  
  
  if(plot_res | plot_store){
    sample_size <- NULL                # To avoid the note: "no visible binding for global variable 'sample_size'"
    dd <- data.frame(rho=rho_seq, sample_size=SS_array_c)
    gg1 <- ggplot(dd,aes(x=rho,y=sample_size)) + 
      geom_line(color='darkblue',size=1.3) +
      xlab(expression(rho)) + ylab('Sample size CE')
  }
  
  ##-- Output data.frame
  df <- data.frame(Endpoint=c('--------','Endpoint 1','Endpoint 2','Composite endpoint'),
                   "Total sample size"=c("-----------------",SS_array_1[1],SS_array_2[1],SS_array_c[1]),
                   check.names = FALSE)
  print(df, row.names = FALSE,right=FALSE)
  
  return_object <- list('ss_E1'     = SS_array_1[1],
                        'ss_E2'     = SS_array_2[1],
                        'ss_Ec'     = SS_array_c[1],
                        'gg_object' = NA)
  
  ## Print graphic
  if(plot_res) print(gg1)
  
  ## Store plot in the output
  if(plot_store) return_object$gg_object <- gg1
  
  
  ##-- Returned list
  return(invisible(return_object))
}
