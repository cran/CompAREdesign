#' Simulation of composite time-to-event endpoints
#'
#' @description This function simulates time-to-event components and their pertinent composite 
#' endpoint via copulas.
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
#' @param sample_size sample size for each arm (treated and control)
#' 
#' @import rootSolve
#' @rawNamespace import(copula, except = c(profile,coef,logLik,confint))
#' @rawNamespace import(numDeriv, except = hessian)
#' 
#' @export 
#'
#' @return A data.frame with 7 colums: 
#' 
#' \describe{
#'     \item{\code{time_e1}}{time to event for endpoint 1}
#'     \item{\code{status_e1}}{The status indicator of endpoint 1, 0=censored, 1=event}
#'     \item{\code{time_e2}}{time to event for endpoint 2}
#'     \item{\code{status_e2}}{The status indicator of endpoint 2, 0=censored, 1=event}
#'     \item{\code{time_ce}}{time to event for composite endpoint}
#'     \item{\code{status_ce}}{The status indicator of the composite endpoint, 0=censored, 1=event}
#'     \item{\code{treated}}{0 if control arm and 1, otherwise}
#' }
#' 
#'     
#' @details If \code{sample_size} is not an integer, it is rounded to the nearest integer.
#'
#'
#'
simula_tte <- function(p0_e1, p0_e2, HR_e1, HR_e2, beta_e1=1, beta_e2=1, 
                           case, copula = 'Frank', rho=0.3, rho_type='Spearman',
                           followup_time=1,sample_size){
 
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
    stop("The correlation (rho) must be a number between -1 and 1 and a number different from 0")
  }else if(!rho_type %in% c('Spearman','Kendall')){
    stop("The correlation type (rho_type) must be one of 'Spearman' or 'Kendall'")
  }else if(!(is.numeric(followup_time) && followup_time>0)){
    stop("The followup_time must be a positive numeric value")      
  }else if(case==4 && p0_e1 + p0_e2 > 1){
    stop("The sum of the proportions of observed events in both endpoints in case 4 must be lower than 1")
  }else if(!is.numeric(sample_size)){
    stop("The sample_size should be numeric")
  }
  
  ############################################################
  # Sample size for each group
  ############################################################
  sample_size <- round(sample_size)

  ###################################################
  ##-- Estimate parameters for distributions
  ###################################################
  ##-- Find parameters
  theta <- CopulaSelection(copula=copula,rho=rho,rho_type=rho_type)[[2]]
  MarginSelec <- MarginalsSelection(beta_e1,beta_e2,HR_e1,HR_e2,
                                    p0_e1,p0_e2,case,rho=rho,theta=theta,copula=copula)
  par_shape <- c(beta_e1,beta_e2,beta_e1,beta_e1)              # Weibull shape parameters
  par_scale <- c(MarginSelec[[5]][[2]],                        # Weibull scale parameters
                 MarginSelec[[6]][[2]],
                 MarginSelec[[7]][[2]],
                 MarginSelec[[8]][[2]])   
      
  ##-- Select copula
  if(copula=='Frank')  cop <- frankCopula(param=theta, dim = 2)
  if(copula=='Clayton')cop <- claytonCopula(param=theta, dim = 2)
  if(copula=='Gumbel') cop <- gumbelCopula(param=theta, dim = 2)
  
  ############################################################
  # Generate data
  ############################################################
  ##-- Control arm
  MVDC0 <- mvdc(copula = cop, 
                margins = c('weibull','weibull'), 
                paramMargins =  list(list(shape = par_shape[1], scale = par_scale[1]), 
                                     list(shape = par_shape[2], scale = par_scale[2])), 
                marginsIdentical = FALSE,
                check = TRUE, fixupNames = TRUE)
  BI0 <- rMvdc(sample_size, MVDC0)
  T10 <- followup_time * BI0[,1]
  T20 <- followup_time * BI0[,2]
      
  ##-- Treated arm
  MVDC1 <- mvdc(copula = cop, 
                margins = c('weibull','weibull'), 
                paramMargins =  list(list(shape = par_shape[3], scale = par_scale[3]), 
                                     list(shape = par_shape[4], scale = par_scale[4])), 
                marginsIdentical = FALSE,
                check = TRUE, fixupNames = TRUE)
  BI1 <- rMvdc(sample_size, MVDC1)  # Unit time
  T11 <- followup_time * BI1[,1]
  T21 <- followup_time * BI1[,2]
      
  ##-- Time for composite
  TC0 <- pmin(T10,T20,followup_time)
  TC1 <- pmin(T11,T21,followup_time)
  
  ##-- Status
  # Endpoint 1
  T10 <- pmin(T10,followup_time)
  T11 <- pmin(T11,followup_time)
  status_10 <- as.numeric(T10<followup_time)
  status_11 <- as.numeric(T11<followup_time)  
  
  # Endpoint 2
  T20 <- pmin(T20,followup_time)
  T21 <- pmin(T21,followup_time)
  status_20 <- as.numeric(T20<followup_time)
  status_21 <- as.numeric(T21<followup_time) 
  
  # Composite endpoint
  status_C0 <- as.numeric(TC0<followup_time)
  status_C1 <- as.numeric(TC1<followup_time)
  
  ##-- Output data.frame
  df_all <- as.data.frame(cbind(rbind(cbind(time_e1=T10,status_e1=status_10),
                                      cbind(time_e1=T11,status_e1=status_11)),
                                rbind(cbind(time_e2=T20,status_e2=status_20),
                                      cbind(time_e2=T21,status_e2=status_21)),
                                rbind(cbind(time_ce=TC0,status_ce=status_C0),
                                      cbind(time_ce=TC1,status_ce=status_C1))))
  df_all$treated <- c(rep(0,sample_size),rep(1,sample_size))

  return(invisible(df_all))
}
