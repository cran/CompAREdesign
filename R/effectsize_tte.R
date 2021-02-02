#' Effect for composite time-to-event endpoints
#'
#' @description This function calculates different effect measures for time-to-event composite outcomes. 
#' The composite endpoint is assumed to be a time-to-event endpoint formed by a combination of two events (E1 and E2).
#' The effect size is calculated on the basis of anticipated information on the composite components and the correlation between them.
#' Marginal distributions are assumed for m
#' The function allows to compute the effect size in terms of the geometric average hazard ratio, the average hazard ratio, 
#' the ratio of restricted mean survival times and the median survival time ratio.
#' 
#' @param p0_e1 numeric parameter between 0 and 1, expected proportion of observed events for the endpoint E1
#' @param p0_e2 numeric parameter between 0 and 1, expected proportion of observed events for the endpoint E2
#' @param HR_e1 numeric parameter between 0 and 1, expected cause specific hazard Ratio the endpoint E1
#' @param HR_e2 numeric parameter between 0 and 1, expected cause specific hazard Ratio the endpoint E2
#' @param beta_e1 numeric positive parameter, shape parameter (\eqn{\beta_1}) for a Weibull distribution for the endpoint E1 in the control group. See details for more info.
#' @param beta_e2 numeric positive parameter, shape parameter (\eqn{\beta_2}) for a Weibull distribution for the endpoint E2 in the control group. See details for more info.
#' @param case integer parameter in {1,2,3,4}
#'             1: none of the endpoints is death
#'             2: endpoint 2 is death
#'             3: endpoint 1 is death
#'             4: both endpoints are death by different causes  
#' @param copula character indicating the copula to be used: "Frank" (default), "Gumbel" or "Clayton". See details for more info.
#' @param rho numeric parameter between -1 and 1, Spearman's correlation coefficient o Kendall Tau between the marginal distribution of the times to the two events E1 and E2. See details for more info.
#' @param rho_type character indicating the type of correlation to be used: "Spearman" (default) or "Tau". See details for more info.
#' @param subdivisions integer parameter greater than or equal to 10. Number of subintervals to estimate the effect size. The default is 1000. 
#' @param plot_HR logical. If the HR over time should be displayed. The default is FALSE
#' 
#' @import ggplot2
#' @import rootSolve
#' @rawNamespace import(copula, except = c(profile,coef,logLik,confint))
#' @rawNamespace import(numDeriv, except = hessian)
#' 
#' @export 
#'
#' @return A list formed by two lists: \code{effect_size}, which contains the expected treatment effect measures and \code{measures_by_group}, which contains some measures for each group
#' 
#' \code{effect_size} list:
#' 
#' \describe{
#'     \item{\code{gAHR}}{geometric Average Hazard Ratio}
#'     \item{\code{AHR}}{Average Hazard Ratio}
#'     \item{\code{RMST_ratio}}{Restricted Mean Survival Time Ratio}
#'     \item{\code{Median_ratio}}{Median Survival Time Ratio}
#' }
#' 
#' \code{measures_by_group} list:
#' 
#' \describe{
#'     \item{\code{pstar}}{array with the probability of observing the composite event for each group}
#'     \item{\code{RMST}}{array with the restricted mean survival time for each group}
#'     \item{\code{Median}}{array with the median surival time for each group}
#' }
#'   
#' @details Some parameters might be difficult to anticipate, especially the shape parameters of Weibull distributions and those referred to the relationship between the marginal distributions. 
#' For the shape parameters (beta_e1, beta_e2) of the Weibull distribution, we recommend to use \eqn{\beta_j=0.5}, \eqn{\beta_j=1} or \eqn{\beta_j=2} if a decreasing, constant or increasing rates over time are expected, respectively.
#' For the correlation (rho) between both endpoints, generally a positive value is expected as it has no sense to design an study with two endpoints negatively correlated. We recommend to use \eqn{\rho=0.1}, \eqn{\rho=0.3} or \eqn{\rho=0.5} for weak, mild and moderate correlations, respectively.
#' For the type of correlation (rho_type), although two different type of correlations are implemented, we recommend the use of the Spearman's correlation.
#' In any case, if no information is available on these parameters, we recommend to use the default values provided by the function.
#' 
#' All returned expected effect sizes for the composite endpoint should be interpreted in relative terms (treated to control). 
#' gAHR and AHR represent the risk reduction that will be achieved with the new therapy, while RMST_ratio and Median_ratio represent the gain in time gain terms until the event.
#'
#'
#' @references Schemper, M., Wakounig, S., Heinze, G. (2009). The estimation of average hazard ratios by weighted Cox regression. Stat. in Med. 28(19): 2473--2489. doi:10.1002/sim.3623
#'
#'
effectsize_tte <- function(p0_e1, p0_e2, HR_e1, HR_e2, beta_e1=1, beta_e2=1, case, copula = 'Frank', rho=0.3, rho_type='Spearman', subdivisions=1000, plot_HR=FALSE){
 
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
  }else if(!(is.numeric(subdivisions) && subdivisions>=10)){
    stop("The number of subdivisions must be an integer greater than or equal to 10")
  }else if(!is.logical(plot_HR)){
    stop("The parameter plot_HR must be logical")
  }else if(case==4 && p0_e1 + p0_e2 > 1){
    stop("The sum of the proportions of observed events in both endpoints in case 4 must be lower than 1")
  }
  
  ##-- Time (t) values to assess the treatment effect
  t <- c(0.0001,seq(0.001,1,length.out = subdivisions))
  
  ##-- Copula
  copula0 <- CopulaSelection(copula,rho,rho_type)
  which.copula <- copula0[[1]]
  theta <- copula0[[2]]   
  
  ##-- Marginal distribution and parameters
  MarginSelec <- MarginalsSelection(beta_e1,beta_e2,HR_e1,HR_e2,p0_e1,p0_e2,case,theta,copula=copula)
  
  ##-- Weibull marginal distributions
  T1dist   <- MarginSelec[[1]]
  T2dist   <- MarginSelec[[2]]
  T1pdist  <- MarginSelec[[3]]
  T2pdist  <- MarginSelec[[4]]
  
  ##-- Parameters of marginal weibull distribution
  T10param <- MarginSelec[[5]]
  T20param <- MarginSelec[[6]]
  T11param <- MarginSelec[[7]]
  T21param <- MarginSelec[[8]]
  
  ##-- Scale parameters
  b10 <- T10param[[2]]
  b20 <- T20param[[2]]
  b11 <- T11param[[2]]
  b21 <- T21param[[2]]
  
  ##-- Bivariate distribution in control and treatment groups
  distribution0 <- mvdc(copula = which.copula, margins = c(T1dist, T2dist), paramMargins = list(T10param, T20param))
  distribution1 <- mvdc(copula = which.copula, margins = c(T1dist, T2dist), paramMargins = list(T11param, T21param))
  
  
  ##-- Densities for both endpoints
  fT10 <- (beta_e1/b10) * ((t/b10)^(beta_e1-1)) * (exp(-(t/b10)^beta_e1))
  fT11 <- (beta_e1/b11) * ((t/b11)^(beta_e1-1)) * (exp(-(t/b11)^beta_e1))
  fT20 <- (beta_e2/b20) * ((t/b20)^(beta_e2-1)) * (exp(-(t/b20)^beta_e2))
  fT21 <- (beta_e2/b21) * ((t/b21)^(beta_e2-1)) * (exp(-(t/b21)^beta_e2))
  
  ##-- Survival for both endpoints
  ST10 <- exp(-(t/b10)^beta_e1)
  ST11 <- exp(-(t/b11)^beta_e1)
  ST20 <- exp(-(t/b20)^beta_e2)
  ST21 <- exp(-(t/b21)^beta_e2)
  
  ##-- Survival for the composite endpoint
  if(copula=='Frank'){
    Sstar0 <- (-log(1+(exp(-theta*ST10)-1)*(exp(-theta*ST20)-1)/(exp(-theta)-1))/theta)
    Sstar1 <- (-log(1+(exp(-theta*ST11)-1)*(exp(-theta*ST21)-1)/(exp(-theta)-1))/theta)  
  }else if(copula=='Clayton'){
    Sstar0 <- (ST10^(-theta) + ST20^(-theta) - 1)^{-1/theta}
    Sstar1 <- (ST11^(-theta) + ST21^(-theta) - 1)^{-1/theta}
  }else if(copula=='Gumbel'){
    Sstar0 <- exp(-((-log(ST10))^theta + (-log(ST20))^theta)^(1/theta))
    Sstar1 <- exp(-((-log(ST11))^theta + (-log(ST21))^theta)^(1/theta))      
  }
  
  ##-- Density, hazards and hazard ratio for the composite
  if(copula=='Frank'){
    fstar0 <- (exp(-theta*ST10)*(exp(-theta*ST20)-1)*fT10 + exp(-theta*ST20)*(exp(-theta*ST10)-1)*fT20)/(exp(-theta*Sstar0)*(exp(-theta)-1))
    fstar1 <- (exp(-theta*ST11)*(exp(-theta*ST21)-1)*fT11 + exp(-theta*ST21)*(exp(-theta*ST11)-1)*fT21)/(exp(-theta*Sstar1)*(exp(-theta)-1))
  }else if(copula=='Clayton'){
    fstar0 <- (ST10^(theta+1) * fT10 + ST20^(theta+1) * fT20)/(Sstar0*(ST10^(-theta) + ST20^(-theta) - 1))
    fstar1 <- (ST11^(theta+1) * fT11 + ST20^(theta+1) * fT21)/(Sstar1*(ST11^(-theta) + ST21^(-theta) - 1))
  }else if(copula=='Gumbel'){
    fstar0 <- Sstar0 * log(Sstar0) * ((-log(ST10))^(theta-1) * fT10 * (-ST10)^(-1) + (-log(ST20))^(theta-1)  * fT20 * (-ST20)^(-1))/((-log(ST10))^theta + (-log(ST20))^theta)
    fstar1 <- Sstar1 * log(Sstar1) * ((-log(ST11))^(theta-1) * fT11 * (-ST11)^(-1) + (-log(ST21))^(theta-1)  * fT21 * (-ST21)^(-1))/((-log(ST11))^theta + (-log(ST21))^theta)
  }
  
  ##-- Hazards and hazard ratio for the composite
  Lstar0 <- (fstar0/Sstar0)
  Lstar1 <- (fstar1/Sstar1)
  HRstar <- (Lstar1/Lstar0)
  
  ##-- Summary measures for the HR* (see Schempfer 2009)
  HRstar_int <- rowMeans(cbind(HRstar[-1],rev(rev(HRstar)[-1]))) # Mean of HRs for each interval
  fstar0_int <- rowMeans(cbind(fstar0[-1],rev(rev(fstar0)[-1]))) # Mean of fstar0 for each interval
  fstar1_int <- rowMeans(cbind(fstar1[-1],rev(rev(fstar1)[-1]))) # Mean of fstar1 for each interval
  Lstar0_int <- rowMeans(cbind(Lstar0[-1],rev(rev(Lstar0)[-1]))) # Mean of Lstar0 for each interval
  Lstar1_int <- rowMeans(cbind(Lstar1[-1],rev(rev(Lstar1)[-1]))) # Mean of Lstar1 for each interval
  
  ##-- Unweighted treatment effect measures (not returned)
  nHR <- (HR_e1+HR_e2)/2                                         # naive HR (arithmetic mean of the 2 average HR)
  mHR <- mean(HRstar)                                            # Unweighted mean of average HR
  
  ##-- Weighted sAHR (not returned)
  sAHR_0 <- sum(HRstar_int*fstar0_int)/sum(fstar0_int)           # sHR "_0" indicates that f_0 is used instead of f_1 
  
  ##-- gAHR weighted by f0 and by f0+f1
  gAHR_0 <- exp(sum(log(HRstar_int)*fstar0_int)/sum(fstar0_int))                         # gHR "_0" indicates that f_0 is used instead of f_0 + f_1
  gAHR   <- exp(sum(log(HRstar_int)*(fstar0_int+fstar1_int))/sum(fstar0_int+fstar1_int)) # Alternative gAHR weighted by f_0 + f_1
  
  ##-- AHR weighted by f0 and by f0+f1
  AHR_0_num <- sum(Lstar1_int/(Lstar0_int + Lstar1_int)*fstar0_int)/sum(fstar0_int)
  AHR_0_den <- sum(Lstar0_int/(Lstar0_int + Lstar1_int)*fstar0_int)/sum(fstar0_int)
  AHR_0 <- AHR_0_num/AHR_0_den                                                           # AHR "_0" indicates that f_0 is used instead of f_1

  AHR_num <- sum(Lstar1_int/(Lstar0_int + Lstar1_int)*(fstar0_int+fstar1_int))/sum(fstar0_int+fstar1_int)
  AHR_den <- sum(Lstar0_int/(Lstar0_int + Lstar1_int)*(fstar0_int+fstar1_int))/sum(fstar0_int+fstar1_int)
  AHR <- AHR_num/AHR_den                                                                 # Alternative AHR weighted by f_0 + f_1
    
  ##-- RMST (Restricted Mean Survival Time)
  RMST_0 <- integrate(Sstar, dist1=T1pdist,dist2=T2pdist,param1=T10param,param2=T20param,dist_biv= distribution0, lower=0,upper=1,subdivisions=subdivisions)$value
  RMST_1 <- integrate(Sstar, dist1=T1pdist,dist2=T2pdist,param1=T11param,param2=T21param,dist_biv= distribution1, lower=0,upper=1,subdivisions=subdivisions)$value
  
  ##-- Probability of event during follow_up
  pstar_0 <- 1 - Sstar(1,dist1=T1pdist,dist2=T2pdist,param1=T10param,param2=T20param,dist_biv= distribution0)
  pstar_1 <- 1 - Sstar(1,dist1=T1pdist,dist2=T2pdist,param1=T11param,param2=T21param,dist_biv= distribution1)
  # Alternative calculation
  # pstar0 <- pMvdc(c(1,Inf), distribution0) + pMvdc(c(Inf,1), distribution0) - pMvdc(c(1,1), distribution0)
  # pstar1 <- pMvdc(c(1,Inf), distribution1) + pMvdc(c(Inf,1), distribution1) - pMvdc(c(1,1), distribution1)
  
  ##-- Median
  limits <- c(0,10)                                                                         # The first and the last values must be in opposite signs for the function
  Med_0 <- uniroot(Sstar_func_perc, interval=limits,extendInt="yes", perc=0.5,dist1=T1pdist,dist2=T2pdist,param1=T10param,param2=T20param,dist_biv= distribution0)$root  # Find the root (value which equals the function to zero). extendInt="yes": 'interval limits' is extended automatically if necessary
  Med_1 <- uniroot(Sstar_func_perc, interval=limits,extendInt="yes", perc=0.5,dist1=T1pdist,dist2=T2pdist,param1=T11param,param2=T21param,dist_biv= distribution1)$root  # Find the root (value which equals the function to zero). extendInt="yes": 'interval limits' is extended automatically if necessary
  
  ##-- Probabilities p11,p21 (Not returned)
  # p11 <- 1-exp(-(1/b11)^beta_e1)
  # p21 <- 1-exp(-(1/b21)^beta_e2)
  
  if(plot_HR){
    dd <- data.frame(t=t, HRstar=HRstar)
    ymin <- min(HRstar,0.5)
    ymax <- max(HRstar,1)
    gg1 <- ggplot(dd, aes(x=t,y=HRstar)) + geom_line(color='blue',size=2) + 
      ylim(ymin,ymax) + 
      ggtitle('HR(t) of the composite endpoint') + 
      xlab('Proportion of follow-up time') + ylab('HR*')
    print(gg1)
  }
  
  return(list(effect_size=list('gAHR' = round(gAHR,4),
                               'AHR' = round(AHR,4),
                               'RMST_ratio' = round(RMST_1/RMST_0,4),
                               'Median_Ratio' = round(Med_1/Med_0,4)),
              measures_by_group=list('pstar'= c('Reference'=pstar_0,'Treated'=pstar_1),
                                     'RMST'= c('Reference'=RMST_0,'Treated'=RMST_1),
                                     'Median'= c('Reference'=Med_0,'Treated'=Med_1))))
}
