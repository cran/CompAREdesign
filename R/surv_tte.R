#' Survival function for composite time-to-event endpoints
#'
#' @description It provides the survival function for time-to-event composite outcomes. 
#' The composite endpoint is assumed to be a time-to-event endpoint formed by a combination of two events (E1 and E2).
#' The effect size is calculated on the basis of anticipated information on the composite components and the correlation between them.
#' Marginal distributions are assumed weibull for both endpoints.
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
#' @param followup_time numeric parameter indicating the maximum follow up time (in any unit). Default is 1.
#' @param plot_res logical indicating if the survival curves should be displayed. The default is TRUE
#' @param plot_store logical indicating if the plot of the survival curve for composite endpoint is stored for future customization. The default is FALSE
#' 
#' @import ggplot2
#' @import rootSolve
#' @rawNamespace import(copula, except = c(profile,coef,logLik,confint))
#' @rawNamespace import(numDeriv, except = hessian)
#' 
#' @export 
#'
#' @return For each group, if  \code{plot_res=TRUE}, the function returns a plot 
#' of the survival functions for composite endpoint as well as the plots of the 
#' survival function for each component.
#'   
#' @details Some parameters might be difficult to anticipate, especially the shape parameters of Weibull distributions and those referred to the relationship between the marginal distributions. 
#' For the shape parameters (beta_e1, beta_e2) of the Weibull distribution, we recommend to use \eqn{\beta_j=0.5}, \eqn{\beta_j=1} or \eqn{\beta_j=2} if a decreasing, constant or increasing rates over time are expected, respectively.
#' For the correlation (rho) between both endpoints, generally a positive value is expected as it has no sense to design an study with two endpoints negatively correlated. We recommend to use \eqn{\rho=0.1}, \eqn{\rho=0.3} or \eqn{\rho=0.5} for weak, mild and moderate correlations, respectively.
#' For the type of correlation (rho_type), although two different type of correlations are implemented, we recommend the use of the Spearman's correlation.
#' In any case, if no information is available on these parameters, we recommend to use the default values provided by the function.
#' 
#'
#'

surv_tte <- function(p0_e1, p0_e2, HR_e1, HR_e2, beta_e1=1, beta_e2=1, case, 
                     copula = 'Frank', rho=0.3, rho_type='Spearman',followup_time=1, 
                     plot_res=TRUE, plot_store=FALSE){
  
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
  }else if(!is.logical(plot_res)){
     stop("The parameter plot_res must be logical")
  }else if(!is.logical(plot_store)){
    stop("The parameter plot_store must be logical")    
  }else if(case==4 && p0_e1 + p0_e2 > 1){
    stop("The sum of the proportions of observed events in both endpoints in case 4 must be lower than 1")
  }
  

  
  ##-- This would be calculated only one time
  ## To change: constant plot independent of the time of follow-up (change x tick marks but not the plot itself)
  copula0 <- CopulaSelection(copula,rho=rho,rho_type)
  theta <- copula0[[2]] 
  MS <- MarginalsSelection(beta1=beta_e1,beta2=beta_e2,HR1=HR_e1,HR2=HR_e2,p1=p0_e1,p2=p0_e2,
                           case=case,rho=rho,theta=theta,copula=copula)
  
  ##################################################
  # Survival function
  ##################################################

  ##-- Survival function
  sweibull <- function(...) 1 - pweibull(...)
  
  ##-- Time (t) values to assess the survival function
  t <- c(0.0001,seq(0.001,1,length.out = 1000)) # subdivisions=1000
  
  ##-- Scale parameters
  b0_e1 <- MS[[5]][[2]]  # Scale parameter endpoint 1, control arm
  b0_e2 <- MS[[6]][[2]]  # Scale parameter endpoint 2, treated arm
  b1_e1 <- MS[[7]][[2]]  # Scale parameter endpoint 1, control arm
  b1_e2 <- MS[[8]][[2]]  # Scale parameter endpoint 2, treated arm
    
  ##-- Survival for both endpoints
  ST10 <- exp(-(t/b0_e1)^beta_e1)
  ST20 <- exp(-(t/b0_e2)^beta_e2)
  ST11 <- exp(-(t/b1_e1)^beta_e1)
  ST21 <- exp(-(t/b1_e2)^beta_e2)
  
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
  
  ##################################################
  # Plots
  ##################################################
  if(plot_res | plot_store){
    xmax <- max(1,as.numeric(followup_time),na.rm=TRUE)
    
    theme.plot <- theme(legend.position="bottom",
                        #legend.text=element_text(size=15,face='bold'),
                        legend.title =element_blank(),
                        legend.key.size = unit(0.8, "cm"))
                        #axis.text=element_text(size=14,hjust=0.5, face='bold'),
                        #axis.title.x=element_text(size=15,face="bold"),
                        #axis.title = element_text(size=10))
                        # plot.title = element_text(face='bold'))
    
    # Data for plot
    t_plot <- c(0.0001, seq(0.04,0.99,0.05))   # x points where to calculate HR*
    n_t_plot <- length(t_plot)                 # number of x points to calculate survival function
    f_time <- as.numeric(followup_time)        # follow_up_time  
    x <- NULL                                  # To avoid the note: "no visible binding for global variable 'x'"
    
    ## Endpoint 1
    gg1  <- ggplot(data = data.frame(x = 0), mapping = aes(x = x)) +
      stat_function(fun = sweibull,args = list(shape=MS[[5]]$shape,scale=MS[[5]]$scale),aes(color='lightblue'),size=1.3,linetype='longdash') +
      stat_function(fun = sweibull,args = list(shape=MS[[7]]$shape,scale=MS[[7]]$scale),aes(color='darkcyan'),size=1.3,linetype='longdash')  +
      xlab('Time') + ylab('Survival E1') + 
      scale_y_continuous(limits=c(0,1),minor_breaks=NULL,expand=c(0,0)) + 
      scale_x_continuous(limits=c(0,1),
                         breaks=pretty(0:1*f_time)/f_time,
                         labels=pretty(0:1*f_time),expand=c(0,0.01)) +
      scale_color_identity(name = "Group",
                           breaks = c('darkcyan','lightblue'), 
                           labels = c("Treated", "Control"),
                           guide = "legend") + theme.plot
    
    ## Endpoint 2
    gg2  <- ggplot(data = data.frame(x = 0), mapping = aes(x = x)) +
      stat_function(fun = sweibull,args = list(shape=MS[[6]]$shape,scale=MS[[6]]$scale),aes(color='lightblue'),size=1.3,linetype='longdash') +
      stat_function(fun = sweibull,args = list(shape=MS[[8]]$shape,scale=MS[[8]]$scale),aes(color='darkcyan'),size=1.3,linetype='longdash')  +
      xlab('Time') + ylab('Survival E2') + 
      scale_y_continuous(limits=c(0,1),minor_breaks=NULL,expand=c(0,0)) +
      scale_x_continuous(limits=c(0,1),breaks=pretty(0:1*f_time)/f_time,labels=pretty(0:1*f_time),expand=c(0,0.01)) +
      scale_color_identity(name = "Group",
                           breaks = c('darkcyan','lightblue'),
                           labels = c("Treated", "Control"),
                           guide = "legend") + theme.plot
    
    
    ## Composite endpoint
    gg3 <- ggplot() +
      geom_line(mapping=aes(x=t,y=Sstar0),linetype ='longdash',size=1.3,color='lightblue') +
      geom_line(mapping=aes(x=t,y=Sstar1),linetype ='longdash',size=1.3,color='darkcyan') +
      scale_y_continuous(limits=c(0,1),minor_breaks=NULL,expand=c(0,0)) +
      scale_x_continuous(limits=c(0,1),breaks=pretty(0:1*f_time)/f_time,
                         labels=pretty(0:1*f_time),expand=c(0,0.01)) +
      xlab('Time') + ylab('Survival CE') + 
      theme.plot
    
    gg_all <- ggarrange(gg1,gg2,gg3,nrow=1,ncol=3,common.legend = TRUE)
  }

  return_object <- list(gg_object=NA)
  
  ## Print graphic
  if(plot_res) print(gg_all)
  
  ## Store plot in the output
  if(plot_store) return_object$gg_object <- gg3
  
  return(invisible(return_object))
}
