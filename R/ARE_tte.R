#' ARE method for composite time to event endpoints
#'
#' @description The composite endpoint is assumed to be a time to event endpoint formed by a combination of two events (E1 and E2). We assume that the endpoint 1 is more relevant for the clinical question than endpoint 2. 
#' This function calculates the ARE (Assymptotic Relative Efficiency) method for time to event endpoints. The method quantifies the differences in efficiency of using the composite or the relevant as primary endpoint to lead the trial and, moreover, provides a decision rule to choose the primary endpoint. If the ARE is larger than 1, the composite endpoint may be considered the best option as primary endpoint. Otherwise, the relevant endpoint is preferred. 
#' 
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
#' @param subdivisions integer parameter greater than or equal to 10. Number of points used to plot the ARE according to correlation. The default is 50. Ignored if plot_res=FALSE and plot_store=FALSE. 
#' @param plot_res logical indicating if the ARE according to the correlation should be displayed. The default is FALSE
#' @param plot_store logical indicating if the plot of ARE according to the correlation is stored for future customization. The default is FALSE
#' 
#' @rawNamespace import(copula, except = c(profile,coef,logLik,confint))
#' @import utils
#' @export 
#'
#' @return Returns the ARE value along with the fixed correlation. If the ARE 
#' value is larger than 1 then the composite endpoint is preferred over the 
#' relevant endpoint. Otherwise, the endpoint 1 is preferred as the primary 
#' endpoint of the study. In addition, if \code{plot_store=TRUE} an object of 
#' class \code{ggplot} with the ARE according to the correlation is stored in the output.
#'
#' @details Some parameters might be difficult to anticipate, especially the shape parameters of Weibull distributions and those referred to the relationship between the marginal distributions. 
#' For the shape parameters (beta_e1, beta_e2) of the Weibull distribution, we recommend to use \eqn{\beta_j=0.5}, \eqn{\beta_j=1} or \eqn{\beta_j=2} if a decreasing, constant or increasing rates over time are expected, respectively.
#' For the correlation (rho) between both endpoints, generally a positive value is expected as it has no sense to design an study with two endpoints negatively correlated. We recommend to use \eqn{\rho=0.1}, \eqn{\rho=0.3} or \eqn{\rho=0.5} for weak, mild and moderate correlations, respectively.
#' For the type of correlation (rho_type), although two different type of correlations are implemented, we recommend the use of the Spearman's correlation.
#' In any case, if no information is available on these parameters, we recommend to use the default values provided by the function.
#'
#'
#' @references Gomez Melis, G. and Lagakos, S.W. (2013). Statistical considerations when using a composite endpoint for comparing treatment groups. Statistics in Medicine. Vol 32(5), pp. 719-738. https://doi.org/10.1002/sim.5547
#'
#' @examples
#' # ARE for a specific study where the composite endpoint is recommended
#' ARE_tte(p0_e1=0.1, p0_e2=0.1, HR_e1=0.9, HR_e2=0.8, beta_e1 = 1, beta_e2 = 1, 
#' case=1, copula = "Frank", rho = 0.3, rho_type = "Spearman")
#' # ARE for a specific study where the composite endpoint is not recommended
#' ARE_tte(p0_e1=0.1, p0_e2=0.05, HR_e1=0.6, HR_e2=0.8, beta_e1 = 1, beta_e2 = 1, 
#' case=1, copula = "Frank", rho = 0.3, rho_type = "Spearman")


ARE_tte <- function(p0_e1, p0_e2, HR_e1, HR_e2, beta_e1=1, beta_e2=1, 
                    case, copula = 'Frank', rho=0.3, rho_type='Spearman', 
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
    stop("The case (case) must be a number in {1,2,3,4}. See ?ARE_tte")
  }else if(!copula %in% c('Frank','Gumbel','Clayton')){
    stop("The copula (copula) must be one of 'Frank','Gumbel' or 'Clayton'")
  }else if(rho < -1 || rho > 1){
    stop("The correlation (rho) must be a number between -1 and 1")
  }else if(!rho_type %in% c('Spearman','Kendall')){
    stop("The correlation type (rho_type) must be one of 'Spearman' or 'Kendall'")
  }else if(!(is.numeric(subdivisions) && subdivisions>=10)){
    stop("The number of subdivisions must be an integer greater than or equal to 10")    
  }else if(!is.logical(plot_res)){
    stop("The parameter plot_res must be logical")
  }else if(!is.logical(plot_store)){
    stop("The parameter plot_store must be logical")      
  }else if(case==4 && p0_e1 + p0_e2 > 1){
    stop("The sum of the proportions of observed events in both endpoints in case 4 must be lower than 1")
  }
  
  # Values of rho where to calculate ARE
  rho_sel <- rho
  if(plot_res | plot_store){
    rho_seq <- unique(c(rho,seq(0.01,0.98,length=subdivisions)))
  }else{
    rho_seq <- rho
  }
  
  # Storage
  ARE_array <- c()
  
  # Calculate ARE for each rho
  pb = txtProgressBar(min = 0, max = length(rho_seq), initial = 0)
  for(rho in rho_seq){
    setTxtProgressBar(pb,which(rho_seq==rho))
    
    # Copula
    copula0 <- CopulaSelection(copula=copula,rho=rho,rho_type=rho_type)
    which.copula <- copula0[[1]]
    theta <- copula0[[2]]   
    
    # Marginal distribution and parameters
    MarginSelec <- MarginalsSelection(beta_e1,beta_e2,HR_e1,HR_e2,p0_e1,p0_e2,case,rho,theta,copula)
    T1dist   <- MarginSelec[[1]]
    T2dist   <- MarginSelec[[2]]
    T1pdist  <- MarginSelec[[3]]
    T2pdist  <- MarginSelec[[4]]
    T10param <- MarginSelec[[5]]
    T20param <- MarginSelec[[6]]
    T11param <- MarginSelec[[7]]
    T21param <- MarginSelec[[8]]
    
    # Bivariate distribution in control and treatment groups
    distribution0 <- mvdc(copula = which.copula, margins = c(T1dist, T2dist), paramMargins = list(T10param, T20param))
    distribution1 <- mvdc(copula = which.copula, margins = c(T1dist, T2dist), paramMargins = list(T11param, T21param))
    
    if(case==1|case==3) {
      
      # Inside the integral in the numerator --> Does not work in a function apart (I do not why)
      inside_integral <- function(t){
        
        Sstar0 <- Sstar(x = t,dist1 = T1pdist, dist2 = T2pdist, param1 = T10param, param2 = T20param, dist_biv = distribution0)
        Sstar1 <- Sstar(x = t,dist1 = T1pdist, dist2 = T2pdist, param1 = T11param, param2 = T21param, dist_biv = distribution1)
        
        fstar0 <- (-grad(Sstar,x = t, dist1 = T1pdist,dist2 = T2pdist, param1 = T10param, param2 = T20param,dist_biv = distribution0))
        fstar1 <- (-grad(Sstar,x = t, dist1 = T1pdist,dist2 = T2pdist, param1 = T11param, param2 = T21param,dist_biv = distribution1))
        
        Lstar0 <- (fstar0/Sstar0)
        Lstar1 <- (fstar1/Sstar1)
        
        HRstar <- (Lstar1/Lstar0)
        
        logHRstar <- log(HRstar)
        
        ##-- Numerical issues
        fstar0[fstar0<0] <- 0
        logHRstar[is.na(logHRstar) | logHRstar==Inf | logHRstar== -Inf] <- 0
        
        return(logHRstar*fstar0)
      }
      
      # Numerator
      integral  <- integrate(inside_integral,lower=0,upper=1,subdivisions=10000,stop.on.error = FALSE) 
      numerator <- (integral$value)^2
      
      # Denominator
      Sstar0_1 <- Sstar(x=1,dist1=T1pdist,dist2=T2pdist,param1=T10param,param2=T20param,dist_biv= distribution0) 
      ST10_1 <- 1-do.call(T1pdist,c(q=1,T10param)) 
      denominator <- ((log(HR_e1))^2)*(1-Sstar0_1)*(1-ST10_1)
      
      # ARE value
      AREstarT <- numerator/denominator
      
      # If the integral is not computed, we assign a missing value
      if(integral$message!="OK") {AREstarT <- NA}
      
    } else if(case==2|case==4) {
      
      b10 <- T10param$scale
      b20 <- T20param$scale
      b11 <- T11param$scale
      b21 <- T21param$scale
      
      
      # Only marginal Weibull distributions for fT10, fT20, ST10, ST20.
      # fT10 <- function(t) dweibull(x=t,beta_e1,b10)
      # ST10 <- function(t) 1-pweibull(q=t,beta_e1,b10)
      # 
      # fT20 <- function(t) dweibull(x=t,beta_e2,b20)
      # ST20 <- function(t) 1-pweibull(q=t,beta_e2,b20)
      fT0 <- function(t,beta,b) dweibull(x=t,beta,b)         # density function
      ST0 <- function(t,beta,b) 1-pweibull(q=t,beta,b)       # survival function
      
      
      
      # Sstar0 and fstar0 for any copula
      Sstar0 <- function(t) Sstar(x=t,dist1=T1pdist,dist2=T2pdist,param1=T10param,param2=T20param,dist_biv= distribution0)
      fstar0 <- function(t) (-grad(Sstar,x=t,dist1=T1pdist,dist2=T2pdist,param1=T10param,param2=T20param,dist_biv= distribution0))
      
      ##-- Cause specific (CS) hazards
      # CS endpoint 1
      aux21 <- function(t,y) theta*exp(-theta*(ST0(t,beta_e1,b10)+y))*(1-exp(-theta))/(exp(-theta)-exp(-theta*ST0(t,beta_e1,b10))-exp(-theta*y)+exp(-theta*(ST0(t,beta_e1,b10)+y)))^2
      aux22 <- function(u) {integrate(aux21,0, ST0(u,beta_e2,b20),t=u,subdivisions=10000)$value} # t=u indicates that we are derivating respect to the other variable in aux21(t,y). That is, respect to y.
      lambdaC10 <- function(t) aux22(t)*fT0(t,beta_e1,b10)/Sstar0(t)
      lambdaC11 <- function(t) HR_e1*lambdaC10(t)
      
      # CS endpoint 2
      aux23 <- function(x,t) theta*exp(-theta*(x+ST0(t,beta_e2,b20)))*(1-exp(-theta))/(exp(-theta)-exp(-theta*x)-exp(-theta*ST0(t,beta_e2,b20))+exp(-theta*(x+ST0(t,beta_e2,b20))))^2
      aux24 <- Vectorize(function(u){integrate(aux23,0,ST0(u,beta_e1,b10),t=u,subdivisions=10000)$value}) #t=u indicates that we are derivating respect to the other variable in aux23(x,t). That is, respect to x.
      lambdaC20 <- function(t) aux24(t)*fT0(t,beta_e2,b20)/Sstar0(t)
      lambdaC21 <- function(t) HR_e2*lambdaC20(t)
      
      
      # Evaluation of LambdaC20 before computation 
      # IT MAY FAIL IN CASES 2/4 FOR BETAS = 0.5 BECAUSE IT IS NOT ALWAYS EVALUABLE AT T=0
      # WHENEVER LambdaC20 FAILS, WE INCREASE THE LOWER LIMIT OF INTEGRATION
      LambdaC20_check <- tryCatch(LambdaC20 <- function(t) integrate(lambdaC20,lower=0,upper=t,subdivisions=10000)$value,error = function(e) e)
      lower_LambdaC20 <- 0
      while(inherits(LambdaC20_check, "error")=="TRUE" ){
        lower_LambdaC20 <- lower_LambdaC20 + 0.001
        LambdaC20_check <- tryCatch(LambdaC20 <- function(t) integrate(lambdaC20,lower=0+lower_LambdaC20,upper=t,subdivisions=10000)$value,error = function(e) e)
        print(LambdaC20_check)
      }
      LambdaC20 <- Vectorize(function(t) integrate(lambdaC20,lower=0+lower_LambdaC20,upper=t,subdivisions=10000)$value)
      
      
      # Computation of the hazards for both groups
      Lstar0 <- function(t) lambdaC10(t) + lambdaC20(t)
      Lstar1 <- function(t) lambdaC11(t) + lambdaC21(t)
      
      # Computation of HRstar
      HRstar <- function(t) Lstar1(t)/Lstar0(t)
      logHRstar <- function(t) log(Lstar1(t)/Lstar0(t))
      
      temp3 <- Vectorize(function(t) logHRstar(t)*fstar0(t))
      temp4_check <- tryCatch(temp4 <- integrate(temp3,0,1,subdivisions=10000)$value,error = function(e) e)
      lower_temp4 <- 0
      while(inherits(temp4_check, "error")=="TRUE" & lower_temp4 < 1){ # add "& lower_temp4 < 1"
        lower_temp4 <- lower_temp4 + 0.001
        temp4_check<-tryCatch(temp4 <- integrate(temp3, lower_temp4, 1, subdivisions=10000)$value,error = function(e) e)
      }
      temp4 <- integrate(temp3,0+lower_temp4, 1, subdivisions=10000)$value
      numerator <- (temp4)^2
      
      #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      ## Computation of PROBT1UNC
      PROBT1UNC_temp_num <- function(t) exp(-HR_e2*LambdaC20(t)) * Sstar0(t) * lambdaC10(t)
      PROBT1UNC_temp_den <- function(t) 1/2 * (exp(-LambdaC20(t)) + exp(-HR_e2 * LambdaC20(t)))
      PROBT1UNC_temp <- Vectorize(function(t){PROBT1UNC_temp_num(t)/PROBT1UNC_temp_den(t)})
      PROBT1UNC_int_check <- tryCatch(integrate(PROBT1UNC_temp,lower=0, upper=1,subdivisions=10000)$value, error = function(e) e)
      lower_LambdaC20 <- 0
      while(inherits(PROBT1UNC_int_check, "error")=="TRUE"){
        lower_LambdaC20 <- lower_LambdaC20+0.001
        PROBT1UNC_int_check <- tryCatch(integrate(PROBT1UNC_temp,lower=lower_LambdaC20, upper=1,theta=theta,HR2=HR_e2,subdivisions=10000)$value, error = function(e) e)
        print(lower_LambdaC20)
      }
      lower_PROBT1UNC_int <- lower_LambdaC20
      PROBT1UNC_int <- integrate(PROBT1UNC_temp,lower=lower_PROBT1UNC_int, upper=1,subdivisions=10000)$value
      ## End of computation of PROBT1UNC
      #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      
      ############################################
      # ARE VALUE:
      AREstarT <- numerator/((log(HR_e1)^2) * PROBT1UNC_int * (1-Sstar0(1)))
      AREstarT 
    }
    
    ARE_array <- c(ARE_array,AREstarT)
  }
  close(pb)
  
  # if(ARE_array[1]>1){
  #   cat("The use of the composite endpoint as primary endpoint is recommended over the use of the relevant endpoint since ARE =",
  #       formatC(ARE_array[1],digits = 3,big.mark = ','),"> 1.\n")
  # }else{
  #   cat("The use of the first endpoint as primary endpoint is recommended over the use of the composite endpoint since ARE =",
  #       formatC(ARE_array[1],digits = 3,big.mark = ','),"< 1.\n")
  # }
  
  if(plot_res | plot_store){
    ARE <- NULL          # To avoid the note: "no visible binding for global variable 'NULL'"
    dd <- data.frame(rho=rho_seq, ARE=ARE_array)
    gg1 <- ggplot(dd,aes(x=rho,y=ARE)) + 
      geom_line(color='darkblue',size=1.3) +
      xlab(expression(rho)) + ylab('ARE CE') +
      scale_y_log10(limits=c(1/max(ARE_array),max(ARE_array))) +
      geom_hline(yintercept=1,linetype='dashed')
  }
  
  return_object <- list(ARE=ARE_array[1],
                        rho=rho_sel,
                        gg_object=NA)
  
  ## Print graphic
  if(plot_res) print(gg1)
  
  ## Store plot in the output
  if(plot_store) return_object$gg_object <- gg1

  
  return(invisible(return_object))
}


