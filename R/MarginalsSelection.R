#' Marginal distributions for both components
#' 
#' @description   Returns the distribution and parameters of the marginal Weibull distributions
#' 
#' @param beta1	  Shape parameter for a Weibull law for the relevant event
#' @param beta2   Shape parameter for a Weibull law for the additional event 
#' @param HR1     Hazard Ratio for a Weibull law for the relevant event
#' @param HR2     Hazard Ratio for a Weibull law for the additional event
#' @param p1      Proportion of the relevant event expected in group zero
#' @param p2      Proportion of the additional event expected in group zero
#' @param case    Censoring case 
#' @param theta   Dependence parameter for the bivariate distribution in control group
#' @param copula  Copula to use
#'
#' @export 
#' @keywords internal 
#'
#'

MarginalsSelection <- function(beta1,beta2,HR1,HR2,p1,p2,case,rho,theta,copula='Frank') 
{

  dcopula <- get(paste0('d',copula))
  
  ## -- Case 1 --------------------------------------------------------
  if(case==1) {
    b10 <- 1/((-log(1-p1))^(1/beta1))
    b20 <- 1/((-log(1-p2))^(1/beta2))
  
  ## -- Case 2 --------------------------------------------------------
  }else if (case==2){  
    
    Fb10 <- function(b10,p1){
      integral<-integrate(function(u) {
        sapply(u, function(u) {
          integrate(function(v){dcopula(u,v,theta)}, lower=1-exp((b10*(-log(1-u))^(1/beta1))^beta2*log(1-p2)), upper=1)$value
        })
      }, lower=0 , upper=1-exp(-1/b10^beta1))$value
      return(integral-p1) 
    }
    limits <- c(0.00001,10000)                                         # The first and the last values must be in opposite signs for the function
    b10 <- try(uniroot(Fb10, interval=limits,p1=p1)$root,silent=TRUE)  # Find the root (value which equals the function zero)
    b20 <- 1/(-log(1-p2))^(1/beta2)
    if(inherits(b10,'try-error')){
      dcopula <- dFrank
      b10 <- uniroot(Fb10, interval=limits,p1=p1)$root
      dcopula <- get(paste0('d',copula))
      limits <- c(0.8,1.2)*b10 
      b10 <- uniroot(Fb10, interval=limits,p1=p1)$root
    }
  
  ## -- Case 3 --------------------------------------------------------
  } else if (case==3) {
    
    b10 <- 1/((-log(1-p1))^(1/beta1))
    
    Fb20 <- function(b20,p2) {
      integral<-integrate(function(v) {
        sapply(v,function(v) { 
          integrate(function(u){dcopula(u,v,theta)},lower=1-exp((b20*(-log(1-v))^(1/beta2))^beta1*log(1-p1)),upper=1)$value
        })
      }, 
      lower=0 , upper=1-exp(-1/b20^beta2))$value
      return(integral-p2)
    }
    limits <- c(0.00001,10000) 
    b20 <- try(uniroot(Fb20, interval=limits,p2=p2)$root,silent=TRUE)
    if(inherits(b20,'try-error')){
      # First attemp: to approximate for limits based on frank copula
      dcopula <- dFrank
      b20 <- uniroot(Fb20, interval=limits,p2=p2)$root
      dcopula <- get(paste0('d',copula))
      limits <- c(0.5,2)*b20
      b20 <- try(uniroot(Fb20, interval=limits,p2=p2)$root,silent=TRUE)
      
      # Second attemp: search for factible limits
      if(inherits(b20,'try-error')){
        FFB20 <- c()
        VALUES <- unique(c(seq(0.1,10,0.1),seq(10,100,1),seq(100,1000,10),seq(1000,10000,100)))
        for (i in 1:length(VALUES)) {aux <- try(Fb20(VALUES[i],p2),silent=TRUE); FFB20[i] <- ifelse(class(aux)=='try-error',NA,aux)}
        POS_LIMITS <- c(neg=NA,pos=NA)
        POS_MIN <- NA
        VAL_MIN <- 1000
        j <- 1
        while(any(is.na(POS_LIMITS)) & j<=length(FFB20)){
          if(is.na(FFB20[j])){
            POS_LIMITS <- c(neg=NA,pos=NA)
          }else if(FFB20[j]<=0){
            POS_LIMITS[1] <- j
            if(abs(FFB20[j])<VAL_MIN){VAL_MIN <- abs(FFB20[j]);POS_MIN <- j}
          }else{
            POS_LIMITS[2] <- j
            if(abs(FFB20[j])<VAL_MIN){VAL_MIN <- abs(FFB20[j]);POS_MIN <- j}
          }
          
          j <- j + 1
        }
        if(!any(is.na(POS_LIMITS))){
          limits <- VALUES[POS_LIMITS]
          b20 <- uniroot(Fb20, interval=limits,p2=p2)$root
        }else{
          b20 <- VALUES[POS_MIN]
          cat('Imposed value to b20:',b20,'Parameters:',beta1,beta2,HR1,HR2,p1,p2,case,theta,copula,'\n',file = "log_marginalselection.txt",append=TRUE)
        }
        
      }
    }
  
  ## -- Case 4 --------------------------------------------------------
  } else if (case==4) {
    
    x <- c(NA,NA)
    
    # To compute b10
    Fb10_case4 <- function(b10,b20,p1){
      x <- c(b10,b20)
      integral<-integrate(function(u) {
        sapply(u, function(u) {
          integrate(function(v){dcopula(u,v,theta)}, lower=1-exp(-(x[1]*(-log(1-u))^(1/beta1)/x[2])^beta2), upper=1)$value
        })
      }, lower= 0 , upper=1-exp(-1/x[1]^beta1))$value
      return(integral-p1)
    }
    
    # To compute b20
    Fb20_case4 <- function(b10,b20,p2) {
      x <- c(b10,b20)
      integral<-integrate(function(v) {
        sapply(v,function(v) {
          integrate(function(u){dcopula(u,v,theta)},lower=1-exp(-(x[2]*(-log(1-v))^(1/beta2)/x[1])^beta1), upper=1)$value
        })
      },
      lower= 0, upper=1-exp(-1/x[2]^beta2))$value
      return(integral-p2)
    }
    
    model <- function(x){
      c(Fb10_case4(x[1],x[2],p1), Fb20_case4(x[1],x[2],p2))
    }
    
    sol <- multiroot(f = model, start = c(1,1))

    sol<-as.data.frame(sol[1])
    b10<-sol[1,]
    b20<-sol[2,]
    
  }
  
  # Scale parameters for group 1 b11,b21
  b11 <- b10/HR1^(1/beta1)
  b21 <- b20/HR2^(1/beta2)
  
  # Probabilities p11,p21
  #p11 <- 1-exp(-(1/b11)^beta1)
  #p21 <- 1-exp(-(1/b21)^beta2)
  if(case==1){
    p11 <- 1-exp(-(1/b11)^beta1)
    p21 <- 1-exp(-(1/b21)^beta2)
  }else if(case==2){
    p11 <- get_prob1(beta1,beta2,b11,b21,case,rho,copula,endpoint=1)
    p21 <- 1-exp(-(1/b21)^beta2)
  }else if(case==3){
    p11 <- 1-exp(-(1/b11)^beta1)
    # cat(beta1,beta2,b11,b21,case,rho,copula,"\n")
    p21 <- get_prob1(beta1,beta2,b11,b21,case,rho,copula,endpoint=2) # theta or rho?
  }else if(case==4){
    p11 <- get_prob1(beta1,beta2,b11,b21,case,rho,copula,endpoint=1)
    p21 <- get_prob1(beta1,beta2,b11,b21,case,rho,copula,endpoint=2) # theta or rho?
  } 

  T1dist<-"weibull"
  T2dist<-"weibull"
  
  T1pdist<-pweibull
  T2pdist<-pweibull
  
  T10param <- list(shape = beta1, scale = b10)
  T20param <- list(shape = beta2, scale = b20)
  
  T11param <- list(shape = beta1, scale = b11)
  T21param <- list(shape = beta2, scale = b21)
  
  return(list(T1dist,T2dist,T1pdist,T2pdist,T10param,T20param,T11param,T21param,p11,p21))	
}