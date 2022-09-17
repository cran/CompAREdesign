#' Simulation 1-arm
#' 
#' @description Simulation of a trial with a composite endpoint
#' 
#' @param p_e1 numeric parameter, probability of the event E1
#' @param p_e2 numeric parameter, probability of the event E2
#' @param p_ce probability of the composite endpoint
#' @param samplesize sample size per arm
#'
#' @export 
#' @keywords internal 
#'
#' @return simulated data
#'
#' @examples
#' data=f_sim(samplesize=100,p_e1=0.1,p_e2=0.5,p_ce=0.6)
#' head(data)

f_sim <- function(samplesize,p_e1,p_e2,p_ce){ 
  
  # 2x2 table
  s1_group = p_e1 + p_e2 - p_ce
  s2_group = ifelse(p_ce-p_e2>0, p_ce-p_e2, 0)#p_ce-p_e2  
  s3_group = ifelse(p_ce-p_e1>0, p_ce-p_e1, 0)#p_ce-p_e1  
  s4_group = 1- p_ce  
  
  data = rmultinom(n=1,size=round(samplesize),prob=c(s1_group,s2_group,s3_group,s4_group)) 
  
  s1=data.frame(binary_e1=rep(1,data[1]),binary_e2=rep(1,data[1]))
  s2=data.frame(binary_e1=rep(1,data[2]),binary_e2=rep(0,data[2]))
  s3=data.frame(binary_e1=rep(0,data[3]),binary_e2=rep(1,data[3]))
  s4=data.frame(binary_e1=rep(0,data[4]),binary_e2=rep(0,data[4]))
  out=rbind(s1,s2,s3,s4)
  out$binary_ce=ifelse(out$binary_e1+out$binary_e2>0,1,0)
  
  return(out)
}

