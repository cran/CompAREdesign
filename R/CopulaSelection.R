#' Copula object
#'
#' @description  Build a copula class object from the given family and the corresponding dependence parameter from the given correlation
#' @param copula  Copula given from the following list: "Frank" (default), "Gumbel" , "Clayton", "Normal" , "T",
#'                , "Galambos", "HuslerReiss", "Tawn", "Tev", "FGM" or "Plackett"
#' @param rho	    Spearman's coefficient between the 2 marginal distributions
#' @param rhoType type of correlation: "Spearman" (default) or "Kendall"
#'
#'
#' @export 
#' @keywords internal 
#'
#'
#'


CopulaSelection <- function(copula,rho,rhoType='Spearman'){
  
  if(rhoType=='Spearman'){
    theta <- switch(copula,
                    Frank =       iRho(frankCopula(1),rho),
                    Gumbel =      iRho(gumbelCopula(2),rho),
                    Clayton =     iRho(claytonCopula(1),rho),
                    FGM =         iRho(fgmCopula(1),rho),
                    Normal =      iRho(normalCopula(0.5),rho),
                    'T' =         iRho(normalCopula(0.5),rho), 
                    Galambos =    iRho(galambosCopula(0.5),rho),
                    HuslerReiss = iRho(huslerReissCopula(0.5),rho),
                    Tawn =        iRho(tawnCopula(0.5),rho),
                    Tev =         iRho(tevCopula(0.5),rho),
                    Plackett =    iRho(plackettCopula(0.5),rho))
  }else{
    theta <- switch(copula,
                    Frank =       iTau(frankCopula(1),rho),
                    Gumbel =      iTau(gumbelCopula(2),rho),
                    Clayton =     iTau(claytonCopula(1),rho),
                    FGM =         iTau(fgmCopula(1),rho),
                    Normal =      iTau(normalCopula(0.5),rho),
                    'T' =         iTau(normalCopula(0.5),rho), 
                    Galambos =    iTau(galambosCopula(0.5),rho),
                    HuslerReiss = iTau(huslerReissCopula(0.5),rho),
                    Tawn =        iTau(tawnCopula(0.5),rho),
                    Tev =         iTau(tevCopula(0.5),rho),
                    Plackett =    iTau(plackettCopula(0.5),rho))
  }
  
  which.copula <- switch(copula,
                         Frank =       archmCopula(family = "frank", dim = 2, param = theta),
                         Gumbel =      archmCopula(family = "gumbel", dim = 2, param = theta),
                         Clayton =     archmCopula(family = "clayton", dim = 2, param = theta),
                         FGM =         fgmCopula(dim = 2, param = theta),
                         Normal =      normalCopula(dim = 2, param = theta),
                         'T' =         tCopula(dim = 2, param = theta),
                         Galambos =    galambosCopula(param = theta),
                         HuslerReiss = huslerReissCopula(param = theta),
                         Tawn =        tawnCopula(param = theta),
                         Tev =         tevCopula(param = theta),
                         Plackett =    plackettCopula(param = theta ))
  
  
  return(c(which.copula,theta))
}