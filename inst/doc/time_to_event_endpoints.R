## ----setup, include=FALSE, echo=FALSE-----------------------------------------
library('CompAREdesign')            # Load library
library('ggplot2')
knitr::opts_chunk$set(echo = TRUE)  # Show code

## -----------------------------------------------------------------------------
## Probabilities of observing the event in control arm during follow-up
p0_e1 <- 0.59   # Death
p0_e2 <- 0.74   # Disease Progression  

## Effect size (Cause specific hazard ratios) for each endpoint
HR_e1 <- 0.91   # Death
HR_e2 <- 0.77   # Disease Progression

## Hazard rates over time
beta_e1 <- 2    # Death --> Increasing risk over time
beta_e2 <- 1    # Disease Progression --> Constant risk over time

## Correlation
rho      <- 0.1        # Correlation between components
rho_type <- 'Spearman' # Type of correlation measure
copula   <- 'Frank'    # Copula used to get the joint distribution

## Additional parameter
case  <- 3  # 1: No deaths;                  2: Death is the secondary event; 
            # 3: Death is the primary event; 4: Both events are death by different causes

## -----------------------------------------------------------------------------
ARE_tte(p0_e1   = p0_e1    , p0_e2    = p0_e2, 
        HR_e1   = HR_e1    , HR_e2    = HR_e2, 
        beta_e1 = beta_e1  , beta_e2  = beta_e2,  
        rho     = rho      , rho_type = rho_type, 
        copula  = copula   , case     = case,
        plot_print = TRUE) 

## ----warning=FALSE------------------------------------------------------------
effectsize_tte(p0_e1   = p0_e1    , p0_e2    = p0_e2, 
               HR_e1   = HR_e1    , HR_e2    = HR_e2, 
               beta_e1 = beta_e1  , beta_e2  = beta_e2,  
               rho     = rho      , rho_type = rho_type, 
               copula  = copula   , case     = case,
               plot_print = TRUE) 

## -----------------------------------------------------------------------------
samplesize_tte(p0_e1   = p0_e1    , p0_e2    = p0_e2, 
               HR_e1   = HR_e1    , HR_e2    = HR_e2, 
               beta_e1 = beta_e1  , beta_e2  = beta_e2,  
               rho     = rho      , rho_type = rho_type, 
               copula  = copula   , case     = case, 
               alpha   = 0.025    , power    = 0.90,
               ss_formula = 'schoenfeld', 
               plot_print = TRUE) 

## ----warning=FALSE------------------------------------------------------------
## Hazards' rates over time Scenario 1
beta_e1 <- 1    # Death               --> CONSTANT over time
beta_e2 <- 2    # Disease Progression --> INCREASE over time
effectsize_tte(p0_e1   = p0_e1    , p0_e2    = p0_e2, 
               HR_e1   = HR_e1    , HR_e2    = HR_e2, 
               beta_e1 = beta_e1  , beta_e2  = beta_e2,  
               rho     = rho      , rho_type = rho_type, 
               copula  = copula   , case     = case, 
               plot_print = TRUE) 

## Hazards' rates over time Scenario 2
beta_e1 <- 1    # Death               --> CONSTANT over time
beta_e2 <- 1    # Disease Progression --> CONSTANT over time
effectsize_tte(p0_e1   = p0_e1    , p0_e2    = p0_e2, 
               HR_e1   = HR_e1    , HR_e2    = HR_e2, 
               beta_e1 = beta_e1  , beta_e2  = beta_e2,  
               rho     = rho      , rho_type = rho_type, 
               copula  = copula   , case     = case, 
               plot_print = TRUE) 

## ----warning=FALSE, fig.width=6, fig.height=6---------------------------------
plot_tte(p0_e1   = p0_e1    , p0_e2    = p0_e2, 
         HR_e1   = HR_e1    , HR_e2    = HR_e2, 
         beta_e1 = beta_e1  , beta_e2  = beta_e2,  
         rho     = rho      , rho_type = rho_type, 
         copula  = copula   , case     = case, 
         summary = TRUE) 

## -----------------------------------------------------------------------------
plot_tte(p0_e1   = p0_e1    , p0_e2    = p0_e2, 
         HR_e1   = HR_e1    , HR_e2    = HR_e2, 
         beta_e1 = beta_e1  , beta_e2  = beta_e2,  
         rho     = rho      , rho_type = rho_type, 
         copula  = copula   , case     = case, 
         summary = FALSE    , type     = 'ARE') + 
  ggtitle('Asymptotic Relative Efficiency') + theme_bw()

## ----echo = FALSE-------------------------------------------------------------
rm(list = ls())

