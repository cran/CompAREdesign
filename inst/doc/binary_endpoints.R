## ----setup, include=FALSE, echo=FALSE-----------------------------------------
library('CompAREdesign')            # Load library
knitr::opts_chunk$set(echo = TRUE)  # Show code

## -----------------------------------------------------------------------------
## Probabilities of observing the event in control arm at the end of follow-up
p0_e1 <- 0.059    # Ischemia-driven target-lesion revascularization 
p0_e2 <- 0.032    # Cardiac death or target-vessel MI 

## Effect size (absolute reduction) for each endpoint
AR_e1 <- -0.0196  # Ischemia-driven target-lesion revascularization 
AR_e2 <- -0.0098  # Cardiac death or target-vessel MI

## Correlation
rho   <- 0.4

## -----------------------------------------------------------------------------
ARE_cbe(p0_e1   = p0_e1  , p0_e2   = p0_e2, 
        eff_e1  = AR_e1  , eff_e2  = AR_e2, 
        effm_e1 = "diff" , effm_e2 = "diff", effm_ce = "or",
        rho     = rho) 

## -----------------------------------------------------------------------------
effectsize_cbe(p0_e1   = p0_e1  , p0_e2   = p0_e2, 
               eff_e1  = AR_e1  , eff_e2  = AR_e2, 
               effm_e1 = "diff" , effm_e2 = "diff", effm_ce = "or",
               rho     = rho) 

## -----------------------------------------------------------------------------
samplesize_cbe(p0_e1   = p0_e1  , p0_e2   = p0_e2, 
               eff_e1  = AR_e1  , eff_e2  = AR_e2, 
               effm_e1 = "diff" , effm_e2 = "diff", effm_ce = "or",
               rho     = rho,
               alpha   = 0.05, beta = 0.2) 

## ----echo = FALSE-------------------------------------------------------------
rm(list = ls())

