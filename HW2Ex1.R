
library(tidyverse)

{
  
  #Convenience function for transforming marginals (used to get marginals right after copula has been sampled from))
  TransformMarginals <- function(U, df = 3){
    return(map_dbl(.x = U, .f = qt, df = df))
  }

  #Empirical expected shortfall function
  EmpEs <- function(losses, alpha){
    SortedLosses <- sort(losses, decreasing = TRUE)
    cutoff <- floor(length(losses) * (1-alpha) + 1)
    return(sum(SortedLosses[1:cutoff]) / (cutoff))
  }
  
  #Convenience function for calculating losses
  CalculateLoss <- function(X, S = NULL, DefaultStock = 100){
    if(is.null(S)){
      S <- rep(DefaultStock, length(X))
    }
    else if (length(S) != length(X)){
      S <- rep(S[1], length(X))
    }
    return(-sum(S * (exp(X) - 1)))
  }
  
  #Plotting function to check that the copulas behave as expected
  PairwiseDependencePlot <- function(X, coord1 = x1, coord2 = x2){
    X %>%
      as.data.frame() %>% t() %>% as.data.frame() %>%
      magrittr::set_colnames(seq(1,ncol(.)) %>% map_chr(.f = function(x){paste0("x", x)})) %>%
      ggplot(aes(x = {{ coord1 }}, y = {{ coord2 }})) +
      geom_point()
  }
  
} #General helper functions

{
  #Find the suitable linear correlation for having Kendall's tau = 0.4
  #Using (Theorem 7.42, MFE)
  rhoTarget <- (0.4 * pi / 2) %>% sin()
  
  #Function for constructing the covariance matrix
  MakeCovarianceMatrix <- function(covar = rhoTarget, dim = 50){
    Sigma <- matrix(rep(covar, dim^2), nrow = dim)
    for (i in 1:dim){
      Sigma[i,i] <- 1
    }
    return(Sigma)
  }
  
  #Convenience function for applying the standard normal DF to each coordinate of a vector
  NormalizeMarginals <- function(Z){
    return(map_dbl(.x = Z, .f = pnorm))
  }
  
  #Function for simulating from a Gaussian copula with t-marginals
  SimulateGaussianCopula <- function(nsim = 10, dim = 50, covar = rhoTarget){
    Z <- MASS::mvrnorm(nsim, rep(0,dim), MakeCovarianceMatrix(covar = covar, dim = dim))
    Y <- map(.x = Z %>% t() %>% as.data.frame() %>% as.list(), .f = NormalizeMarginals)
    return(Y)
  }
  
  {
    if(FALSE){
      #Plot the results along two coordinates axes - One cant see tail dependence, but there are some large events due to the usage of the t-distribution
      sample <- MetaGaussianSimulation(nsim = 10000, tdf = 3) %>%
        map(.f = `[`, c(1,2)) %>%
        as.data.frame() %>%
        t() %>%
        as.data.frame()
      colnames(sample) <- c("x","y")
      
      ggplot(data = sample, aes(x = x, y = y)) +
        geom_point()
    }
    
  } #Visualization Example
  
  
} #Helper functions for Gaussian copula simulation

{
  #The theta ensuring the kendall tau is 0.4
  #Using (Table 7.5, MFE)
  ThetaTarget <- 2 * 0.4/(1-0.4)
  
  #Parametrization is taken from p.260 MFE
  PhiClayton <- function(t, theta = ThetaTarget){
    return(
      (1+theta * t)^(-1/theta)
      )
  }
  
  ggplot(mapping = aes(x = seq(0.01, 10, 0.1), y = seq(0.001, 10, 0.1) %>% map_dbl(.f = PhiClayton))) +
    geom_point()
  
  #Function for simulating from the Clayton copula (see Algorithm 7.52, MFE)
  SimulateClayton <- function(nsim = 10, dim = 50, theta = ThetaTarget){
    V <- rgamma(nsim, 1/theta, 1) %>% as.list()
    Y <- 
      rgamma(nsim * dim, 1, 1) %>%
      split(cut(seq_along(.), nsim, labels = F))
    
    U <- map2(.x = Y, .y = V, .f = function(Y,v){
      (Y/v) %>%
        map_dbl(.f = PhiClayton, theta = theta) %>% 
        return()
      })
    
    return(U)
  }
  
} #Helper functions for Clayton copula simulation

{
  set.seed(28122022)
  GaussCopula <- SimulateGaussianCopula(nsim = 10^4) %>%
    map(.f = TransformMarginals) %>% 
    map(.f = function(X){return(X / 100 / sqrt(3))})
  set.seed(NULL)
  GaussLosses <- map_dbl(.x = GaussCopula, .f = CalculateLoss)
  
  #Just for sanity checking
  PairwiseDependencePlot(GaussCopula) #Low tail dependence
  hist(GaussLosses[GaussLosses > 0], breaks = 50)
  
  quantile(GaussLosses, 0.99, type = 1) #VaR
  EmpEs(GaussLosses, 0.99) #ES
  
  set.seed(28122022)
  ClaytonCopula <- SimulateClayton(nsim = 10^4) %>% 
    map(.f = TransformMarginals) %>%
    map(.f = function(X){return(X / 100 / sqrt(3))})
  set.seed(NULL)
  ClaytonLosses <- map_dbl(.x = ClaytonCopula, .f = CalculateLoss)
  
  #Just for sanity checking
  PairwiseDependencePlot(ClaytonCopula) #More tail dependence
  hist(ClaytonLosses[ClaytonLosses > 0], breaks = 50)
  
  quantile(ClaytonLosses, 0.99, type = 1) #VaR
  EmpEs(ClaytonLosses, 0.99) #ES
} #Solution to the exercise




