
library(tidyverse)

Bmw <- read_tsv("Data/bmw_returns.txt")

Dax <- read_tsv("Data/dax_returns.txt")

StockData <- inner_join(Bmw, Dax, by = "Day", suffix = c("Bmw", "Dax")) %>% 
  rename(LRBmw = LogReturnBmw, LRDax = LogReturnDax)

{
  
  #M-estimation algorithm based on (Algorithm 6.29, MFE)
  M_Estimation <- function(data,
                           maxiter = 10, tol = 10^-6,
                           nu = 3,
                           w_2 = function(x, data){(ncol(data) + nu)/(x + nu)},
                           w_1 = function(x, data){w_2(x^2, data)},
                           ReturnHistory = F){
    
    mu <- vector(mode = "list", maxiter)
    sigma <- vector(mode = "list", maxiter)
    D <- vector(mode = "list", maxiter)
    
    mu[[1]] <- colMeans(data)
    sigma[[1]] <- cov(data)
    D[[1]] <- mahalanobis(data, center = mu[[1]], cov = sigma[[1]])
    
    for (k in 2:maxiter){
      w1 <- map_dbl(.x = D[[k-1]], .f = w_1, data = data)
      w2 <- map_dbl(.x = D[[k-1]], .f = w_2, data = data)
      
      mu[[k]] <- colSums(data * matrix(rep(w1, ncol(data)), ncol = ncol(data))) / 
        colSums(matrix(rep(w1, ncol(data)), ncol = ncol(data)))
      
      sigma[[k]] <-
      (data - (mu[[k-1]] %>% rep(nrow(data)) %>% matrix(ncol = ncol(data), byrow = T))) %>%
        t() %>% as_tibble() %>% as.list() %>% 
        map2(.y = w2, .f = function(x,y){y * (x %*% t(x)) / nrow(data)}) %>% 
        reduce(.f = `+`)
      
      D[[k]] <- mahalanobis(data, center = mu[[k]], cov = sigma[[k]])
      
      if(
        norm(
          (mu[[k]] - mu[[k-1]]) %>% as.matrix()
          , type = "2") < tol 
        || 
        k == maxiter
        ){
        if(ReturnHistory){
          Results <- vector(mode = "list", maxiter) %>% 
            imap(.f = function(x,y){list(mu = mu[[y]], sigma = sigma[[y]])}) %>% 
            `[`(1:k)
        }
        else {
          Results <- vector(mode = "list", maxiter) %>% 
            imap(.f = function(x,y){list(mu = mu[[y]], sigma = sigma[[y]])}) %>% 
            `[[`(k)
        }
        
        return(Results)
      }
      
    }
    
  }
  
  #Test
  #M_est_results <- M_Estimation(StockData %>% select(LRBmw, LRDax), maxiter = 100)
  
  
  
} #Estimating the mean vector and dispersion matrix

{
  #Function for iteratively finding the best df (nu) for the elliptical data fit
  #OBS - The current implementation uses stochastic search and is quite sensity to initial conditions
  ProfileNu <- function(data,
                        InitialNu = 3,
                        StepsizeFunc = function(k){1/log(k+1)},
                        maxiter = 10, tol = 10^-6,
                        innerMaxiter = 10, innerTol = 10^-6){
    
    nu <- rep(NA_real_, maxiter)
    EllipticalParams <- vector(mode = "list", maxiter)
    LikelihoodScores <- rep(NA_real_, maxiter)
    nu[[1]] <- InitialNu
    
    best_index <- 1
    
    for(k in 1:maxiter){
      #Fit the elliptical distribution parameters using M-estimation
      EllipticalParams[[k]] <- M_Estimation(data,
                                            nu = nu[[k]],
                                            maxiter = innerMaxiter, tol = innerTol)
      
      #Calculate t-likelihood for the given data
      LikelihoodScores[[k]] <- 
        data %>% 
        TibbleToVectorlist() %>% 
        TransformToSpherical(mean = EllipticalParams[[k]][["mu"]],
                             dispersion = EllipticalParams[[k]][["sigma"]]) %>% 
        map_dbl(.f = norm, type = "2") %>% 
        tLikelihood(nu = nu[[k]])
      
      #Update (nu) to reflect the observed likelihood
      if(k == 1){
        nu[[k+1]] <- nu[[k]] + 1
      } else if (F) {
        nu[[k+1]] <-
          (LikelihoodScores[[k]] - LikelihoodScores[[k-1]])/(nu[[k]] - nu[[k-1]]) * StepsizeFunc(k) +
          nu[[k]]
      }
      else {
        if(LikelihoodScores[[k]] > LikelihoodScores[[best_index]]){
          best_index <- k
        }
        nu[[k+1]] <- (2 * rbinom(1, 1, 0.5) - 1) * StepsizeFunc(k) + nu[[best_index]]
      }
      
      
    }
    
    
    
    return(list(nu = nu, Elliptical = EllipticalParams, Likelihood = LikelihoodScores))
    
  }
  
  #Convert the data from tibble format to list format
  TibbleToVectorlist <- function(df){
    df %>% as.matrix() %>% t() %>% as_tibble() %>% as.list()
  }
  
  #Calculate the t-distribution loglikelihood
  tLikelihood <- function(tData, nu){
    tData %>% 
      map_dbl(.f = dt, df = nu) %>%
      log() %>% 
      sum()
  }
  
  #Given mean and dispersion, transform the data into spherical form (assumes original data to be elliptical)
  TransformToSpherical <- function(data, mean, dispersion){
    data %>% 
      map(.f = function(x){x - mean}) %>% 
      map(.f = function(x){solve(chol(dispersion)) %*% x})
  }
} #Estimating the best (in a likelihood sense) degree of freedom for a t-distribution

{
  #Fit an elliptical distribution based on a t-distribution, but optimize over the degree's of freedom
  TDistributions <-
  ProfileNu(StockData %>% select(-Day),
            maxiter = 100,
            StepsizeFunc = function(k){20 / k})

  #Plots for checking that the algorithm converged
  TDistributions[["Likelihood"]] %>%
    ggplot(data = NULL, mapping = aes(x = seq_along(.), y = .)) +
    geom_point()

  TDistributions[["nu"]] %>%
    ggplot(data = NULL, mapping = aes(x = seq_along(.), y = .)) +
    geom_point()
  
  #Get the best value of (nu, mu and sigma) - Right now I'm being sloppy and just taking the last one
  nu_optimal <- TDistributions[["nu"]] %>% tail(1)
  mu_optimal <- TDistributions[["Elliptical"]] %>% tail(1) %>% unlist(recursive = F) %>% `[[`("mu")
  sigma_optimal <- TDistributions[["Elliptical"]] %>% tail(1) %>% unlist(recursive = F) %>% `[[`("sigma")
  
  #Simulate a large number of losses from the implied distribution
  Nsim <- 10^4
  radialSim <- runif(Nsim, min = 0, max = 2 * pi)
  magnitudeSim <- rt(Nsim, df = nu_optimal)
  lrSim <- radialSim %>% 
    map2(.y = magnitudeSim, .f = function(t,r){c(cos(t), sin(t)) * r}) %>% 
    map(.f = function(x){chol(sigma_optimal) %*% x + mu_optimal})
  
  M_est_fake_results <- M_Estimation(StockData %>% select(-Day),
               maxiter = 2, nu = 5)
  mu_fake <- M_est_fake_results[["mu"]]
  sigma_fake <- M_est_fake_results[["sigma"]]
  
  radialSimFake <- runif(Nsim, min = 0, max = 2 * pi)
  magnitudeSimFake <- rt(Nsim, df = 5)
  lrSimFake <- radialSimFake %>% 
    map2(.y = magnitudeSimFake, .f = function(t,r){c(cos(t), sin(t)) * r}) %>% 
    map(.f = function(x){chol(sigma_fake) %*% x + mu_fake})
  
  #Generally speaking, the plots show that optimizing wrt. loglikelihood yields much too
  #light tails (all the blue points are in the center of the distribution)
  plotdata <- lrSim %>% reduce(cbind) %>% t() %>% magrittr::set_colnames(c("LRBmw", "LRDax")) %>%
    as_tibble()
  
  plotdata_fake <- lrSimFake %>% reduce(cbind) %>% t() %>% magrittr::set_colnames(c("LRBmw", "LRDax")) %>%
    as_tibble()
  
  ggplot() +
    geom_point(mapping = aes(x = plotdata_fake$LRBmw, y = plotdata_fake$LRDax), color = "red") +
    geom_point(mapping = aes(x = plotdata$LRBmw, y = plotdata$LRDax), color = "blue") +
    geom_point(mapping = aes(x = StockData$LRBmw, y = StockData$LRDax))
  
} #Solving the exercise





