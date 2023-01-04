
library(tidyverse)
library(MASS)

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
  #The eyeballed df = 5 fit looks better, but the covariance seems to be off
  plotdata <- lrSim %>% reduce(cbind) %>% t() %>% magrittr::set_colnames(c("LRBmw", "LRDax")) %>%
    as_tibble()
  
  plotdata_fake <- lrSimFake %>% reduce(cbind) %>% t() %>% magrittr::set_colnames(c("LRBmw", "LRDax")) %>%
    as_tibble()
  
  ggplot() +
    geom_point(mapping = aes(x = plotdata_fake$LRBmw, y = plotdata_fake$LRDax), color = "red") +
    geom_point(mapping = aes(x = plotdata$LRBmw, y = plotdata$LRDax), color = "blue") +
    geom_point(mapping = aes(x = StockData$LRBmw, y = StockData$LRDax))
  
  
  if(exists("CalculateLoss") == F || exists("EmpEs") == F){
    warning("The function(s) CalculateLoss/EmpEs are not loaded. They can be found in the Ex1 code.")
  } else {
    #Calculate simulated losses
    simLosses <- 
      lrSimFake %>% 
      map_dbl(.f = CalculateLoss)
    
    #Calculate empirical VaR and ES
    quantile(simLosses, 0.99, type = 1) %>% paste("VaR:", .) %>%  print()
    EmpEs(simLosses, 0.99) %>% paste("ES:", .) %>% print()
  }
  
} #Solving the exercise

{
  
  m_hat <- 
    StockData %>% 
    select(-Day) %>% 
    colMeans()
  
  S_hat <- StockData %>% 
    select(-Day) %>% 
    cov()
  
  mu_hat <- t(-c(100, 100)) %*% m_hat
  
  sigma_hat <- t(-c(100, 100)) %*% S_hat %*% (-c(100, 100))
  
  VaR_vc <- mu_hat + sigma_hat * qnorm(0.99)
  
  
}#Variance-Covariance implementation

{
  
  StockData %>% 
    ggplot(aes(x = LRDax)) +
    geom_histogram()
  
  StockData %>% 
    ggplot(aes(x = LRBmw)) +
    geom_histogram()
  
  BmwFit <- 
  StockData[["LRBmw"]] %>% 
    fitdistr(densfun = "t")
  
  DaxFit <- 
  StockData[["LRDax"]] %>% 
    fitdistr(densfun = "t")
  
  xForPlot <- seq(-0.1, 0.1, 0.001)
  yBmw <- xForPlot %>%
    as.list() %>% 
    map_dbl(.f = function(x){
      (x - BmwFit$estimate["m"])/BmwFit$estimate["s"]}
    ) %>%
    map_dbl(.f = dt, df = BmwFit$estimate["df"]) %>% 
    `*`(1700)
  
  yDax <- xForPlot %>%
    as.list() %>% 
    map_dbl(.f = function(x){
      (x - DaxFit$estimate["m"])/DaxFit$estimate["s"]}
    ) %>%
    map_dbl(.f = dt, df = DaxFit$estimate["df"]) %>% 
    `*`(1700)
  
  #Looking pretty good
  ggplot() +
    geom_histogram(mapping = aes(x = StockData$LRBmw)) +
    geom_line(mapping = aes(x = xForPlot, y = yBmw), color = "red")
  
  ggplot() +
    geom_histogram(mapping = aes(x = StockData$LRDax)) +
    geom_line(mapping = aes(x = xForPlot, y = yDax), color = "red")
    
  StockData <- 
  StockData %>% 
    mutate(UDax = pt((LRDax - DaxFit$estimate["m"])/DaxFit$estimate["s"], df = DaxFit$estimate["df"]),
           UBmw = pt((LRBmw - BmwFit$estimate["m"])/BmwFit$estimate["s"], df = BmwFit$estimate["df"]))
  
  ggplot(StockData, aes(x = UDax, y = UBmw)) +
    geom_point()
  
  
  #Likelihood for a Generalized Clayton copula
  #theta >= 0, delta >= 1
  
  phi_gc <- function(t,delta,theta){
    (1 + theta * t^(1/delta))
  }
  
  phi_gc_inv <- function(u,delta,theta){
    (u^(-theta) / theta - 1)^(delta)
  }
  
  phi_gc_prime <- function(t, delta, theta){
    -1/delta * (1 + theta * t^(1/delta)) * t^(1/delta - 1)
  }
  
  phi_gc_primeprime <- function(t,delta,theta){
    -delta * (1 + theta * t^(1/delta))^(-1/theta - 2)*(-1/theta - 1) * t^(1/delta - 1) * theta * (1/delta) * t^(1/delta - 1) +
      -1/delta * (1 + theta * t^(1/delta))^(-1/theta - 1) * (1/delta - 1)*t^(1/delta - 2)
  }
  
  likelihood <- function(u_1, u_2, delta, theta){
    return(
    phi_gc_primeprime(
      phi_gc_inv(u_1, delta, theta) + phi_gc_inv(u_2, delta, theta)
      , delta, theta) / 
      phi_gc_prime(phi_gc_inv(u_1, delta, theta), delta, theta) /
      phi_gc_prime(phi_gc_inv(u_2, delta, theta), delta, theta)
    )
  }
  
  TileSearch <- function(lowerLimits, upperLimits, objective, subdivisions = 1, iterations = 4){
    if(iterations == 4){
      #browser()
    }
    d <- length(lowerLimits)
    repLengthForThisRun <- subdivisions + 1
    alpha_midpoint <- seq(1/(subdivisions + 1), subdivisions / (subdivisions + 1), 1/(subdivisions + 1))
    
    evals <- 
    1:d %>% 
      map(.f = function(i){
        map2(.x=c(0, alpha_midpoint), .y = c(alpha_midpoint, 1), .f = function(x,y){(x+y)/2}) %>% 
          map_dbl(.f = function(alpha){
            lowerLimits[i] * (1-alpha) + upperLimits[i] * (alpha)
          })
      })
    
    #Noget med at loope over alle punkter og tage det bedste og køre optimizeren rekursivt herfra
    
    results <- 
      (0:((subdivisions + 1)^d-1)) %>% 
      map(.f = BaseD, repLength = repLengthForThisRun) %>% 
      map_dbl(.f = function(index){
        map2_dbl(.x = evals, .y = index, .f = .subset2) %>% objective()
      })
    
    bestResultIndex <- which.max(results)
    
    diag <- (upperLimits - lowerLimits) / (subdivisions + 1)
    
    newLowerBound <- map2_dbl(.x = evals, .y = BaseD(bestResultIndex - 1, d, repLengthForThisRun), .f = .subset2) - diag/2
    newUpperBound <- map2_dbl(.x = evals, .y = BaseD(bestResultIndex - 1, d, repLengthForThisRun), .f = .subset2) + diag/2
    
    if(iterations > 0){
      return(TileSearch(newLowerBound, newUpperBound, objective, subdivisions, iterations - 1))
    }
    else {
      return((newLowerBound + newUpperBound) / 2)
    }
    
  }
  
  testFunc <- function(x){
    -(x[1] - 2)^2 - (x[2] - 3)^2
  }
  
  TileSearch(c(2,2), c(3,3), testFunc, iterations = 20)
  
  objectiveFunction <- function(params){
    Observations %>% 
    map_dbl(.f = function(obs){
      likelihood(obs[1], obs[2], params[1], params[2])
    }) %>% 
      prod()
  }
  
  likelihood(5, 2, 2, 5)
  
  objectiveFunction(c(4,5))
  
  Observations %>% 
    map_dbl(.f = function(obs){
      likelihood(obs[1], obs[2], 3, 3)
    }) %>% 
    prod()
  
  Observations <- 
  StockData %>% 
    dplyr::select(UDax, UBmw) %>% 
    as.list() %>% 
    transpose() %>% 
    map(.f = unlist)
  
  
  
  #Produces the base-d representation of a given integer
  BaseD <- function(x, d = 2, repLength = NULL, rIndex = T){
    
    #Decide the length of the base-d representation
    if(is.null(repLength)){
      finalPower <- log(x, base = d) %>% ceiling()
    }
    else {
      finalPower <- repLength
    }
    
    #Calculate the base-d representation
    baseRep <- map_dbl(.x = 1:finalPower, .f = function(k){
      (x %% (d^k) - x %% (d^(k-1))) / d^(k-1)
      })
    if(rIndex == F){
      return(baseRep)
    }
    else {
      return(baseRep + 1)
    }
    
  }
  
} #Copula approach


