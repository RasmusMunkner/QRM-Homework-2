#Disproving that C_tilde is a copula

#The function we wish to prove is not a copula
C_tilde <- function(x,y,z){
  max(x*y + z - 1,0)
}

#Evaluate the corresponding (Stieltjes) measure on [a_1, b_1] x [a_2, b_2] x [a_3, b_3]
CubeEval <- function(a,b,f = C_tilde){
  target <- 0
  
  for(i in 1:8){
    
    s <- 1
    
    if(i %in% c(1,2,3,4)){
      x <- a[1]
      s <- (-1) * s
    } else {
      x <- b[1]
    }
    
    if(i %in% c(2,4,6,8)){
      y <- a[2]
      s <- (-1) * s
    } else {
      y <- b[2]
    }
    
    if(i %in% c(1,2,5,6)){
      z <- a[3]
      s <- (-1) * s
    } else {
      z <- b[3]
    }
    target <- target + s * f(x,y,z)
    
  }
  
  return(target)
  
}

#Check that stuff breaks down for sqrt(1/3)
CubeEval( a = c(sqrt(1/3), sqrt(1/3), sqrt(1/3)), b = c(1,1,1))

#A lille extra for visualization
wrapper <- function(x, f){
  CubeEval(a = rep(x, 3), b = rep(1,3), f = f)
}

library(tidyverse)
w <- seq(0.001, 1, 0.0001)
fw <- map_dbl(.x = w, .f = wrapper, f = C_tilde)
ggplot() +
  geom_point(mapping = aes(x = w, y = fw)) +
  geom_vline(xintercept = (-1 + sqrt(5))/2)


