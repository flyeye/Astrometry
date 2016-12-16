# m - parameters qty 
# n - equation qty
# s0 - vector of error in data and variables
TLS_make_test_data <- function (n, m, X, s0)
{

  #X <- rep(0,10)
  equation <- matrix(0,n, m+1)
  
  for (i in 1:n)  
  {
     g <- rnorm(1)*10;
     #g <- 0;
  
     if (m>=10) 
       equation[i, 10] <- 1;
     if (m>=9) 
       equation[i, 9] <- sin(g*5/pi)^4 
     if (m>=8)
       equation[i, 8] <- sin(g*25/pi)^2 
     if (m>=7) 
       equation[i, 7] <- cos(g*100/pi)^3
     if (m>=6) 
       equation[i, 6] <- g^5; 
     if (m>=5) 
       equation[i, 5] <- g^4; 
     if (m>=4) 
       equation[i, 4] <- g^3; 
     if (m>=3)
       equation[i, 3] <- g^2; 
     if (m>=2)
       equation[i, 2] <- g;
     equation[i, 1] <- cos(g/pi);  
     
     equation[i, m + 1] <- sum(equation[i,1:m]*X[1:m]) + rnorm(1, 0, s0[m+1])
     
     equation[i,1:m] <- equation[i,1:m] + rnorm(m, 0, s0)
     
     #message(equation[i,])
  }
  
  equation
}

MakeTestData_OLS_test_lm <- function(m, s0)
{
  X_0 <- c(-2, -1, 2, 1, 2, -1, 2, -6, +5, -3);
  print(X_0)
  e <- TLS_make_test_data(100, m, X_0, c(rep(0, m), s0))
  f <- lm(e[,m+1] ~ 0 +e[,1] + e[,2] + e[,3] + e[,4] + e[,5] + e[,6] + e[,7] + e[,8] + e[,9] + e[,10])
  print(f$coefficients)
  print(sd(f$residuals))
}


MakeTestData_OLS_test_lsfit <- function(m, s0)
{
  X_0 <- c(-2, -1, 2, 1, 2, -1, 2, -6, +5, -3);
  print(X_0)
  e <- TLS_make_test_data(100, m, X_0, c(rep(0, m), s0))
  f <- lsfit(e[,1:m], e[,m+1], intercept = FALSE)
  print(f$coefficients)
  print(sd(f$residuals))
}