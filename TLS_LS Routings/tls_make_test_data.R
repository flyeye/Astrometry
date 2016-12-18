# m - parameters qty 
# n - equation qty
# s0 - vector of error in data and variables
TLS_make_test_data <- function (n, m, X, s0)
{

  #X <- rep(0,10)
  equation <- matrix(0,n, m+1)
  equation_c <- matrix(0,n, m+1)
  
  for (i in 1:n)  
  {
     g <- (i-n/2)/(n/50) + rnorm(1);
     #g <- 0;
  
     if (m>=10) 
       equation[i, 10] <- 1;
     if (m>=9) 
       equation[i, 9] <- sin(g*5)^4 
     if (m>=8)
       equation[i, 8] <- sin(g*25)^2 
     if (m>=7) 
       equation[i, 7] <- cos(g*100)^3
     if (m>=6) 
       equation[i, 6] <- 3*cos(g/10+pi/2)^5; 
     if (m>=5) 
       equation[i, 5] <- sin(g*3-3*pi/2)^6;
     if (m>=4) 
       #equation[i, 4] <- 3*log(abs(g+20)); 
       equation[i, 4] <- sin(7*g);  
     if (m>=3)
       #equation[i, 3] <- -1*(g+1)^2; 
       equation[i, 3] <- cos(6*g-3*pi/4);  
     if (m>=2)
       #equation[i, 2] <- 0.2*g;
       equation[i, 2] <- sin(g);  
     equation[i, 1] <- cos(g);  
     
     equation[i, m + 1] <- sum(equation[i,1:m]*X[1:m]);

     equation_c[i,] <- equation[i,]
     
     
     equation[i,1:m] <- equation[i,1:m] + rnorm(m, 0, s0[1:m])
     equation[i, m + 1] <- sum(equation[i,1:m]*X[1:m]) + rnorm(1, 0, s0[m+1])
     
  }
  data <- list(A = equation[,1:m], B = equation[,m+1], A0 = equation_c[,1:m], B0 = equation_c[,m+1], S_0 = s0, X_0 = X)
  data
}

MakeTestData_OLS_test_lm <- function(m, s0)
{
  X_0 <- c(-2, -1, 2, 1, 2, -1, 2, -6, +5, -3);
  print(X_0)
  data<- TLS_make_test_data(100, m, X_0, c(rep(0, m), s0))
  f <- lm(data$B ~ 0 + data$A[,1] + data$A[,2] + data$A[,3] + data$A[,4] + data$A[,5] + data$A[,6] + data$A[,7] + data$A[,8] + data$A[,9] + data$A[,10])
  print(f$coefficients)
  print(sd(f$residuals))
}


MakeTestData_OLS_test_lsfit <- function(m, s0)
{
  X_0 <- c(-2, -1, 2, 1, 2, -1, 2, -6, +5, -3);
  print(X_0)
  data <- TLS_make_test_data(100, m, X_0, c(rep(0, m), s0))
  f <- lsfit(data$A[,1:m], data$B, intercept = FALSE)
  print(f$coefficients)
  print(sd(f$residuals))
}