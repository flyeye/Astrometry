#-----------------------------------------------------------------
# l, b - degrees, r - parsec
# mu_l, mu_b - "/year
# k = 4.74 "km/sec*parsec"
GetOM_R <- function (l_d,b_d,r)
{
  l <- NISTdegTOradian(l_d)
  b <- NISTdegTOradian(b_d)
  
  result <- vector("numeric", 12)
  
  result[1] <- -cos(l)*cos(b)/r;
  result[2] <- -sin(l)*cos(b)/r;
  result[3] <- -sin(b)/r;
  result[4] <- 0;
  result[5] <- 0;
  result[6] <- 0;
  result[7] <- sin(2*b)*cos(l);                    #M13
  result[8] <- sin(2*b)*sin(l);                    #M23
  result[9] <- cos(b)*cos(b)*sin(2*l);            #M12
  result[10] <- 0.5*cos(b)*cos(b)*cos(2*l);       #M11* = M11-M22
  result[11] <- sin(b)*sin(b)-(1/3);               #x = M33-0.5*(M11+M22)
  result[12] <- 1;                                   #y = 1/3(M11+M22+M33)
  
  result;
}
#-----------------------------------------------------------------
GetOM_L <- function(l_d,b_d,r)
{
  l <- NISTdegTOradian(l_d)
  b <- NISTdegTOradian(b_d)
  
  result <- vector("numeric", 12)
 
  result[1] <- sin(l)/r;                # U
  result[2] <- -cos(l)/r;               # V
  result[3] <- 0/r;                     # W
  result[4] <- -sin(b)*cos(l);          # W1
  result[5] <- -sin(b)*sin(l);          # W2
  result[6] <- cos(b);                  # W3
  result[7] <- -sin(b)*sin(l);          # M13
  result[8] <- sin(b)*cos(l);           # M23
  result[9] <- cos(b)*cos(2*l);         # M12
  result[10] <- -0.5*cos(b)*sin(2*l);   # M11* = M11-M22
  result[11] <- 0;                        # x = M33-0.5*(M11+M22)
  result[12] <- 0;                        # y = 1/3(M11+M22+M33)
  
  result;
}
#-----------------------------------------------------------------
GetOM_B <- function(l_d,b_d,r)
{
  l <- NISTdegTOradian(l_d)
  b <- NISTdegTOradian(b_d)
  
  result <- vector("numeric", 12)
  
  result[1] <- cos(l)*sin(b)/r;
  result[2] <- sin(l)*sin(b)/r;
  result[3] <- -cos(b)/r;
  result[4] <- sin(l);
  result[5] <- -cos(l);
  result[6] <- 0;
  result[7] <- cos(2*b)*cos(l);             #M13
  result[8] <- cos(2*b)*sin(l);             #M23
  result[9] <- -0.5*sin(2*b)*sin(2*l);      #M12
  result[10] <- -0.25*sin(2*b)*cos(2*l);    #M11
  result[11] <- 0.5*sin(2*b);               #x
  result[12] <- 0;                          #y
  
  result;
}

#-----------------------------------------------------------------
# stars - matrix(n,3), where 
# stars[,1] - l in degrees, [,2] - b in degrees, [,3] px - kPc
MakeOMCoef <- function(stars, use_vr = TRUE)
{
  n <- nrow(stars)
  
  if (use_vr == TRUE)
    a0 <- matrix(0, n*3, 12)
  else 
    a0 <- matrix(0, n*2, 12)
  
  for (i in 1:n)
  {
    a0[i,] <- GetOM_L(stars[i,1], stars[i,2], stars[i,3])
    a0[n+i,] <- GetOM_B(stars[i,1], stars[i,2], stars[i,3])
    if (use_vr == TRUE)
      a0[2*n+i,] <- GetOM_R(stars[i,1], stars[i,2], stars[i,3])
  }
  
  if (use_vr == FALSE)
    a0 <- a0[,-12]
  
  return(a0);
}

Calc_OM_Model <- function(stars, use_vr = TRUE, mode = 1, scaling = 0, ef = 0)
{
  #  calculate equation of conditions
  # l, b, px, mu_l, mu_b, vr
  
  a <- MakeOMCoef(stars, use_vr)
  
  if (use_vr == TRUE)
    b <- matrix(rbind(stars[,4]*4.74, stars[,5]*4.74, stars[,6]), nrow(stars)*3, 1)  #*cos(stars[,2])
  else 
    b <- matrix(rbind(stars[,4]*4.74, stars[,5]*4.74), nrow(stars)*2, 1)  #*cos(stars[,2])
  
  res <- TLS_Gen(a, b, mode, scaling, ef)
  
  return(res)
}

#=====================================================================
#----------------------    Test functions    -------------------------
#---------------------------------------------------------------------

GetOM_Default <- function ()
{
  result <- c(10.3, 15.2, 8.0, -2, 1, -15, -1, 15, 0.5, -0.5, 0.5,-0.5);
  return(result)
}

#--------------------------------

MakeTestStars <- function (n)
{
  stars <- matrix(0, n, 3)
  stars[,1]  <- runif(n, min = 0, max = 360)  # l
  stars[,2]  <- runif(n, min = -90, max = 90)  # b
  stars[,3]  <- runif(n, min = 0.1, max = 3)  # px
  
  return(stars)
}

#--------------------------------

Make_OM_Test <- function()
{
  n <- 1000
  
  #  create stars
  stars <- MakeTestStars(n)
  
  #  calculate A0
  a0 <- MakeOMCoef(stars)
  
  #  calculate B0
  OM_0 <- GetOM_Default()
  b0 <- rowSums(t(t(a0)*OM_0))
  
  #  add noise to positions and parallax
  #     l, b, px, vr, ml, mb
  s0 <- c(0.001, 0.001, 0.01, 1, 1, 1)
  #s0 <- c(0, 0, 0, 0, 0, 0)
  noise <- matrix(0, n, 3)
  noise[,1] <- rnorm(n, 0, s0[1])
  noise[,2] <- rnorm(n, 0, s0[2])
  noise[,3] <- rnorm(n, 0, s0[3])
  
  #  calculate A
  a <- MakeOMCoef(stars+noise)
  
  #  calculate B
  b<- vector("numeric", 3*n)
  b[1:n] <- b0[1:n]  + rnorm(n, 0, s0[4])
  b[(n+1):(2*n)] <- b0[(n+1):(2*n)] + rnorm(n, 0, s0[5])
  b[(2*n+1):(3*n)] <- b0[(2*n+1):(3*n)] + rnorm(n, 0, s0[6])
  
  #--------
  

  data <- list(A = a, B = b, A0 = a0, B0 = b0, S_0 = s0, X_0 = OM_0)
  
  m<-12
  ef <- 0
  scaling <- 0
  TLS_Gen_Solve(data, n*3, m, ef, scaling)
  
}
