# ===============================================================================
# --------------------    Astrometric routings    -------------------------------
# -------------------------------------------------------------------------------
# Переход от экваториальных координат, выраженных в радианах,
# в новые галактические, выраженных в радианах. 
# e - matrix, e[,1] - RA, e[,2] - DE, radians
# -------------------------------------------------------------------------------
get_galaxy <- function(e)
{
  
  Leo <- 4.936829261   # 282.85948083°
  L0  <- 0.57477039907 # 32.931918056° 
  si  <- 0.88998807641 # sin 62.871748611° 
  ci  <- 0.45598379779 # cos 62.871748611° 
  
  
  aa <- e[,1]-Leo;
  sa <- sin(aa);   
  ca <- cos(aa);
  
  sd<-sin(e[,2])  
  cd<-cos(e[,2])
  
  g <- cbind( (atan2(sd*si+cd*ci*sa,cd*ca)+L0), asin(sd*ci-cd*si*sa))
  
  g[(g[,1]<0),1] <- g[(g[,1]<0),1]+2*pi;
  
  return(g)
}

# -------------------------------------------------------------------------------
# Перевод собственных движений из экваториальных в галактические
#    pm.a = mu*cos(d)       GalPm.l = mu(l)*cos(b)  
#    pm.d = mu'             GalPm.b = mu(b)              
#    
#  pm - собственные движения в экваториальной СК
# -------------------------------------------------------------------------------
get_galaxy_mu <- function(pm, gal, dec)
{
  
  L0  <- 0.57477039907; # 32.931918056° 
  si  <- 0.88998807641; # sin 62.871748611° 
  ci  <- 0.45598379779; # cos 62.871748611° 
  
  cd <- cos(dec);
  
  sfi <- si*cos(gal[,1]-L0)/cd;
  
  cfi <- ( cos(gal[,2])*ci - (sin(gal[,2])*si)*sin(gal[,1] - L0) )/cd;
  
  GalPm <- cbind( ( cfi*pm[,1]+sfi*pm[,2]), 
                  (-sfi*pm[,1]+cfi*pm[,2]))
  
  return(GalPm)
}

get_galaxy_mu2 <- function(pm, gal, ecv)
{
  
  
  A <- pi*192.85948/180
  D <- pi*27.12825/180
  
  
  sfi <- cos(D) * sin(ecv[,1] - A) / cos(gal[,2]) 
  cfi <- ( sin(D) - sin(gal[,2])*sin(ecv[,2]) ) / (cos(gal[,2]) * cos(ecv[,2]))
  
  
  GalPm <- cbind( ( cfi*pm[,1]+sfi*pm[,2]), 
                  (-sfi*pm[,1]+cfi*pm[,2]))
  
  return(GalPm)
}

# -------------------------------------------------------------------------------
cat_eq2gal<-function(cat_data)
{
  ecv <- cbind(cat_data$RA, cat_data$DE);
  gal <- get_galaxy(ecv)
  
  cat_data <- cat_data %>% mutate(gl = gal[,1], gb = gal[,2])
  
  gal_mu <- get_galaxy_mu(cbind(cat_data$pmRA, cat_data$pmDE), gal, cat_data$DE)
  #gal_mu <- get_galaxy_mu2(cbind(cat_data$pmRA, cat_data$pmDE), gal, ecv)
  cat_data <- cat_data %>% mutate(pm_l = gal_mu[,1], pm_b = gal_mu[,2])
  
  
  #g <-glactc_pm(cat_data$RA*180/pi, cat_data$DE*180/pi, cat_data$pmRA, cat_data$pmDE, year = 2000, j = 1, degree = TRUE, mustar = TRUE)
  #print(mean(g$mu_gl-gal_mu[,1]))
  #cat_data <- cat_data %>% mutate(pm_l = g$mu_gl, pm_b = g$mu_gb)
  
  
  return(cat_data)
}


# ===============================================================================
# -------------------------------------------------------------------------------
#                 Вычисление коэффициентов уравнений модели ОМ
# -------------------------------------------------------------------------------
# l, b - degrees, r - parsec
# mu_l, mu_b - "/year
# k = 4.74 "km/sec*parsec"
GetOM_R <- function ( l, b, r)
{
  #l <- l_d*pi/180
  #b <- b_d*pi/180
  
  result <- vector("numeric", 12)
  
  result[1] <- -cos(l)*cos(b)*4.74064/r;
  result[2] <- -sin(l)*cos(b)*4.74064/r;
  result[3] <- -sin(b)*4.74064/r;
  result[4] <- 0;
  result[5] <- 0;
  result[6] <- 0;
  result[7] <- sin(2*b)*cos(l);                    #M13
  result[8] <- sin(2*b)*sin(l);                    #M23
  result[9] <- cos(b)*cos(b)*sin(2*l);             #M12 = A
  result[10] <- cos(b)*cos(b)*cos(l)*cos(l);       #M11* = M11-M22
  result[11] <- sin(b)*sin(b);                     #M33* = M33-M22
  result[12] <- 1;                                 #M22
  
  return(result);
}
#-----------------------------------------------------------------
GetOM_L <- function( l, b, r)
{
  #l <- l_d*pi/180
  #b <- b_d*pi/180
  
  result <- vector("numeric", 12)
  
  result[1] <- sin(l)*4.74064/r;        # U
  result[2] <- -cos(l)*4.74064/r;       # V
  result[3] <- 0;                       # W
  result[4] <- -sin(b)*cos(l);          # W1
  result[5] <- -sin(b)*sin(l);          # W2
  result[6] <- cos(b);                  # W3 = B
  result[7] <- -sin(b)*sin(l);          # M13
  result[8] <- sin(b)*cos(l);           # M23
  result[9] <- cos(b)*cos(2*l);         # M12 = A
  result[10] <- -0.5*cos(b)*sin(2*l);   # M11* = M11-M22
  result[11] <- 0;                      # M33* = M33-M22
  result[12] <- 0;                      # M22
  
  return(result);
}
#-----------------------------------------------------------------
GetOM_B <- function( l, b, r)
{
  #l <- l_d*pi/180
  #b <- b_d*pi/180
  
  result <- vector("numeric", 12)
   
  result[1] <- cos(l)*sin(b)*4.74064/r;             # X
  result[2] <- sin(l)*sin(b)*4.74064/r;             # Y
  result[3] <- -cos(b)*4.74064/r;                   # Z
  result[4] <- sin(l);                      # W1
  result[5] <- -cos(l);                     # W2 
  result[6] <- 0;                           # W3 = B
  result[7] <- cos(2*b)*cos(l);             # M13
  result[8] <- cos(2*b)*sin(l);             # M23
  result[9] <- -0.5*sin(2*b)*sin(2*l);      # M12 = A
  result[10] <- -0.5*sin(2*b)*cos(l)*cos(l);# M11* = M11-M22
  result[11] <- 0.5*sin(2*b);               # M33* = M33-M22
  result[12] <- 0;                          # M22
  
  return(result);
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
    a0[(n+i),] <- GetOM_B(stars[i,1], stars[i,2], stars[i,3])
    if (use_vr == TRUE)
      a0[(2*n+i),] <- GetOM_R(stars[i,1], stars[i,2], stars[i,3])
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
  
  n <- nrow(stars)
  
  if (use_vr == TRUE){
    b <- matrix(0, n*3, 1)  
    b[1:n,1] <- stars[,4]*4.74064;
    b[(n+1):(2*n),1] <- stars[,5]*4.74064;
    b[(2*n+1):(3*n)] <- stars[,6];
  }
  else 
    {
    b <- matrix(0, n*2, 1)  
    b[1:n,1] <- stars[,4]*4.74064;
    b[(n+1):(2*n),1] <- stars[,5]*4.74064;
  }
  
  #b <- rowSums(t(t(a)*GetOM_Default()))
  
  res <- TLS_Gen(a, b, mode, scaling, ef);
  
  if (use_vr == TRUE)
  {
    names(res$X) <- c("U", "V", "W","W1", "W2", "W3(B)", "M13", "M23", "M12(A)", "M11*", "M33*", "M22")
  }
  else 
  {
    names(res$X) <- c("U", "V", "W","W1", "W2", "W3(B)", "M13", "M23", "M12(A)", "M11*", "M33*")
  }
  
  return(res)
}

#=====================================================================
#----------------------    Test functions    -------------------------
#---------------------------------------------------------------------

GetOM_Default <- function ()
{
  result <- c(10.3, 15.2, 8.0, -2, 1, -15,  -1, 15, 0.5, -0.5, 0.5,-0.5);
  return(result)
}

#--------------------------------

MakeTestStars <- function (n)
{
  stars <- matrix(0, n, 3)
  stars[,1]  <- runif(n, min = 0, max = 360) * pi /180 # l
  stars[,2]  <- runif(n, min = -90, max = 90) * pi /180  # b
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
