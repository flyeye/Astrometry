# ===============================================================================
# --------------------    Astrometric routings    -------------------------------
# -------------------------------------------------------------------------------
# Переход от экваториальных координат, выраженных в радианах,
# в новые галактические, выраженных в радианах. 
# e - matrix, e[,1] - RA, e[,2] - DE, radians
# -------------------------------------------------------------------------------
get_galaxy <- function(e)
{

  Leo <- 4.936829261   # 282.85948083?
  L0  <- 0.57477039907 # 32.931918056?
  si  <- 0.88998807641 # sin 62.871748611?
  ci  <- 0.45598379779 # cos 62.871748611?


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

  L0  <- 0.57477039907; # 32.931918056?
  si  <- 0.88998807641; # sin 62.871748611?
  ci  <- 0.45598379779; # cos 62.871748611?

  cd <- cos(dec);

  sfi <- si*cos(gal[,1]-L0)/cd;

  cfi <- ( cos(gal[,2])*ci - (sin(gal[,2])*si)*sin(gal[,1] - L0) )/cd;

  GalPm <- cbind( ( cfi*pm[,1]+sfi*pm[,2]),
                  (-sfi*pm[,1]+cfi*pm[,2]))

  return(GalPm)
}

calc_galaxy_mu_DT <- function(data)
{
  
  L0  <- 0.57477039907; # 32.931918056?
  si  <- 0.88998807641; # sin 62.871748611?
  ci  <- 0.45598379779; # cos 62.871748611?
  
  data[, cs := cos(DE)]
  data[, sfi := si*cos(gl-L0)/cs]
  data[,  cfi := (cos(gb)*ci - (sin(gb)*si)*sin(gl - L0) )/cs]
  data[, pm_l := cfi*gpmRA+sfi*gpmDE]
  data[, pm_b := -sfi*gpmRA+cfi*gpmDE]
  data[, cs := NULL]
  data[, sfi := NULL]
  data[, cfi := NULL]
  gc()
  
  return(data)
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
# l, b - degrees, r - kiloparsec
# mu_l, mu_b - "/year
# k = 4.74 "km/sec*kiloparsec"
GetOM_R <- function ( l, b, r, type = 0)
{
  #l <- l_d*pi/180
  #b <- b_d*pi/180

  result <- vector("numeric", 12)

  result[1] <- -cos(l)*cos(b)/r;
  result[2] <- -sin(l)*cos(b)/r;
  result[3] <- -sin(b)/r;
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

OM_U_R <- function ( l, b, r)         # U
{
  return( -cos(l)*cos(b)/r );
}

OM_V_R <- function ( l, b, r)         # V
{
  return( -sin(l)*cos(b)/r );
}

OM_W_R <- function ( l, b, r)         # W
{
  return( -sin(b)/r);
}

OM_w1_R <- function ( l, b, r)         # w1
{
  return( 0 );
}

OM_w2_R <- function ( l, b, r)         # w2
{
  return( 0 );
}

OL_B_R <- function ( l, b, r)               # B
{
  return( OM_w3_R(l, b, r));
}

OM_w3_R <- function ( l, b, r)         # w3
{
  return( 0 );
}

OM_M13_R <- function ( l, b, r)
{
  return( sin(2*b)*cos(l) );                    #M13
}
 
OM_M23_R <- function ( l, b, r)
{
  return( sin(2*b)*sin(l) );                    #M23
}

OL_A_R <- function ( l, b, r)               # A
{
  return( OM_M12_R(l, b, r));
}

OM_M12_R <- function ( l, b, r)
{
  return( cos(b)*cos(b)*sin(2*l) );             #M12 = A
}

OM_M11s_R <- function ( l, b, r)
{
  return( cos(b)*cos(b)*cos(l)*cos(l) );       #M11* = M11-M22
}

OM_M33s_R <- function ( l, b, r)
{
  return( sin(b)*sin(b) );                     #M33* = M33-M22
}

OM_M22_R <- function ( l, b, r)
{
  return( 1 );                                 #M22
}


OL_C_R <- function ( l, b, r)
{
  return( 0 );                                 # C, не реализовано
}

OL_K_R <- function ( l, b, r)
{
  return( 0 );                                 # K, не реализовано
}

#-----------------------------------------------------------------
GetOM_L <- function( l, b, r, type = 0)
{
  #l <- l_d*pi/180
  #b <- b_d*pi/180

  result <- vector("numeric", 12)

  result[1] <- sin(l)/r;        # U
  result[2] <- -cos(l)/r;       # V
  result[3] <- 0;                       # W
  result[4] <- -sin(b)*cos(l);          # W1
  result[5] <- -sin(b)*sin(l);          # W2
  result[6] <- cos(b);                  # W3 = B
  result[7] <- -sin(b)*sin(l);          # M13
  result[8] <- sin(b)*cos(l);           # M23
  result[9] <- cos(b)*cos(2*l);         # M12 = A
  if(type==0)
  {
    result[10] <- -0.5*cos(b)*sin(2*l);   # M11* = M11-M22
  } else {
    result[10] <- -cos(b)*sin(2*l);   # C
  }
  result[11] <- 0;                      # M33* = M33-M22
  result[12] <- 0;                      # M22

  return(result);
}

OM_U_L <- function( l, b, r)
{
  return( sin(l)/r);        # U
}  
  
OM_V_L <- function( l, b, r)
{
  return( -cos(l)/r );       # V
}

OM_W_L <- function( l, b, r)
{
  return( 0 );                       # W
}

OM_w1_L <- function( l, b, r)
{
  return( -sin(b)*cos(l) );          # W1
}

OM_w2_L <- function( l, b, r)
{
  return( -sin(b)*sin(l) );          # W2
}

OL_B_L <- function( l, b, r)
{
  return(OM_w3_L(l, b, r))
}

OM_w3_L <- function( l, b, r)
{
  return( cos(b) );                  # W3 = B
}

OM_M13_L <- function( l, b, r)
{
  return( -sin(b)*sin(l) );          # M13
}

OM_M23_L <- function( l, b, r)
{
  return( sin(b)*cos(l) );           # M23
}

OL_A_L <- function(l, b, r)
{
  return(OM_M12_L(l, b, r))
}

OM_M12_L <- function( l, b, r)
{
  return( cos(b)*cos(2*l) );         # M12 = A
}

OM_M11s_L <- function( l, b, r)
{
  return( -0.5*cos(b)*sin(2*l));   # M11* = M11-M22
}

OL_C_L <- function( l, b, r)
{
  return( -cos(b)*sin(2*l));   # C
}

OL_K_L <- function( l, b, r)
{
  return(0);   # K, не определяется
}
 
OM_M33s_L <- function( l, b, r)
{
  return(0);                      # M33* = M33-M22, не определяется 
}

OM_M22_L <- function( l, b, r)
{
  return(0);                      # M22, не определяется
}

#-----------------------------------------------------------------
GetOM_B <- function( l, b, r, type = 0)
{
  #l <- l_d*pi/180
  #b <- b_d*pi/180

  result <- vector("numeric", 12)

  result[1] <- cos(l)*sin(b)/r;             # X
  result[2] <- sin(l)*sin(b)/r;             # Y
  result[3] <- -cos(b)/r;                   # Z
  result[4] <- sin(l);                      # W1
  result[5] <- -cos(l);                     # W2
  result[6] <- 0;                           # W3 = B
  result[7] <- cos(2*b)*cos(l);             # M13
  result[8] <- cos(2*b)*sin(l);             # M23
  result[9] <- -0.5*sin(2*b)*sin(2*l);      # M12 = A
  if (type == 0)
  {
    result[10] <- -0.5*sin(2*b)*cos(l)*cos(l);# M11* = M11-M22
    result[11] <- 0.5*sin(2*b);               # M33* = M33-M22
  } else 
  {
    result[10] <- -0.5*sin(2*b)*cos(2*l);     # C
    result[11] <- -0.5*sin(2*b);               # K
  }
  result[12] <- 0;                          # M22

  return(result);
}

OM_U_B <- function( l, b, r)
{
  return( cos(l)*sin(b)/r );             # U
}

OM_V_B <- function( l, b, r)
{
  return( sin(l)*sin(b)/r );             # V
}

OM_W_B <- function( l, b, r)
{
  return( -cos(b)/r );                   # W
}
  
OM_w1_B <- function( l, b, r)
{
  return( sin(l) );                      # W1
}

OM_w2_B <- function( l, b, r)
{
  return( -cos(l) );                     # W2
}
  
OL_B_B <- function(l, b, r)
{
  return(OM_w3_B(l, b, r))
}

OM_w3_B <- function( l, b, r)
{
  return( 0 );                           # W3 = B
}
  
OM_M13_B <- function( l, b, r)
{
  return( cos(2*b)*cos(l) );             # M13
}

OM_M23_B <- function( l, b, r)
{
  return( cos(2*b)*sin(l) );             # M23
}

OL_A_B <- function (l, b, r)
{
  return( OM_M12_B(l, b, r))
}

OM_M12_B <- function( l, b, r)
{
  return( -0.5*sin(2*b)*sin(2*l) );      # M12 = A
}

OM_M11s_B <- function( l, b, r)
{
  return( -0.5*sin(2*b)*cos(l)*cos(l) ); # M11* = M11-M22
}

OM_M33s_B <- function( l, b, r)
{
  return( 0.5*sin(2*b) );               # M33* = M33-M22
}
  
OL_C_B <- function( l, b, r)
{
  return( -0.5*sin(2*b)*cos(2*l) );     # C
}
  
OL_K_B <- function( l, b, r)
{
  return( -0.5*sin(2*b) );               # K
}

OM_M22_B <- function( l, b, r)
{
  return( 0 );                          # M22, не определяется
}
  

#-----------------------------------------------------------------

OL_GxA_L <- function( l, b, r)
{
  return(-sin(abs(b)) * sin(l))
}


OL_GyA_L <- function( l, b, r)
{
  return(-sin(abs(b)) * cos(l))
}

OL_Gx_L <- function( l, b, r)
{
  return(-sin(b) * sin(l))
}


OL_Gy_L <- function( l, b, r)
{
  return(-sin(b) * cos(l))
}

OL_Gx_B <- function( l, b, r)  # not implemented
{
  return(0)
}


OL_Gy_B <- function( l, b, r)  # not implemented
{
  return(0)
}

OL_Gx_R <- function( l, b, r)  # not implemented
{
  return(0)
}


OL_Gy_R <- function( l, b, r)  # not implemented
{
  return(0)
}


#------------------------  Bottlinger model ---------------------

Bottlinger_W0_R <- function( l, b, r, R0, R)
{
  return( 0 )
}

Bottlinger_W0_L <- function( l, b, r, R0, R)
{
  return( -cos(b))
}

Bottlinger_W0_B <- function( l, b, r, R0, R)
{
  return(0)
}

Bottlinger_W1_R <- function( l, b, r, R0, R)
{
  return( R0*(R-R0)*sin(l)*cos(b) )
}

Bottlinger_W1_L <- function( l, b, r, R0, R)
{
  return( (R-R0)*(R0*cos(l) - r*cos(b))/r )
}

Bottlinger_W1_B <- function( l, b, r, R0, R)
{
  return( -R0*(R-R0)*sin(l)*sin(b)/r )
}

Bottlinger_W2_R <- function( l, b, r, R0, R)
{
  return( 0.5*R0*((R-R0)^2)*sin(l)*cos(b))
}

Bottlinger_W2_L <- function( l, b, r, R0, R)
{
  return( 0.5*((R-R0)^2)*(R0*cos(l) - r*cos(b))/r )
}

Bottlinger_W2_B <- function( l, b, r, R0, R)
{
  return( -0.5*R0*((R-R0)^2)*sin(l)*sin(b)/r )
}

Bottlinger_K_R <- function( l, b, r, R0, R)
{
  return(r*cos(b)*cos(b))
}

Bottlinger_K_L <- function( l, b, r, R0, R)
{
  return(0)
}

Bottlinger_K_B <- function( l, b, r, R0, R)
{
  return(-cos(b)*sin(b))
}



#-----------------------------------------------------------------
# stars - matrix(n,3), where 
# stars[,1] - l in degrees, [,2] - b in degrees, [,3] px - kPc
# model - кинематическая модель, 
#   1 - Огородникова-Милна, 
#       type - модификация модели ОМ
#          0 - классический вариант с М11 и М33
#          1 - модифицированная модель Огородникова-Милна, где M11 и М33 заменены на K и С
#   2 - Оорта-Линдблада, 
#          0 - классический вариант с A и B
#          1 - расширенный вариант с С и К
#          2 - расширенный вариант с С, К, |Gx| и |Gy|}
#          3 - расширенный вариант с С, К, Gx и  Gy
#   3 - Эри-Ковальского, 
#   4 - Модель Боттлингера с W(W0), W'(W1), W" (W2) и K. 
# use - массив флагов используемых скоростей
#   1 - mu_l
#   2 - mu_b
#   3 - v_r
MakeOMCoef <- function(stars, use = c(TRUE, TRUE, TRUE), model = 1, type = 0, R0 = 8.0)
{
  n <- nrow(stars)
  
  a1 <- matrix(0, nrow = nrow(stars), ncol = 12)
  colnames(a1) <- c("U", "V", "W","", "", "", "", "", "", "", "", "")
  a2 <- matrix(0, nrow = nrow(stars), ncol = 12)
  a3 <- matrix(0, nrow = nrow(stars), ncol = 12)
  
  if (use[1] == TRUE)
  {
    a1[,1] <- OM_U_L(stars[,1], stars[,2], stars[,3])
    a1[,2] <- OM_V_L(stars[,1], stars[,2], stars[,3])
    a1[,3] <- OM_W_L(stars[,1], stars[,2], stars[,3])
  }
  
  if (use[2] == TRUE){
    a2[,1] <- OM_U_B(stars[,1], stars[,2], stars[,3])
    a2[,2] <- OM_V_B(stars[,1], stars[,2], stars[,3])
    a2[,3] <- OM_W_B(stars[,1], stars[,2], stars[,3])
  }
  
  if (use[3] == TRUE)
  {
    a3[,1] <- OM_U_R(stars[,1], stars[,2], stars[,3])
    a3[,2] <- OM_V_R(stars[,1], stars[,2], stars[,3])
    a3[,3] <- OM_W_R(stars[,1], stars[,2], stars[,3])
  }
  
  if (model == 2)  # Модель Оорта-Линдаблада
  {
    if (use[1] == TRUE)
    {
      #a1[,4] <- OL_B_L(stars[,1], stars[,2], stars[,3])
      a1[,4] <- OM_w3_L(stars[,1], stars[,2], stars[,3])
      #a1[,5] <- OL_A_L(stars[,1], stars[,2], stars[,3])
      a1[,5] <- OM_M12_L(stars[,1], stars[,2], stars[,3])
      colnames(a1)[4:5] <- c("B", "A")
    }
    
    if (use[2] == TRUE){
      #a2[,4] <- OL_B_B(stars[,1], stars[,2], stars[,3])
      a2[,4] <- OM_w3_B(stars[,1], stars[,2], stars[,3])
      #a2[,5] <- OL_A_B(stars[,1], stars[,2], stars[,3])
      a2[,5] <- OM_M12_B(stars[,1], stars[,2], stars[,3])
    }
    
    if (use[3] == TRUE)
    {
      #a3[,4] <- OL_B_R(stars[,1], stars[,2], stars[,3])
      a3[,4] <- OM_w3_R(stars[,1], stars[,2], stars[,3])
      #a3[,5] <- OL_A_R(stars[,1], stars[,2], stars[,3])
      a3[,5] <- OM_M12_R(stars[,1], stars[,2], stars[,3])
    }
    
    if ((type == 1) | (type == 2) | (type == 3))
    {
      if (use[1] == TRUE)
      {
        a1[,6] <- OL_C_L(stars[,1], stars[,2], stars[,3])
        a1[,7] <- OL_K_L(stars[,1], stars[,2], stars[,3])
        colnames(a1)[6:7] <- c("C", "K")
      }
      
      if (use[2] == TRUE){
        a2[,6] <- OL_C_B(stars[,1], stars[,2], stars[,3])
        a2[,7] <- OL_K_B(stars[,1], stars[,2], stars[,3])
      }
      
      if (use[3] == TRUE)
      {
        a3[,6] <- OL_C_R(stars[,1], stars[,2], stars[,3])
        a3[,7] <- OL_K_R(stars[,1], stars[,2], stars[,3])
      } 
    }
    
    if (type == 2)
    {
      if (use[1] == TRUE)
      {
        a1[,8] <- OL_GxA_L(stars[,1], stars[,2], stars[,3])
        a1[,9] <- OL_GyA_L(stars[,1], stars[,2], stars[,3])
        colnames(a1)[8:9] <- c("Gx", "Gy")
      }
      
      if (use[2] == TRUE){
        a2[,8] <- OL_Gx_B(stars[,1], stars[,2], stars[,3])
        a2[,9] <- OL_Gy_B(stars[,1], stars[,2], stars[,3])
      }
      
      if (use[3] == TRUE)
      {
        a3[,8] <- OL_Gx_R(stars[,1], stars[,2], stars[,3])
        a3[,9] <- OL_Gy_R(stars[,1], stars[,2], stars[,3])
      } 
      
    }
    
    if (type == 3)
    {
      if (use[1] == TRUE)
      {
        a1[,8] <- OL_Gx_L(stars[,1], stars[,2], stars[,3])
        a1[,9] <- OL_Gy_L(stars[,1], stars[,2], stars[,3])
        colnames(a1)[8:9] <- c("Gx", "Gy")
      }
      
      if (use[2] == TRUE){
        a2[,8] <- OL_Gx_B(stars[,1], stars[,2], stars[,3])
        a2[,9] <- OL_Gy_B(stars[,1], stars[,2], stars[,3])
      }
      
      if (use[3] == TRUE)
      {
        a3[,8] <- OL_Gx_R(stars[,1], stars[,2], stars[,3])
        a3[,9] <- OL_Gy_R(stars[,1], stars[,2], stars[,3])
      } 
      
    }
    
  } else if (model == 1)  # Модель Огородникова-Милна
  {
    if (use[1] == TRUE)
    {
      a1[,4] <- OM_w1_L(stars[,1], stars[,2], stars[,3])
      a1[,5] <- OM_w2_L(stars[,1], stars[,2], stars[,3])
      a1[,6] <- OM_w3_L(stars[,1], stars[,2], stars[,3])
      a1[,7] <- OM_M13_L(stars[,1], stars[,2], stars[,3])
      a1[,8] <- OM_M23_L(stars[,1], stars[,2], stars[,3])
      a1[,9] <- OM_M12_L(stars[,1], stars[,2], stars[,3])
      colnames(a1)[4:9] <- c("Wx", "Wy", "Wz(B)", "M13", "M23", "M12(A)")
      if (type == 1)
      {
        a1[,10] <- OL_C_L(stars[,1], stars[,2], stars[,3])
        a1[,11] <- OL_K_L(stars[,1], stars[,2], stars[,3])
        colnames(a1)[10:11] <- c("C", "K")
      } else if (type == 0)
      {
        a1[,10] <- OM_M11s_L(stars[,1], stars[,2], stars[,3])
        a1[,11] <- OM_M33s_L(stars[,1], stars[,2], stars[,3])
        colnames(a1)[10:11] <- c("M11*", "M33*")
      }
      a1[,12] <- OM_M22_L(stars[,1], stars[,2], stars[,3])
      colnames(a1)[12] <- "M22"
    }
    
    if (use[2] == TRUE){
      a2[,4] <- OM_w1_B(stars[,1], stars[,2], stars[,3])
      a2[,5] <- OM_w2_B(stars[,1], stars[,2], stars[,3])
      a2[,6] <- OM_w3_B(stars[,1], stars[,2], stars[,3])
      a2[,7] <- OM_M13_B(stars[,1], stars[,2], stars[,3])
      a2[,8] <- OM_M23_B(stars[,1], stars[,2], stars[,3])
      a2[,9] <- OM_M12_B(stars[,1], stars[,2], stars[,3])
      if (type == 1)
      {
        a2[,10] <- OL_C_B(stars[,1], stars[,2], stars[,3])
        a2[,11] <- OL_K_B(stars[,1], stars[,2], stars[,3])
      } else 
      {
        a2[,10] <- OM_M11s_B(stars[,1], stars[,2], stars[,3])
        a2[,11] <- OM_M33s_B(stars[,1], stars[,2], stars[,3])
      }
      a2[,12] <- OM_M22_B(stars[,1], stars[,2], stars[,3])
    }
    
    if (use[3] == TRUE)
    {
      a3[,4] <- OM_w1_R(stars[,1], stars[,2], stars[,3])
      a3[,5] <- OM_w2_R(stars[,1], stars[,2], stars[,3])
      a3[,6] <- OM_w3_R(stars[,1], stars[,2], stars[,3])
      a3[,7] <- OM_M13_R(stars[,1], stars[,2], stars[,3])
      a3[,8] <- OM_M23_R(stars[,1], stars[,2], stars[,3])
      a3[,9] <- OM_M12_R(stars[,1], stars[,2], stars[,3])
      if (type == 1)
      {
        a3[,10] <- OL_C_R(stars[,1], stars[,2], stars[,3])
        a3[,11] <- OL_K_R(stars[,1], stars[,2], stars[,3])
      } else 
      {
        a3[,10] <- OM_M11s_R(stars[,1], stars[,2], stars[,3])
        a3[,11] <- OM_M33s_R(stars[,1], stars[,2], stars[,3])
      }
      a3[,12] <- OM_M22_R(stars[,1], stars[,2], stars[,3])
    } 
    
  } else if (model == 4)  # Модель Боттлингера
  {
    
    R <- sqrt( (stars[,3]*cos(stars[,2]))^2 - 2*R0*stars[,3]*cos(stars[,2])*cos(stars[,1]) + R0^2 )
    
    colnames(a1)[4:7] <- c("W0", "W1", "W2", "K")
    
    if (use[1] == TRUE)
    {
      a1[,4] <- Bottlinger_W0_L(stars[,1], stars[,2], stars[,3], R0, R)
      a1[,5] <- Bottlinger_W1_L(stars[,1], stars[,2], stars[,3], R0, R)
      a1[,6] <- Bottlinger_W2_L(stars[,1], stars[,2], stars[,3], R0, R)
      a1[,7] <- Bottlinger_K_L(stars[,1], stars[,2], stars[,3], R0, R)
    }
    
    if (use[2] == TRUE){
      a2[,4] <- Bottlinger_W0_B(stars[,1], stars[,2], stars[,3], R0, R)
      a2[,5] <- Bottlinger_W1_B(stars[,1], stars[,2], stars[,3], R0, R)
      a2[,6] <- Bottlinger_W2_B(stars[,1], stars[,2], stars[,3], R0, R)
      a2[,7] <- Bottlinger_K_B(stars[,1], stars[,2], stars[,3], R0, R)
    }
    
    if (use[3] == TRUE)
    {
      a3[,4] <- Bottlinger_W1_R(stars[,1], stars[,2], stars[,3], R0, R)
      a3[,5] <- Bottlinger_W1_R(stars[,1], stars[,2], stars[,3], R0, R)
      a3[,6] <- Bottlinger_W2_R(stars[,1], stars[,2], stars[,3], R0, R)
      a3[,7] <- Bottlinger_K_R(stars[,1], stars[,2], stars[,3], R0, R)
    }
      
  }
  
  
  # for (i in 1:n)
  # {
  #   if (use[1] == TRUE)
  #     a1[i,] <- GetOM_L(stars[i,1], stars[i,2], stars[i,3], type)
  #   
  #   if (use[2] == TRUE)
  #     a2[i,] <- GetOM_B(stars[i,1], stars[i,2], stars[i,3], type)
  #   
  #   if (use[3] == TRUE)
  #     a3[i,] <- GetOM_R(stars[i,1], stars[i,2], stars[i,3], type)
  # }
  
  a0 <- matrix(0, nrow = 0, ncol = 12)
  if(use[1]==TRUE)
    a0 <- rbind(a0, a1)
  if(use[2]==TRUE)
    a0 <- rbind(a0, a2)
  if(use[3]==TRUE)
    a0 <- rbind(a0, a3)
  colnames(a0) <- colnames(a1)
  
  # if (model == 1)  #полная модель Огородникова-Милна
  # {
  #   if(use[3] == FALSE)
  #     a0 <- a0[,-12]
  # }  else if (model == 2)  # Модель плоского вращения Оорта-Линдблада
  # {
  #   a0 <- a0[,c(-4, -5, -7, -8, -10:-12)]
  # } else if (model == 3)  # только солнечные члены, уравнение Эри-Ковальского
  # {
  #   a0 <- a0[,1:3]
  # }
  # 
  # if ((use[1] == TRUE) & (use[2] == FALSE) & (use[3] = FALSE))  # если используется только решение по mu_l, выкидываем W
  #   a0 <- a0[,-3]
  
  
  return(RemoveEmptyCols(a0));
  
}

RemoveEmptyCols <- function(a)
{
  for (i in ncol(a):1)
  {
    if (sum(a[,i] != 0) == 0)
      a <- a[,-i]
  }
  return(a)
}

PrepareOMRightSide <- function(stars, use = c(TRUE, TRUE, TRUE))
{
  n <- nrow(stars)
  b <- matrix(0, nrow = 0, ncol = 1)
  if (use[1] == TRUE){
    b1 <- matrix(0, n, 1)
    b1[1:n,1] <- stars[,4]*4.74064;
    b <- rbind(b, b1)
  }
  if (use[2] == TRUE){
    b2 <- matrix(0, n, 1)
    b2[1:n,1] <- stars[,5]*4.74064;
    b <- rbind(b, b2)
  }
  if (use[3] == TRUE){
    b3 <- matrix(0, n, 1)
    b3[1:n,1] <- stars[,6];
    b <- rbind(b, b3)
  }
  return(b);
}


PrepareOMRightSide_old <- function(stars, use_vr = TRUE)
{
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
  return(b)
}

#-----------------------------------------------------------------

Calc_VSF_Coef <- function(stars, use = c(TRUE, FALSE, FALSE), mode = 2, scaling = 0, ef = -1, R0 = 8.0, J = 5)
  # вычисление гармоник сферических функций для заданного множества звезд 
  # stars - матрица положений и скоростей звезд (l, b, px, mu_l, mu_b, v_r)
  # use - флаг (mu_l, mu_b, v_r) - использовать соответствующую скорость или нет
  # mode - способ решения, 1 - TLS через SVD, 2 TLS-LS через собственные числа, см. TLS_Gen()
  # scaling - способ масштабирования, см. TLS_Gen()
  # ef - количество переменных, не содержащих ошибки, см. TLS_Gen
{
  
  #  calculate equation of conditions
  # l, b, px, mu_l, mu_b, vr
  a <- GetSphFuncK_matrix(J, stars[,1], stars[,2])
  
  
  b <- PrepareOMRightSide(stars = stars, use = use)
  
  if (ef == -1)
    ef <-  ncol(a)
  res <- TLS_Gen(a, b, mode = mode, scaling = scaling, ef = ef);
  
  return(res)
}


GetOM_B_SphK <- function(sphK, ver)
{
  if (ver == 1)
  {
    return( sphK[GetJbyNKP(0, 0, 1)+1]/2.784)
  } else if (ver == 2)
  {
    return (sphK[GetJbyNKP(2, 0, 1)+1]/-0.778)
  } else if (ver == 3)
  {
    return (sphK[GetJbyNKP(4, 0, 1)+1]/-0.130)
  } else if (ver == 4)
  {
    return (sphK[GetJbyNKP(6, 0, 1)+1]/-0.049)
  } else if (ver == 5)
  {
    return (sphK[GetJbyNKP(8, 0, 1)+1]/-0.024)
  }
  return(0);
}


GetOM_A_SphK <- function(sphK, ver)
{
  if (ver == 1)
  {
    return( sphK[GetJbyNKP(2, 2, 1)+1]/2.022)
  } else if (ver == 2)
  {
    return (sphK[GetJbyNKP(4, 2, 1)+1]/0.292)
  } else if (ver == 3)
  {
    return (sphK[GetJbyNKP(6, 2, 1)+1]/0.1065461)
  } else if (ver == 4)
  {
    return (sphK[GetJbyNKP(8, 2, 1)+1]/0.05275831)
  } 
  return(0);
}

GetOM_C_SphK <- function(sphK, ver)
{
  if (ver == 1)
  {
    return( sphK[GetJbyNKP(2, 2, 0)+1]/-2.022)
  } else if (ver == 2)
  {
    return (sphK[GetJbyNKP(4, 2, 0)+1]/-0.292)
  } else if (ver == 3)
  {
    return (sphK[GetJbyNKP(6, 2, 0)+1]/-0.1065461)
  } else if (ver == 4)
  {
    return (sphK[GetJbyNKP(8, 2, 0)+1]/-0.05275831)
  } 
  return(0);
}

GetOM_K_SphK <- function(sphK, ver)
{
  return(0)
}
  

#     V_110     V_310
# U/R 0.5530657 -0.5912082
# Gx  0.3259447 -1.3933557

#    V_510     V_710
# U/R 7.326007 -5.999747
# Gx  5.494505 -8.810155

GetOM_UR_SphK <- function(sphK, ver)
{
  if (ver == 1)
  {
    return( 0.5530657*sphK[GetJbyNKP(1, 1, 0)+1] - 0.5912082*sphK[GetJbyNKP(3, 1, 0)+1])
  } else if (ver == 2)
  {
    return( 7.326007*sphK[GetJbyNKP(5, 1, 0)+1] - 5.999747*sphK[GetJbyNKP(7, 1, 0)+1])
    #return( -15.749*sphK[GetJbyNKP(5, 1, 0)+1] + 19.194*sphK[GetJbyNKP(7, 1, 0)+1])
  } 
  
  return(0);
}

#       V_111      V_311
#V/R -0.5530657  0.5912082
#Gy   0.3259447 -1.3933557

#      V_511       V_711
#V/R -7.326007  5.999747
#Gy  5.494505 -8.810155

GetOM_VR_SphK <- function(sphK, ver)
{
  if (ver == 1)
  {
    return( -0.5530657*sphK[GetJbyNKP(1, 1, 1)+1] + 0.5912082*sphK[GetJbyNKP(3, 1, 1)+1])
  } else if (ver == 2)
  {
    return( -7.326007*sphK[GetJbyNKP(5, 1, 1)+1] + 5.999747*sphK[GetJbyNKP(7, 1, 1)+1])
  } 
  
  return(0);
}

GetOM_WR_SphK <- function(sphK, ver)
{
  return(0)
}

GetOM_Gx_SphK <- function(sphK, ver)
{
  if (ver == 1)
  {
    return( 0.3259447*sphK[GetJbyNKP(1, 1, 0)+1] - 1.3933557*sphK[GetJbyNKP(3, 1, 0)+1])
  } else if (ver == 2)
  {
    return( 5.494505*sphK[GetJbyNKP(5, 1, 0)+1] - 8.810155*sphK[GetJbyNKP(7, 1, 0)+1])
    #return( -28.389*sphK[GetJbyNKP(5, 1, 0)+1] + 28.1845*sphK[GetJbyNKP(7, 1, 0)+1])
  } 
  
  return(0);
}

GetOM_Gy_SphK <- function(sphK, ver)
{
  if (ver == 1)
  {
    return( 0.3259447*sphK[GetJbyNKP(1, 1, 1)+1] - 1.3933557*sphK[GetJbyNKP(3, 1, 1)+1])
  } else if (ver == 2)
  {
    return( 5.494505*sphK[GetJbyNKP(5, 1, 1)+1] - 8.810155*sphK[GetJbyNKP(7, 1, 1)+1])
  } 
  
  return(0);
}

GetOm_VSF <- function(sphK, ver)
{
  
}


#-----------------------------------------------------------------
# stars - matrix(n,3), where 
# stars[,1] - l in degrees, [,2] - b in degrees, [,3] px - kPc
# model - кинематическая модель, 
#   1 - Огородникова-Милна, 
#       type - модификация модели ОМ
#          0 - классический вариант с М11 и М33
#          1 - модифицированная модель Огородникова-Милна, где M11 и М33 заменены на K и С
#   2 - Оорта-Линдблада, 
#          0 - классический вариант с A и B
#          1 - расширенный вариант с С и К
#          2 - расширенный вариант с С, К, Gx и  Gy
#   3 - Эри-Ковальского, 

MakeOMCoef_old <- function(stars, use = c(TRUE, TRUE, TRUE), model = 1, type = 0)
{
  n <- nrow(stars)

  a1 <- matrix(0, nrow = nrow(stars), ncol = 12)
  a2 <- matrix(0, nrow = nrow(stars), ncol = 12)
  a3 <- matrix(0, nrow = nrow(stars), ncol = 12)
  
  for (i in 1:n)
  {
    if (use[1] == TRUE)
      a1[i,] <- GetOM_L(stars[i,1], stars[i,2], stars[i,3], type)
    
    if (use[2] == TRUE)
      a2[i,] <- GetOM_B(stars[i,1], stars[i,2], stars[i,3], type)
    
    if (use[3] == TRUE)
      a3[i,] <- GetOM_R(stars[i,1], stars[i,2], stars[i,3], type)
  }

  a0 <- matrix(0, nrow = 0, ncol = 12)
  if(use[1]==TRUE)
    a0 <- rbind(a0, a1)
  if(use[2]==TRUE)
    a0 <- rbind(a0, a2)
  if(use[3]==TRUE)
    a0 <- rbind(a0, a3)
  
  if (model == 1)  #полная модель Огородникова-Милна
  {
    if(use[3] == FALSE)
      a0 <- a0[,-12]
  }  else if (model == 2)  # Модель плоского вращения Оорта-Линдблада
  {
    a0 <- a0[,c(-4, -5, -7, -8, -10:-12)]
  } else if (model == 3)  # только солнечные члены, уравнение Эри-Ковальского
  {
    a0 <- a0[,1:3]
  }


  return(a0);
}


# вычисление параметров заданной кинематической модели
# stars - матрица положений и скоростей звезд (l, b, px, mu_l, mu_b, v_r)
# use - флаг (mu_l, mu_b, v_r) - использовать соответствующую скорость или нет
# mode - способ решения, 1 - TLS через SVD, 2 TLS-LS через собственные числа, см. TLS_Gen()
# scaling - способ масштабирования, см. TLS_Gen()
# ef - количество переменных, не содержащих ошибки, см. TLS_Gen
# model - кинематическая модель, 1 - Огородникова Милна, 2 - Оорта-Линдблада, 3 - Эри-Ковальского, 4 - Боттлингера
# type - вариант модели, см. в MakeOmCoef()
Calc_OM_Model <- function(stars, use = c(TRUE, TRUE, TRUE), mode = 1, scaling = 0, ef = -1, model = 1, type = 0, R0 = 8.0)
{
  #  calculate equation of conditions
  # l, b, px, mu_l, mu_b, vr
  a <- MakeOMCoef(stars = stars, use = use, model = model, type = type, R0 = R0)
  
  b <- PrepareOMRightSide(stars = stars, use = use)
  #b <- rowSums(t(t(a)*GetOM_Default()))

  if (ef == -1)
    ef <-  ncol(a)
  res <- TLS_Gen(a, b, mode, scaling, ef);

  names(res$X) <- colnames(a)
  names(res$s_X) <- paste0("e", colnames(a))
  
  if (model == 1)
  {
  #   if (use_vr == TRUE)
  #   {
  #     if (type == 0)
  #     {
  #       names(res$X) <- c("U", "V", "W","Wx", "Wy", "Wz(B)", "M13", "M23", "M12(A)", "M11*", "M33*", "M22")
  #       names(res$s_X) <- c("eU", "eV", "eW","eWx", "eWy", "eWz(B)", "eM13", "eM23", "eM12(A)", "eM11*", "eM33*", "eM22")
  #     }
  #     else 
  #     {
  #      names(res$X) <- c("U", "V", "W","Wx", "Wy", "Wz(B)", "M13", "M23", "M12(A)", "C", "K", "M22")
  #      names(res$s_X) <- c("eU", "eV", "eW","eWx", "eWy", "eWz(B)", "eM13", "eM23", "eM12(A)", "eC*", "eK*", "eM22")
  #     }
  #   }
  #   else
  #   {
  #     if (type == 0 )
  #     {
  #       names(res$X) <- c("U", "V", "W","Wx", "Wy", "Wz(B)", "M13", "M23", "M12(A)", "M11*", "M33*")
  #       names(res$s_X) <- c("eU", "eV", "eW","eWx", "eWy", "eWz(B)", "eM13", "eM23", "eM12(A)", "eM11*", "eM33*")
  #     } else {
  #       names(res$X) <- c("U", "V", "W","Wx", "Wy", "Wz(B)", "M13", "M23", "M12(A)", "C", "K")
  #       names(res$s_X) <- c("eU", "eV", "eW","eWx", "eWy", "eWz(B)", "eM13", "eM23", "eM12(A)", "eC", "eK")
  #       
  #     }
  #   }
    
    if (type == 0)
    {
      res$Oort <- c(res$X["M12(A)"], res$X["Wz(B)"], 0.5*res$X["M11*"], 0.5*(res$X["M11*"]-2*res$X["M33*"]), NA, NA)
      res$s_Oort <- c(res$s_X["eM12(A)"], res$s_X["eWz(B)"], 0.5*res$s_X["eM11*"], sqrt((0.5*res$s_X["eM11*"])**2 + res$s_X["eM33*"]**2), NA, NA)
    } else
    {
      res$Oort <- c(res$X["M12(A)"], res$X["Wz(B)"], res$X["C"], res$X["K"], NA, NA)
      res$s_Oort <- c(res$s_X["eM12(A)"], res$s_X["eWz(B)"], res$s_X["eC"], res$s_X["eK"], NA, NA)
    }
    
  } else if (model == 2)
  {
    if (type == 0){
      res$Oort <- c(res$X["A"], res$X["B"], NA, NA, NA, NA)
      res$s_Oort <- c(res$s_X["eA"], res$s_X["eB"], NA, NA, NA, NA)  
    } else if (type == 1)
    {
      res$Oort <- c(res$X["A"], res$X["B"], res$X["C"], res$X["K"], NA, NA)
      res$s_Oort <- c(res$s_X["eA"], res$s_X["eB"], res$s_X["eC"], res$s_X["eK"], NA, NA)
    } else if ((type == 2)|(type == 3))
    {
      res$Oort <- c(res$X["A"], res$X["B"], res$X["C"], res$X["K"], res$X["Gx"], res$X["Gy"])
      res$s_Oort <- c(res$s_X["eA"], res$s_X["eB"], res$s_X["eC"], res$s_X["eK"], res$s_X["eGx"], res$s_X["eGy"])
    }
    
  } else if (model == 4)
  {
    res$Oort <- c(-0.5*R0*res$X["W1"], -res$X["W0"] - 0.5*R0*res$X["W1"], NA, res$X["K"], NA, NA)
    res$s_Oort <- c(0.5*R0*res$s_X["eW1"], sqrt(res$s_X["eW0"]^2 + (0.5*R0*res$s_X["eW1"])^2), NA, res$s_X["eK"], NA, NA)
    res$X["A"] <- res$Oort[1]
    res$s_X["eA"] <- res$s_Oort[1]
    res$X["B"] <- res$Oort[2] 
    res$s_X["eB"] <- res$s_Oort[2] 
  } else if (model == 3)
  {
    res$Oort <- c(NA, NA, NA, NA, NA, NA)
    res$s_Oort <- c(NA, NA, NA, NA, NA, NA)
  }
  
  
  names(res$Oort) <- c("A", "B", "C", "K", "Gx", "Gy")
  names(res$s_Oort) <- c("eA", "eB", "eC", "eK", "eGx", "eGy")
  
  # else if (model == 2)
  # {
  #   names(res$X) <- c("U", "V", "W", "Wz(B)", "M12(A)")
  #   names(res$s_X) <- c("eU", "eV", "eW", "eWz(B)", "eM12(A)")
  # } else if(model == 3)
  # {
  #   names(res$X) <- c("U", "V", "W")
  #   names(res$s_X) <- c("eU", "eV", "eW")
  # }

  return(res)
}

#=====================================================================
#----------------------    Test functions    -------------------------
#---------------------------------------------------------------------

GetOM_Default <- function ()
{
  result <- c(10.3, 15.2, 8.0, -2,  1, -15,  -1, 15, 0.5, -0.5, 0.5, -0.5);
  return(result)
}

GetSphCoefDefault <- function(n = 81, U = 10.3, V = 15.2, W = 8.0, A = 15.0, B = -15.0, C = -5.0, Gx = -10.0, Gy = 30.0, r = 1.0)
{
  result <- c(rep(0, max(81,n)))
  result[1:81] <- c( 
   2.784*B,  # 001
   0,        # 101
   2.411*U/r - 1.023*Gx, # 110
   -2.411*V/r - 1.023*Gy, # 111
   -0.778*B, #201
   0,          # 210
   0,          # 211
   -2.022*C,   #220
   2.022*A,   #221
   0,          # 301
   0.564*U/r - 0.957*Gx,    # 310
   -0.564*V/r - 0.957*Gy,   # 311
   0,   # 320
   0,   # 321
   0,   # 330
   0,   # 331
   -0.130*B,   # 401
   0,  # 410
   0,  # 411
   -0.292*C, # 420
   0.292*A, # 421
   0,    # 430
   0,    # 431
   0,    # 440
   0,    # 441
   0,    # 501
   0.279*U/r - 0.190*Gx,    # 510
   -0.279*V/r - 0.190*Gy,    # 511
   0,    # 520 
   0,    # 521 
   0,    # 530    
   0,    # 531 
   0,    # 540 
   0,    # 541
   0,    # 550
   0,    # 551
   -0.049*B,    # 601
   0,    # 610 
   0,    # 611
   0.1065461*C,    # 620 
   -0.1065461*A,    # 621
   0,    # 630 
   0,    # 631
   0,    # 640
   0,    # 641
   0,    # 650
   0,    # 651
   0,    # 660
   0,    # 661
   0,    # 701
   0.174*U/r - 0.232*Gx,    # 710
   -0.174*V/r - 0.232*Gy,    # 711
   0,    # 720
   0,    # 721
   0,    # 730
   0,    # 731
   0,    # 740
   0,    # 741
   0,    # 750
   0,    # 751
   0,    # 760
   0,    # 761
   0,    # 770
   0,    # 771
   -0.024*B,    # 801
   0,    # 810 
   0,    # 811
   0.05275831*C,    # 820 
   -0.05275831*A,    # 821
   0,    # 830 
   0,    # 831
   0,    # 840
   0,    # 841
   0,    # 850
   0,    # 851
   0,    # 860
   0,    # 861
   0,    # 870
   0,    # 871
   0,    # 880
   0    # 881
   ) 
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

Make_Sph_Test <- function()
{
  tgas_sample <- filter_tgs_px(tgas, r_lim = c(950, 1050))
  tgas_sample <- tgas_calc_gpm(tgas_sample)
  stars <- tgas_get_stars(tgas_sample)
  
  
  a0 <- GetSphFuncK_matrix(24, stars[,1], stars[,2])
  
  Sph_0 <- GetSphCoefDefault()
  b0 <- rowSums(t(t(a0)*Sph_0))
  
  res <- TLS_Gen(a0, b0, mode = 1, scaling = 0, ef = -1);

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
