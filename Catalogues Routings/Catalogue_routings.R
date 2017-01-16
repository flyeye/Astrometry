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
#    pm.a = mu*cos(d)       GalPm.l = mu(l)*cos(l)  ?!!!
#    pm.d = mu'             GalPm.b = mu(b)              }
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
  cfi <- (cos(gal[,2])*ci-sin(gal[,2])*si*sin(gal[,1]-L0))/cd;

  GalPm <- cbind( (cfi*pm[,1]+sfi*pm[,2]), (-sfi*pm[,1]+cfi*pm[,2]))

  return(GalPm)
}


#================================================================================
#-------------------------   Hipparcos 2 Routings  ------------------------------
# -------------------------------------------------------------------------------

# Byte-by-byte Description of file: hip2.dat
# -------------------------------------------------------------------------------
#  Bytes Format Units    Label   Explanations
# -------------------------------------------------------------------------------
#    1-  6  I6    ---      HIP     Hipparcos identifier
#    8- 10  I3    ---      Sn      [0,159] Solution type new reduction (1)
#    12  I1    ---      So      [0,5] Solution type old reduction (2)
#    14  I1    ---      Nc      Number of components
#    16- 28 F13.10 rad      RArad   Right Ascension in ICRS, Ep=1991.25
#    30- 42 F13.10 rad      DErad   Declination in ICRS, Ep=1991.25
#    44- 50  F7.2  mas      Plx     Parallax
#    52- 59  F8.2  mas/yr   pmRA    Proper motion in Right Ascension
#    61- 68  F8.2  mas/yr   pmDE    Proper motion in Declination
#    70- 75  F6.2  mas    e_RArad   Formal error on DErad
#    77- 82  F6.2  mas    e_DErad   Formal error on DErad
#    84- 89  F6.2  mas    e_Plx     Formal error on Plx
#    91- 96  F6.2  mas/yr e_pmRA    Formal error on pmRA
#    98-103  F6.2  mas/yr e_pmDE    Formal error on pmDE
#    105-107  I3    ---      Ntr     Number of field transits used
#    109-113  F5.2  ---      F2      Goodness of fit
#    115-116  I2    %        F1      Percentage rejected data
#    118-123  F6.1  ---      var     Cosmic dispersion added (stochastic solution)
#    125-128  I4    ---      ic      Entry in one of the suppl.catalogues
#    130-136  F7.4  mag      Hpmag   Hipparcos magnitude
#    138-143  F6.4  mag    e_Hpmag   Error on mean Hpmag
#    145-149  F5.3  mag      sHp     Scatter of Hpmag
#    151  I1    ---      VA      [0,2] Reference to variability annex
#    153-158  F6.3  mag      B-V     Colour index
#    160-164  F5.3  mag    e_B-V     Formal error on colour index
#    166-171  F6.3  mag      V-I     V-I colour index
#    172-276 15F7.2 ---      UW      Upper-triangular weight matrix (G1)
#--------------------------------------------------------------------------------
#  Note (1): Solution type.
#The solution type is a number 10xd+s consisting of two parts d and s:
#  - s describes the type of solution adopted:
#  1 = stochastic solution (dispersion is given in the 'var' column)
#3 = VIM solution (additional parameters in file hipvim.dat)
#5 = 5-parameter solution (this file)
#7 = 7-parameter solution (additional parameters in hip7p.dat)
#9 = 9-parameter solution (additional parameters in hip9p.dat)
#- d describes peculiarities, as a combination of values:
#  0 = single star
#1 = double star
#2 = variable in the system with amplitude > 0.2mag
#4 = astrometry refers to the photocenter
#8 = measurements concern the secondary (fainter) in the double system

#Note (2): as follows:
#  0 = standard 5-parameter solution
#1 = 7- or 9-parameter solution
#2 = stochastic solution
#3 = double and multiple stars
#4 = orbital binary as resolved in the published catalog
#5 = VIM (variability-induced mover) solution
#--------------------------------------------------------------------------------
  
#--------------------------------------------------------------------------------
# read_hip2 - function to read Hipparcos catalogue into the dataframe
# path - path and file name of the Hipparcos catalogue
#--------------------------------------------------------------------------------
read_hip2 <- function(path)
{
#  fwf_positions(c(1, 6), c(8,10), c(12,12), c(14,14), c(16,28), c(30,42), c(44,50), c(52,59), c(61,68), c(70,75), 
 #               c(77,82), c(84,89), c(91,96), c(98,103), c(105,107), c(109,113), c(115,116), c(118,123), c(125,128), 
#                c(130,136), c(138,143), c(145,149), c(151,151), c(153,158), c(160,164), c(166,171), 
  hip_data <- read_fwf(path, col_positions = fwf_positions(c(1, 8,  12, 14, 16, 30, 44, 52, 61, 70, 77, 84, 91, 98,  105, 109, 115, 118, 125, 130, 138, 145, 151, 153, 160, 166),
                                                           c(6, 10, 12, 14, 28, 42, 50, 59, 68, 75, 82, 89, 96, 103, 107, 113, 116, 123, 128, 136, 143, 149, 151, 158, 164, 171),
                                                           c("HIP", "Sn", "So", "Nc", "RA", "DE", "Px", "pmRA", "pmDE", "e_RA", "e_DE", "e_Px", "e_pmRA", "e_pmDE", 
                                                             "Ntr", "F2", "F1", "var", "ic", "Hpmag", "e_Hpmag", "sHp", "VA", "B_V", "e_B_V", "V_I")), 
                       col_types = "iiiiddddddddddidididddiddd")
  
  return(hip_data)
}

read_hip2_default <- function()
{
  hip_path <- "X:/Data/Catalogues/Hipparcos 2/hip2.dat"
  hip_data <- read_hip2(hip_path)
}


get_hip2_OB <- function(hip_data)
{
  hip1 <- hip_data %>% filter(Px > 0.5) %>% filter(Px < 4) %>% filter(B_V<0)
  return(hip1)
}

hip2_eq2gal<-function(hip_data)
{
  gal <- get_galaxy(cbind(hip_data$RA, hip_data$DE))
  hip_data <- hip_data %>% mutate(gl = gal[,1], gb = gal[,2])
  gal_mu <- get_galaxy_mu(cbind(hip_data$pmRA, hip_data$pmDE), gal, hip_data$DE)
  hip_data <- hip_data %>% mutate(pm_l = gal_mu[,1], pm_b = gal_mu[,2])
  return(hip_data)
}

hip2_get_stars <- function(hip_data)
{
  stars <- matrix(0,nrow(hip_data), 6)
  stars[,1] <- hip_data$gl*180/pi
  stars[,2] <- hip_data$gb*180/pi
  stars[,3] <- 1000/hip_data$Px
  stars[,4] <- hip_data$pm_l
  stars[,5] <- hip_data$pm_b
  stars[,6] <- 0
  
  return(stars)
}
