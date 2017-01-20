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

# -------------------------------------------------------------------------------
cat_eq2gal<-function(cat_data)
{
  gal <- get_galaxy(cbind(cat_data$RA, cat_data$DE))
  cat_data <- cat_data %>% mutate(gl = gal[,1], gb = gal[,2])
  gal_mu <- get_galaxy_mu(cbind(cat_data$pmRA, cat_data$pmDE), gal, cat_data$DE)
  cat_data <- cat_data %>% mutate(pm_l = gal_mu[,1], pm_b = gal_mu[,2])
  return(cat_data)
}


# ===============================================================================
# -------------------------   Hipparcos 2 Routings  -----------------------------
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
                                                             "Ntr", "F2", "F1", "var", "ic", "Mag", "e_Mag", "sHp", "VA", "B_V", "e_B_V", "V_I")), 
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
  hip1 <- hip_data %>% filter(Px > 0.5) %>% filter(Px < 4)  %>% filter(B_V<0)
  return(hip1)
}

#hip2_eq2gal<-function(hip_data)
#{
  #gal <- get_galaxy(cbind(hip_data$RA, hip_data$DE))
  #hip_data <- hip_data %>% mutate(gl = gal[,1], gb = gal[,2])
  #gal_mu <- get_galaxy_mu(cbind(hip_data$pmRA, hip_data$pmDE), gal, hip_data$DE)
  #hip_data <- hip_data %>% mutate(pm_l = gal_mu[,1], pm_b = gal_mu[,2])
  #return(hip_data)
#}

hip2_get_stars <- function(hip_data)
{
  stars <- matrix(0,nrow(hip_data), 6)
  stars[,1] <- hip_data$gl*180/pi
  stars[,2] <- hip_data$gb*180/pi
  stars[,3] <- (1/hip_data$Px)    #kPc
  stars[,4] <- hip_data$pm_l
  stars[,5] <- hip_data$pm_b
  stars[,6] <- 0
  
  return(stars)
}

hip2_test_OM <- function(hip_data)
{
   hip <- cat_eq2gal(get_hip2_OB(hip_data))
   stars <- hip2_get_stars(hip)
   return(stars)
}

# ===============================================================================
# -------------------------   Tycho-2 2 Routings  -------------------------------
# -------------------------------------------------------------------------------
# Byte-by-byte description of file: catalog.dat
# --------------------------------------------------------------------------------
#   Bytes   Format   Units   Label     Explanations
# --------------------------------------------------------------------------------
# 1-  4   I4.4     ---     TYC1      [1,9537]+= TYC1 from TYC or GSC (1)
# 6- 10   I5.5     ---     TYC2      [1,12121]  TYC2 from TYC or GSC (1)
# 12- 13   I1,1X    ---     TYC3      [1,3]      TYC3 from TYC (1)
# 14- 15   A1,1X    ---     pflag     [ PX] mean position flag (2)
# 16- 28   F12.8,1X deg     mRAdeg    []? Mean Right Asc, ICRS, epoch J2000 (3)
# 29- 41   F12.8,1X deg     mDEdeg    []? Mean Decl, ICRS, at epoch J2000 (3)
# 42- 49   F7.1,1X  mas/yr  pmRA*     [-4418.0,6544.2]? prop. mot. in RA*cos(dec)
# 50- 57   F7.1,1X  mas/yr  pmDE      [-5774.3,10277.3]? prop. mot. in Dec
# 58- 61   I3,1X    mas     e_mRA*    [3,183]? s.e. RA*cos(dec),at mean epoch (5)
# 62- 65   I3,1X    mas     e_mDE     [1,184]? s.e. of Dec at mean epoch (5)
# 66- 70   F4.1,1X  mas/yr  e_pmRA*   [0.2,11.5]? s.e. prop mot in RA*cos(dec)(5)
# 71- 75   F4.1,1X  mas/yr  e_pmDE    [0.2,10.3]? s.e. of proper motion in Dec(5)
# 76- 83   F7.2,1X  yr      mepRA     [1915.95,1992.53]? mean epoch of RA (4)
# 84- 91   F7.2,1X  yr      mepDE     [1911.94,1992.01]? mean epoch of Dec (4)
# 92- 94   I2,1X    ---     Num       [2,36]? Number of positions used 
# 95- 98   F3.1,1X  ---     g_mRA     [0.0,9.9]? Goodness of fit for mean RA (6)
# 99-102   F3.1,1X  ---     g_mDE     [0.0,9.9]? Goodness of fit for mean Dec (6)
# 103-106   F3.1,1X  ---     g_pmRA    [0.0,9.9]? Goodness of fit for pmRA (6)
# 107-110   F3.1,1X  ---     g_pmDE    [0.0,9.9]? Goodness of fit for pmDE (6)
# 111-117   F6.3,1X  mag     BT        [2.183,16.581]? Tycho-2 BT magnitude (7)
# 118-123   F5.3,1X  mag     e_BT      [0.014,1.977]? s.e. of BT (7)
# 124-130   F6.3,1X  mag     VT        [1.905,15.193]? Tycho-2 VT magnitude (7)
# 131-136   F5.3,1X  mag     e_VT      [0.009,1.468]? s.e. of VT (7)
# 137-140   I3,1X    ---     prox      [3,999] proximity indicator (8)
# 141-142   A1,1X    ---     TYC       [ T] Tycho-1 star (9)
# 143-148   I6       ---     HIP       [1,120404]? Hipparcos number
# 149-152   A3,1X    ---     CCDM      CCDM component identifier for HIP stars(10)
# 153-165   F12.8,1X deg     RAdeg     Observed Tycho-2 Right Ascension, ICRS
# 166-178   F12.8,1X deg     DEdeg     Observed Tycho-2 Declination, ICRS
# 179-183   F4.2,1X  yr      epRA      [0.81,2.13]  epoch-1990 of RAdeg
# 184-188   F4.2,1X  yr      epDE      [0.72,2.36]  epoch-1990 of DEdeg
# 189-194   F5.1,1X  mas     e_RA*     s.e.RA*cos(dec), of observed Tycho-2 RA (5)
# 195-200   F5.1,1X  mas     e_DE      s.e. of observed Tycho-2 Dec (5)
# 201-202   A1,1X    ---     posflg    [ DP] type of Tycho-2 solution (11)
# 203-206   F4.1     ---     corr      correlation (RAdeg,DEdeg)
# --------------------------------------------------------------------------------
#   Note (1): The TYC identifier is constructed from the GSC region number
# (TYC1), the running number within the region (TYC2) and a component
# identifier (TYC3) which is normally 1. Some non-GSC running numbers 
# were constructed for the first Tycho Catalogue and for Tycho-2. 
# The recommended star designation contains a hyphen between the 
# TYC numbers, e.g. TYC 1-13-1.
# Note (2):
#   ' ' = normal mean position and proper motion.
# 'P' = the mean position, proper motion, etc., refer to the photo-
#   centre of two Tycho-2 entries, where the BT magnitudes 
# were used in weighting the positions.
# 'X' = no mean position, no proper motion.
# Note (3):
#   The mean position is a weighted mean for the catalogues 
# contributing to the proper motion determination. This mean has
# then been brought to epoch 2000.0 by the computed proper motion. 
# See Note(2) for details. Tycho-2 is one of the several catalogues
# used to determine the mean position and proper motion. The 
# observed Tycho-2 position is given in the fields RAdeg and DEdeg.
# Note (4):
#   The mean epochs are given in Julian years.
# Note (5):
#   The errors are based on error models.
# Note (6):
#   This goodness of fit is the ratio of the scatter-based and the 
# model-based error. It is only defined when Num > 2. Values
# exceeding 9.9 are truncated to 9.9.
# Note (7):
#   Blank when no magnitude is available. Either BT or VT is always 
# given. Approximate Johnson photometry may be obtained as:
#   V   = VT -0.090*(BT-VT)
# B-V = 0.850*(BT-VT)
# Consult Sect 1.3 of Vol 1 of "The Hipparcos and Tycho Catalogues",
# ESA SP-1200, 1997, for details.
# Note (8):
#   Distance in units of 100 mas to the nearest entry in the Tycho-2 
# main catalogue or supplement. The distance is computed for the 
# epoch 1991.25. A value of 999 (i.e. 99.9 arcsec) is given if the
# distance exceeds 99.9 arcsec.
# Note (9):
#   ' ' = no Tycho-1 star was found within 0.8 arcsec (quality 1-8)
# or 2.4 arcsec (quality 9).
# 'T' = this is a Tycho-1 star. The Tycho-1 identifier is given in the 
# beginning of the record. For Tycho-1 stars, resolved in
# Tycho-2 as a close pair, both components are flagged as
# a Tycho-1 star and the Tycho-1 TYC3 is assigned to the
# brightest (VT) component.
# The HIP-only stars given in Tycho-1 are not flagged as Tycho-1 stars.
# Note (10):
#   The CCDM component identifiers for double or multiple Hipparcos 
# stars contributing to this Tycho-2 entry. For photocentre 
# solutions, all components within 0.8 arcsec contribute. For double
# star solutions any unresolved component within 0.8 arcsec 
# contributes. For single star solutions, the predicted signal from 
# close stars were normally subtracted in the analysis of the photon 
# counts and such stars therefore do not contribute to the solution.  
# The components are given in lexical order. 
# Note (11):
#   ' ' = normal treatment, close stars were subtracted when possible.
# 'D' = double star treatment. Two stars were found. The companion is
# normally included as a separate Tycho-2 entry, but may have
# been rejected.
# 'P' = photocentre treatment, close stars were not subtracted. This
# special treatment was applied to known or suspected doubles
# which were not successfully (or reliably) resolved in the
# Tycho-2 double star processing.


# -------------------------------------------------------------------------------
# read_tyc2 - function to read Tycho-2 catalogue into the dataframe
# path - path and file name of the Hipparcos catalogue
# -------------------------------------------------------------------------------
read_tyc2 <- function(path, start = 1, n = Inf, is_short = TRUE)
{
  if (is_short == TRUE){
    start_pos <- c(1, 6,  12, 14, 16, 29, 42, 50, 58, 62, 66, 71, 76, 84, 92, 
                   111, 118, 124, 131, 137, 141, 143, 201, 203);
    end_pos <- c(4, 10, 12, 14, 27, 40, 48, 56, 60, 64, 69, 74, 82, 90, 93, 116, 
                 122, 129, 135, 139, 141, 147, 201, 205);
    var_names <- c("TYC1", "TYC2", "TYC3", "pflag", "RA", "DE", "pmRA", "pmDE", "e_RA", "e_DE", "e_pmRA", "e_pmDE", 
                   "mepRA", "mepDE", "Num", "BT", "e_BT", "VT", "e_VT", "prox", "TYC", "HIP", "posflg", "Corr");
    types <- "iiicddddiiddddiddddicicd";
  }
  else   
  {
    start_pos <- c(1, 6,  12, 14, 16, 29, 42, 50, 58, 62, 66, 71, 76, 84, 92, 95, 99, 103, 107, 
                   111, 118, 124, 131, 137, 141, 143, 149, 153, 166, 179, 184, 189, 195, 201, 203);
    end_pos <- c(4, 10, 12, 14, 27, 40, 48, 56, 60, 64, 69, 74, 82, 90, 93, 97, 101, 105, 109, 116, 
               122, 129, 135, 139, 141, 147, 151, 164, 177, 182, 187, 193, 199, 201, 205);
    var_names <- c("TYC1", "TYC2", "TYC3", "pflag", "RA", "DE", "pmRA", "pmDE", "e_RA", "e_DE", "e_pmRA", "e_pmDE", 
                 "mepRA", "mepDE", "Num", "g_RA", "g_DE", "g_pmRA", "g_pmDE", "BT", "e_BT", "VT", "e_VT", "prox", 
                 "TYC", "HIP", "CCDM", "oRA", "oDE", "epRA", "epDE", "e_oRA", "e_oDE", "posflg", "Corr");
    types <- "iiicddddiiddddiddddddddicicddddddcd";
  }
    
  tyc2_data <- read_fwf(path, col_positions = fwf_positions(start_pos,
                                                            end_pos,
                                                            var_names),
                              col_types = types, skip = start-1, n_max = n)
  
  tyc2_data[is.na(tyc2_data$pflag),4] <- "";
  tyc2_data <- tyc2_data %>% filter(tyc2_data$pflag != "X") %>% 
                             mutate(Px = 1, B_V = (0.850*(BT-VT)), Mag = (VT - 0.090*(BT-VT)), 
                                    RA = RA*pi/180, DE = DE*pi/180) 
  
  return(tyc2_data)  
}

read_tyc2_default <- function(start = 1, n = Inf)
{
  tyc2_path <- "X:/Data/Catalogues/Tycho-2/CATALOG.DAT"
  tyc2_data <- read_tyc2(tyc2_path, start, n)
}


get_tyc2_OB <- function(tyc2_data)
{
  tyc <- tyc2_data %>% filter(Px > 0.5) %>% filter(Px < 4)  %>% filter(B_V<0)
  return(tyc)
}

tyc2_get_stars <- function(tyc2_data)
{
  stars <- matrix(0,nrow(tyc2_data), 6)
  stars[,1] <- tyc2_data$gl*180/pi
  stars[,2] <- tyc2_data$gb*180/pi
  stars[,3] <- (1/tyc2_data$Px)    #kPc
  stars[,4] <- tyc2_data$pm_l
  stars[,5] <- tyc2_data$pm_b
  stars[,6] <- 0
  
  return(stars)
}

tyc2_test_OM <- function(tyc2_data)
{
  tyc <- cat_eq2gal(get_tyc2_OB(tyc2_data))
  stars <- tyc2_get_stars(tyc)
  return(stars)
}