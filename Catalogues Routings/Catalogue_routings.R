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
  hip1 <- hip_data %>% filter(Px > 0.5) %>% filter(Px < 1.5)  %>% filter(B_V<0.75)
  print(paste("stars in sample:", nrow(hip1)))
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
  stars[,4] <- hip_data$pm_l*cos(hip_data$gb)
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
# -------------------------    Tycho-2 Routings   -------------------------------
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
  tyc <- tyc2_data %>% filter(Px > 0.5) %>% filter(Px < 4)  %>% filter(B_V>1.5)
  print(paste("stars in sample:", nrow(tyc)))
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


# ===============================================================================
# -------------------   Tycho-2 Spectral Types Routings  ------------------------
# -------------------------------------------------------------------------------
# Byte-by-byte Description of file: catalog.dat
# --------------------------------------------------------------------------------
#   Bytes Format  Units   Label    Explanations
# -------------- ------------------------------------------------------------------
# 1-  3  A3     ---     ---      [TYC] Tycho-2 label
# 5-  8  I4     ---     TYC1     First part of Tycho-2 identifier
# 10- 14  I5     ---     TYC2     Second part of Tycho-2 identifier
# 16      I1     ---     TYC3     Third part of Tycho-2 identifier
# 18- 29  F12.8  deg     RAdeg    Right Ascension, J2000, decimal deg.
# 31- 42  F12.8  deg     DEdeg    Declination, J2000, decimal deg.
# 44- 49  F6.3   mag     VTmag    ?=99.99 Tycho-2 V_T_ magnitude
# 51- 56  F6.3   mag     BTmag    ?=99.99 Tycho-2 B_T_ magnitude
# 58- 60  A3     ---   r_SpType   Source of spectral type (2)
# 62- 76  A15    ---     Name     Alternate designation for star (3)
# 78- 83  F6.3   arcsec  Dist     Distance between Tycho object and spectral type match (4)
# 85- 90  F6.2   mag     Mag      ?=99.99 Magnitude from SpType catalog (5)
# 92      A1     ---   f_Mag      [VPBX*] Flag indicating type of magnitude (6)
# 94      A1     ---     TClass   Temperature class (7)
# 95      I1     ---     SClass   ? Temperature Subclass (7)
# 97      I1     ---     LClass   ? Luminosity class in numeric form (7)
# 99-103  I5     K       Teff     Effective temperature of the star, based on spectral type (G1)
# 105-124  A20    ---     SpType   Spectral Type (1)
# --------------------------------------------------------------------------------
#   
#   Note (1): This is the spectral type of the star, exactly as it appears
# in the original spectral type catalog.
# 
# Note (2): This column contains a code for the catalog of origin of the
# spectral type:
#   mc1 = Michigan Catalog, Vol. 1, <III/31>
#   mc2 = Michigan Catalog, Vol. 2, <III/51>
#   mc3 = Michigan Catalog, Vol. 3, <III/80>
#   mc4 = Michigan Catalog, Vol. 4, <III/133>
#   mc5 = Michigan Catalog, Vol. 5, <III/214>
#   j64 = Jaschek et al. 1964, <III/18>
#   k83 = Kennedy 1983, <III/78>
#   fI  = FK5, Part I, <I/149>
#   fII = FK5, Part II, <I/175>
#   ppN = PPM North, <I/146>
#   ppS = PPM South, <I/193>
#   sim = SIMBAD Astronomical Database
# 
# Note (3): This is an alternate designation for the star, other than
# its Tycho-2 identifier. It is usually the designation for the star
# that appeared in the spectral type catalog.
# 
# Note (4): This is the distance in arcsec between the Tycho-2 object
# and the star from the spectral type catalog to which it was matched.
# 
# Note (5): This is the magnitude that appears in the spectral type catalog.
# If no magnitude was included, it will have a value of 99.99.
# 
# Note (6): This indicates the type of magnitude that appears in the
# spectral type catalog:
#   V = visual
# P = photographic
# B = blue
# X = unknown
# * = no magnitude was included
# 
# Note (7):
#   In the spectral type reformatting, it was necessary to "choose" a
# concrete spectral type for those that were listed ambiguously, and
# the rule adhered to was to take the first listing. For example, if a
# spectral type was originally listed as K2/3 III, it will be K2 3,
# where K is the temperature class, 2 is the subclass, and 3 is the
# luminosity class.
# 
# Other examples: A9/F2 V, B2.5 V, and G8 IV/V in the original catalog
# become A9 5, B2 5, and G8 4 in the reformatted spectral type
# 
# --------------------------------------------------------------------------------
# -------------------------------------------------------------------------------
# read_tyc2sp - function to read Tycho-2 Spectral Type catalogue into the dataframe
# path - path and file name of the Hipparcos catalogue
# -------------------------------------------------------------------------------
read_tyc2sp <- function(path)
{
  start_pos <- c(1, 5,  10, 16, 18, 31, 44, 51, 58, 62, 78, 85, 92, 94, 95, 
                 97, 99, 105);
  end_pos <- c(3, 8, 14, 16, 29, 42, 49, 56, 60, 76, 83, 90, 92, 94, 95, 97, 
               103, 124);
  var_names <- c("TYC","TYC1", "TYC2", "TYC3", "RA", "DE", "VT", "BT", "SpType_source", "Name", "T2SpDist", "Mag_sp", 
                 "f_Mag_sp", "TClass", "SClass", "LClass", "Teff", "SpType");
  types <- "ciiiddddccddcciiic";
  
  tyc2sp_data <- read_fwf(path, col_positions = fwf_positions(start_pos, end_pos, var_names), col_types = types)
  
  #tyc2sp_data[is.na(tyc2_data$pflag),4] <- "";
  #tyc2_data <- tyc2_data %>% filter(tyc2_data$pflag != "X") %>% 
    #mutate(Px = 1, B_V = (0.850*(BT-VT)), Mag = (VT - 0.090*(BT-VT)), 
     #      RA = RA*pi/180, DE = DE*pi/180) 
  
  # Spectral parallaxes calculation
  
  return(tyc2sp_data)
}

read_tyc2sp_default <- function()
{
  tyc2sp_path <- "X:/Data/Catalogues/Tycho-2 Spectral Type Catalogue/catalog.dat"
  tyc2sp_data <- read_tyc2sp(tyc2sp_path)
}


get_tyc2sp_OB <- function(tyc2sp_data)
{
  tyc <- tyc2sp_data %>% filter(PxSp > 0.5) %>% filter(PxSp < 4)  %>% filter(B_V<0)
  return(tyc)
}


# ===============================================================================
# --------------------------------   TGAS Routings  -----------------------------
# field description: https://gaia.esac.esa.int/documentation/GDR1/datamodel/Ch1/tgas_source.html
# hip, (int), Hipparcos identifier
# tycho2_id, (character, X-Y-Z) Tycho 2 identifier, no spaces, no zeros
# solution_id, long 
# source_id, long 
# random_index, long
# ref_epoch, (double, Time[Julian Years]) Reference epoch 
# ra, (double, Angle[deg]) Right ascension 
# ra_error,  (double, Angle[mas]) Standard error of right ascension
# dec, (double, Angle[deg]) Declination 
# dec_error, (double, Angle[mas]) Standard error of declination 
# parallax, (double, Angle[mas] ) Parallax 
# parallax_error, (double, Angle[mas] ) Standard error of parallax 
# pmra, (double, Angular Velocity[mas/year] ) Proper motion in right ascension direction
# pmra_error, (double, Angular Velocity[mas/year] ) Standard error of proper motion in right ascension direction 
# pmdec, (double, Angular Velocity[mas/year] ) Proper motion in declination direction 
# pmdec_error, (double, Angular Velocity[mas/year] ) Standard error of proper motion in declination direction
# ra_dec_corr, double
# ra_parallax_corr, double
# ra_pmra_corr, double
# ra_pmdec_corr, double
# dec_parallax_corr, double 
# dec_pmra_corr, double 
# dec_pmdec_corr, double 
# parallax_pmra_corr, double 
# parallax_pmdec_corr, double 
# pmra_pmdec_corr,double 
# astrometric_n_obs_al, int
# astrometric_n_obs_ac, int 
# astrometric_n_good_obs_al,int 
# astrometric_n_good_obs_ac, int 
# astrometric_n_bad_obs_al, int 
# astrometric_n_bad_obs_ac,int 
# astrometric_delta_q, (float) Hipparcos/Gaia data discrepancy (Hipparcos subset of TGAS only)
# astrometric_excess_noise, double
# astrometric_excess_noise_sig, double 
# astrometric_primary_flag, (boolean), Primary or seconday
# astrometric_relegation_factor, float 
# astrometric_weight_al, (float, Angle[mbabs-2]) Mean astrometric weight of the source 
# astrometric_weight_ac,  (float, Angle[mbabs-2]) Mean astrometric weight of the source
# astrometric_priors_used,  (int) Type of prior used in the astrometric solution
# matched_observations,  (short) Amount of observations matched to this source
# duplicated_source, (boolean) Source with duplicate sources
# scan_direction_strength_k1, float 
# scan_direction_strength_k2, float 
# scan_direction_strength_k3, float 
# scan_direction_strength_k4, float 
# scan_direction_mean_k1, float 
# scan_direction_mean_k2, float 
# scan_direction_mean_k3, float 
# scan_direction_mean_k4,float 
# phot_g_n_obs,  (int) Number of observations contributing to G photometry
# phot_g_mean_flux, (double, Flux[e-/s]) G-band mean flux 
# phot_g_mean_flux_error, (double, Flux[e-/s]) Error on G-band mean flux 
# phot_g_mean_mag, (double, Magnitude[mag]) G-band mean magnitude 
# phot_variable_flag, (string, Dimensionless[see description]) Photometric variability flag
# l, (double, Angle[deg]) Galactic longitude 
# b, (double, Angle[deg]) Galactic latitude 
# ecl_lon, (double, Angle[deg]) Ecliptic longitude 
# ecl_lat, (double, Angle[deg]) Ecliptic latitude 
# iccccddddddddddd___________________________________dddcdddd
# iccccddddddddddddddddddddd______d______i__________idddcdddd

read_tgas <- function(path, start = 1, n = Inf, is_short = TRUE)
{
  
  if (is_short == TRUE)
  {
     types <- "iccccddddddddddd___________________________________dddcdddd"; 
     var_names <- c("HIP","TYC2","solution_id","source_id","random_index","ref_epoch","RA","ra_error","DE","dec_error",
                    "gPx","parallax_error","pmRA","pmra_error","pmDE","pmdec_error","Gm_flux","phot_g_mean_flux_error",
                    "Gm_mag","phot_variable_flag","l","b","ecl_lon","ecl_lat");
  } else 
  {
     types <- "iccccddddddddddddddddddddd______d______i__________idddcdddd";
     var_names <- c("HIP","TYC2","solution_id","source_id","random_index","ref_epoch","RA","ra_error","DE","dec_error",
                    "gPx","parallax_error","pmRA","pmra_error","pmDE","pmdec_error","ra_dec_corr","ra_parallax_corr","ra_pmra_corr",
                    "ra_pmdec_corr","dec_parallax_corr","dec_pmra_corr","dec_pmdec_corr","parallax_pmra_corr","parallax_pmdec_corr",
                    "pmra_pmdec_corr","astrometric_n_obs_al","astrometric_n_obs_ac","astrometric_n_good_obs_al","astrometric_n_good_obs_ac",
                    "astrometric_n_bad_obs_al","astrometric_n_bad_obs_ac","astrometric_delta_q","astrometric_priors_used",
                    "phot_g_n_obs","Gm_flux","phot_g_mean_flux_error","Gm_mag","phot_variable_flag","l","b","ecl_lon","ecl_lat");
  }
  
  tgas_data <- data.frame()
  
  to_read <- n;
  readed <- 0;
  gi <- 0;
  
  for (i in 0:15)
  {
    if (gi >= (start + n - 1))
      break;
    
    s <- as.character(i);
    if(nchar(s) == 1)
      s <- paste0("0",s);
    filename <- paste0(path, "TgasSource_000-000-0", s,".csv.gz")
    print(paste0("reading: ", filename))
    data <- read_csv(filename, col_names = var_names, col_types = types, skip = 1)
    readed <- nrow(data)
    
    if ((gi + readed) < start)
    {
      print("skip")
      gi <- gi + readed
      next;
    }
      
    if(((start + n)>gi)&((start + n)<(gi+1+readed)))
    {
       print(paste("cut after", as.character(gi), as.character(readed)))
       data <- data[-((start + n - gi):readed),]
    }
    if ((start>gi+1)&(start<=(gi+readed)))
    {
      print(paste("cut before", as.character(gi), as.character(readed)))
      data <- data[-(1:(start - 1 - gi)),]
    }
    gi <- gi + readed

    tgas_data <- rbind(tgas_data, data)
    
  }
  
  tgas_data <- tgas_data %>% mutate(RA = RA*pi/180, DE = DE*pi/180) 

  
  return(tgas_data)
}

read_tgas_default <- function(start = 1, n = Inf, is_short = TRUE)
{
  tgas_path <- "X:/Data/Catalogues/TGAS/"
  tgas_data <- read_tgas(tgas_path, start, n, is_short)
  return (tgas_data)
}

# min_px, max_px - mas
filter_tgs_px <- function(tgs, min_px, max_px)
{
  tgs <- tgs %>% filter(Px > min_px) %>% filter(Px < max_px)
  print(paste("stars in sample:", nrow(tgs)))
  return(tgs)
}

tgas_get_stars <- function(tgas_data)
{
  stars <- matrix(0,nrow(tgas_data), 6)
  stars[,1] <- tgas_data$gl*180/pi
  stars[,2] <- tgas_data$gb*180/pi
  stars[,3] <- (1/tgas_data$Px)    #kPc
  stars[,4] <- tgas_data$pm_l
  stars[,5] <- tgas_data$pm_b
  stars[,6] <- 0
  
  return(stars)
}

tgas_test_OM <- function(tgas_data, min_px = 0, max_px = inf)
{
  tgas_ <- cat_eq2gal(filter_tgs_px(tgas_data, min_px, max_px))
  stars <- tgas_get_stars(tgas_)
  return(stars)
}
