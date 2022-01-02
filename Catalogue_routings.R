# ===============================================================================
# ---------------------------------  libraries ----------------------------------
cr_libs_required <- function()
{
  #library("tidyverse", lib.loc="~/R/win-library/3.4")
  #library("readxl", lib.loc="~/R/win-library/3.4")
  #library("xlsx", lib.loc="~/R/win-library/3.4")
  #library("stargazer", lib.loc="~/R/win-library/3.4")
  #library("scales", lib.loc="~/R/win-library/3.4")
  #library("gdata", lib.loc="~/R/win-library/3.4")
  #library("xtable", lib.loc="~/R/win-library/3.4")
  #library("gridExtra", lib.loc="~/R/win-library/3.4")
  
  require(readr)
  require(readxl)
  require(tidyverse)
  require(xlsx)
  require(ggplot2)
  require(scales)
  require(stargazer)
  require(gdata)
  require(xtable)
  require(gridExtra)
  require(fst)
  require(data.table)
  require(bit64)
  require(R.utils)
}

lib_install <- function()
{
  install.packages("bit64", "data.table", "fst", "gdata", "ggplot2", "gridExtra", "R.utils", "readr", "readxl", "scales", "stargazer", "tidyverse", "xlsx", "xtable")
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
                                                             "Ntr", "F2", "F1", "var", "ic", "hMag", "e_hMag", "sHp", "VA", "hBV", "e_hBV", "hVI")),
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
  hip1 <- hip_data %>% filter(Px > 0.5) %>% filter(Px < 1.5)  %>% filter(hBV<0.75)
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
  stars[,1] <- hip_data$gl
  stars[,2] <- hip_data$gb
  stars[,3] <- (1/hip_data$Px)    #kPc
  stars[,4] <- hip_data$pm_l #*cos(hip_data$gb)
  stars[,5] <- hip_data$pm_b
  stars[,6] <- 0

  return(stars)
}

hip2_test_OM <- function(hip_data)
{
   hip_ <- get_hip2_OB(hip_data)
   #hip$pmRA <- hip$pmRA*cos(hip$DE)
   hip_ <- cat_eq2gal(hip_)
   stars <- hip2_get_stars(hip_)
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
  tyc2_data <- tyc2_data %>% filter(pflag != "X") %>% filter(pflag != "P")

  tyc2_data <- tyc2_data %>% mutate(Px = 1,
                                    tyc_BV = (0.850*(BT-VT)),
                                    tyc_m = (VT - 0.090*(BT-VT)),
                                    RA = RA*pi/180,
                                    DE = DE*pi/180,
                                    TYC = paste(TYC1, TYC2, TYC3, sep = "-"))

  return(tyc2_data)
}

read_tyc2_default <- function(start = 1, n = Inf)
{
  tyc2_path <- "X:/Data/Catalogues/Tycho-2/CATALOG.DAT"
  tyc2_data <- read_tyc2(tyc2_path, start, n)
}


get_tyc2_OB <- function(tyc2_data)
{
  tyc_ <- tyc2_data %>% filter(B_V>-0.5) %>% filter(B_V<0) %>% filter(Mag>9) %>% filter(Mag<11) # %>% filter(Px > 0.5) %>% filter(Px < 4)
  print(paste("stars in sample:", nrow(tyc_)))
  return(tyc_)
}

tyc2_get_stars <- function(tyc2_data)
{
  stars <- matrix(0,nrow(tyc2_data), 6)
  stars[,1] <- tyc2_data$gl
  stars[,2] <- tyc2_data$gb
  stars[,3] <- (1/tyc2_data$Px)    #kPc
  stars[,4] <- tyc2_data$pm_l
  stars[,5] <- tyc2_data$pm_b
  stars[,6] <- 0

  return(stars)
}

tyc2_test_OM <- function(tyc2_data)
{
  tyc_ <- cat_eq2gal(get_tyc2_OB(tyc2_data))
  stars_ <- tyc2_get_stars(tyc_)
  return(stars_)
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

  tyc2sp_data <- tyc2sp_data %>% mutate(TYC = paste(TYC1, TYC2, TYC3, sep = "-"))

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
  tyc <- tyc2sp_data %>% filter(PxSp > 0.5) %>% filter(PxSp < 4)  %>% filter( B_V<0 )
  return(tyc)
}
# 
# ===============================================================================
# --------------------------------   GAIA DR3 Routings  -----------------------------
# field description: https://gea.esac.esa.int/archive/documentation/GEDR3/Gaia_archive/chap_datamodel/sec_dm_main_tables/ssec_dm_gaia_source.html
#
# solution_id : Solution Identifier (long)
# designation : Unique source designation (unique across all Data Releases) (string)
# source_id: Unique source identifier (unique within a particular Data Release) (long)
# random_index : Random index used to select subsets (long)
# ref_epoch : Reference epoch (double, Time[Julian Years])
# ra : Right ascension (double, Angle[deg])
# ra_error : Standard error of right ascension (float, Angle[mas])
# dec : Declination (double, Angle[deg])
# dec_error : Standard error of declination (float, Angle[mas])
# parallax : Parallax (double, Angle[mas] )
# parallax_error : Standard error of parallax (float, Angle[mas] )
# parallax_over_error : Parallax divided by its standard error (float)
# pm : Total proper motion (float, Angular Velocity[mas/year])
# pmra : Proper motion in right ascension direction (double, Angular Velocity[mas/year])
# pmra_error : Standard error of proper motion in right ascension direction (float, Angular Velocity[mas/year] )
# pmdec : Proper motion in declination direction (double, Angular Velocity[mas/year] )
# pmdec_error : Standard error of proper motion in declination direction (float, Angular Velocity[mas/year] )
# ra_dec_corr : Correlation between right ascension and declination (float, Dimensionless[see description])
# ra_parallax_corr : Correlation between right ascension and parallax (float, Dimensionless[see description])
# ra_pmra_corr : Correlation between right ascension and proper motion in right ascension (float, Dimensionless[see description])
# ra_pmdec_corr : Correlation between right ascension and proper motion in declination (float, Dimensionless[see description])
# dec_parallax_corr : Correlation between declination and parallax (float, Dimensionless[see description])
# dec_pmra_corr : Correlation between declination and proper motion in right ascension (float, Dimensionless[see description])
# dec_pmdec_corr : Correlation between declination and proper motion in declination (float, Dimensionless[see description])
# parallax_pmra_corr : Correlation between parallax and proper motion in right ascension (float, Dimensionless[see description])
# parallax_pmdec_corr : Correlation between parallax and proper motion in declination (float, Dimensionless[see description])
# pmra_pmdec_corr : Correlation between proper motion in right ascension and proper motion in declination (float, Dimensionless[see description])
# astrometric_n_obs_al : Total number of observations AL (short)
# astrometric_n_obs_ac : Total number of observations AC (short)
# astrometric_n_good_obs_al : Number of good observations AL (short)
# astrometric_n_bad_obs_al : Number of bad observations AL (short)
# astrometric_gof_al : Goodness of fit statistic of model wrt along-scan observations (float)
# astrometric_chi2_al : AL chi-square value (float)
# astrometric_excess_noise : Excess noise of the source (float, Angle[mas])
# astrometric_excess_noise_sig : Significance of excess noise (float)
# astrometric_params_solved : Which parameters have been solved for? (byte)
#   00000112=3	- position
#   00001112=7	- pos + px
#   00110112=27	- pos + pm
#   00111112=31	- pos + px + pm
#   01111112=63	- pos + px + pm + vr
#   10111112=95	- pos + px + pm + C
# astrometric_primary_flag : Primary or seconday (boolean)
# nu_eff_used_in_astrometry : Effective wavenumber of the source used in the astrometric solution
# pseudocolour : Astrometrically estimated pseudocolour of the source (float, Misc[??m???1])
# pseudocolour_error : Standard error of the pseudocolour of the source (float, Misc[??m???1])
# ra_pseudocolour_corr : Correlation between right ascension and pseudocolour (float, Dimensionless[see description])
# dec_pseudocolour_corr : Correlation between declination and pseudocolour (float, Dimensionless[see description])
# parallax_pseudocolour_corr : Correlation between parallax and pseudocolour (float, Dimensionless[see description])
# pmra_pseudocolour_corr : Correlation between proper motion in right asension and pseudocolour (float, Dimensionless[see description])
# pmdec_pseudocolour_corr : Correlation between proper motion in declination and pseudocolour (float, Dimensionless[see description])
# astrometric_matched_transits : Matched FOV transits used in the AGIS solution (short)
# visibility_periods_used : Number of visibility periods used in Astrometric solution (short)
# astrometric_sigma5d_max : The longest semi-major axis of the 5-d error ellipsoid (float, Angle[mas])
# matched_transits : The number of transits matched to this source (short)
# new_matched_transits : The number of transits newly incorporated into an existing source in the current cycle (short)
# matched_transits_removed : The number of transits removed from an existing source in the current cycle (short)
# ipd_gof_harmonic_amplitude : Amplitude of the IPD GoF versus position angle of scan (float, Dimensionless[see description])
# ipd_gof_harmonic_phase : Phase of the IPD GoF versus position angle of scan (float, Angle[deg])
# ipd_frac_multi_peak : Percent of successful-IPD windows with more than one peak (byte)
# ipd_frac_odd_win : Percent of transits with truncated windows or multiple gate (byte)
# ruwe : Renormalised unit weight error (float)
# scan_direction_strength_k1 : Degree of concentration of scan directions across the source (float)
# scan_direction_strength_k2 : Degree of concentration of scan directions across the source (float)
# scan_direction_strength_k3 : Degree of concentration of scan directions across the source (float)
# scan_direction_strength_k4 : Degree of concentration of scan directions across the source (float)
# scan_direction_mean_k1 : Mean position angle of scan directions across the source (float, Angle[deg])
# scan_direction_mean_k2 : Mean position angle of scan directions across the source (float, Angle[deg])
# scan_direction_mean_k3 : Mean position angle of scan directions across the source (float, Angle[deg])
# scan_direction_mean_k4 : Mean position angle of scan directions across the source (float, Angle[deg])
# duplicated_source : Source with multiple source identifiers (boolean)
# phot_g_n_obs : Number of observations contributing to G photometry (short)
# phot_g_mean_flux : G-band mean flux (double, Flux[e-/s])
# phot_g_mean_flux_error : Error on G-band mean flux (float, Flux[e-/s])
# phot_g_mean_flux_over_error : G-band mean flux divided by its error (float)
# phot_g_mean_mag : G-band mean magnitude (float, Magnitude[mag])
# phot_bp_n_obs : Number of observations contributing to BP photometry (short)
# phot_bp_mean_flux : Integrated BP mean flux (double, Flux[e-/s])
# phot_bp_mean_flux_error : Error on the integrated BP mean flux (float, Flux[e-/s])
# phot_bp_mean_flux_over_error : Integrated BP mean flux divided by its error (float)
# phot_bp_mean_mag : Integrated BP mean magnitude (float, Magnitude[mag])
# phot_rp_n_obs : Number of observations contributing to RP photometry (short)
# phot_rp_mean_flux : Integrated RP mean flux (double, Flux[e-/s])
# phot_rp_mean_flux_error : Error on the integrated RP mean flux (float, Flux[e-/s])
# phot_rp_mean_flux_over_error : Integrated RP mean flux divided by its error (float)
# phot_rp_mean_mag : Integrated RP mean magnitude (float, Magnitude[mag])
# phot_bp_n_contaminated_transits : Number of BP contaminated transits (short)
# phot_bp_n_blended_transits : Number of BP blended transits (short)
# phot_rp_n_contaminated_transits : Number of RP contaminated transits (short)
# phot_rp_n_blended_transits : Number of RP blended transits (short)
# phot_proc_mode : Photometry processing mode (byte)
# phot_bp_rp_excess_factor : BP/RP excess factor (float)
# bp_rp : BP - RP colour (float, Magnitude[mag])
# bp_g : BP - G colour (float, Magnitude[mag])
# g_rp : G - RP colour (float, Magnitude[mag])
# dr2_radial_velocity : Radial velocity from Gaia DR2 (float, Velocity[km/s] )
# dr2_radial_velocity_error : Radial velocity error from Gaia DR2 (float, Velocity[km/s] )
# dr2_rv_nb_transits : Number of transits used to compute radial velocity in Gaia DR2 (short)
# dr2_rv_template_teff : Teff of the template used to compute radial velocity in Gaia DR2 (float,Temperature[K])
# dr2_rv_template_logg : logg of the template used to compute radial velocity in Gaia DR2 (float, GravitySurface[log cgs])
# dr2_rv_template_fe_h : Fe/H of the template used to compute radial velocity in Gaia DR2 (float, Abundances[dex])
# l : Galactic longitude (double, Angle[deg])
# b : Galactic latitude (double, Angle[deg])
# ecl_lon : Ecliptic longitude (double, Angle[deg])
# ecl_lat : Ecliptic latitude (double, Angle[deg])


GDR3_distance_prepare <- function(path = "/Catalogues/Gaia/", filename ="dedr3dis.gz")
{
  var_types <- "c__d__d__i";
  var_names <- c("source_id", "Rg","Rph", "flags")
  
  for (i in 1:15)  
  {
    data.table(fread("/media/flyeye/T7/Gaia Distance/dedr3dis.csv", select = c(1, 4, 7, 10), colClasses = c("integer64", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric","numeric", "numeric", "integer"), header =TRUE, skip = 0, nrows = 100000000)) 
  }
  
}




var_names <- c("source_id","ref_epoch","RA","ra_error","DE","dec_error","gPx","parallax_error","parallax_over_error",
               "gpmRA","pmra_error","gpmDE","pmdec_error","n_good_obs","params","duplicated", "gMag", "bMag", "rMag",
               "Vr", "Vr_error","Teff","FeH","l","b")

Gaia_DR3_prepare <- function(path = "D:/Gaia/")
{

var_types <- "__c_dddddddd_dddd____________i_____i____________________________l____d____d____d_________dd_d_ddd__";
          

var_names <- c("source_id","ref_epoch","RA","ra_error","DE","dec_error","gPx","parallax_error","parallax_over_error",
               "gpmRA","pmra_error","gpmDE","pmdec_error","n_good_obs","params","duplicated", "gMag", "bMag", "rMag",
               "Vr", "Vr_error","Teff","FeH","l","b")

# var_names <- c("solution_id","designation","source_id","random_index","ref_epoch","ra","ra_error","dec","dec_error","parallax","parallax_error",
#                "parallax_over_error","pm","pmra","pmra_error","pmdec","pmdec_error","ra_dec_corr","ra_parallax_corr","ra_pmra_corr","ra_pmdec_corr",
#                "dec_parallax_corr","dec_pmra_corr","dec_pmdec_corr","parallax_pmra_corr","parallax_pmdec_corr","pmra_pmdec_corr","astrometric_n_obs_al",
#                "astrometric_n_obs_ac","astrometric_n_good_obs_al","astrometric_n_bad_obs_al","astrometric_gof_al","astrometric_chi2_al","astrometric_excess_noise",
#                "astrometric_excess_noise_sig","astrometric_params_solved","astrometric_primary_flag","nu_eff_used_in_astrometry","pseudocolour","pseudocolour_error",
#                "ra_pseudocolour_corr","dec_pseudocolour_corr","parallax_pseudocolour_corr","pmra_pseudocolour_corr","pmdec_pseudocolour_corr","astrometric_matched_transits",
#                "visibility_periods_used","astrometric_sigma5d_max","matched_transits","new_matched_transits","matched_transits_removed",
#                "ipd_gof_harmonic_amplitude","ipd_gof_harmonic_phase","ipd_frac_multi_peak","ipd_frac_odd_win","ruwe","scan_direction_strength_k1",
#                "scan_direction_strength_k2","scan_direction_strength_k3","scan_direction_strength_k4","scan_direction_mean_k1","scan_direction_mean_k2",
#                "scan_direction_mean_k3","scan_direction_mean_k4","duplicated_source","phot_g_n_obs","phot_g_mean_flux","phot_g_mean_flux_error",
#                "phot_g_mean_flux_over_error","phot_g_mean_mag","phot_bp_n_obs","phot_bp_mean_flux","phot_bp_mean_flux_error","phot_bp_mean_flux_over_error",
#                "phot_bp_mean_mag","phot_rp_n_obs","phot_rp_mean_flux","phot_rp_mean_flux_error","phot_rp_mean_flux_over_error","phot_rp_mean_mag","phot_bp_n_contaminated_transits",
#                "phot_bp_n_blended_transits","phot_rp_n_contaminated_transits","phot_rp_n_blended_transits","phot_proc_mode","phot_bp_rp_excess_factor",
#                "bp_rp","bp_g","g_rp","dr2_radial_velocity","dr2_radial_velocity_error","dr2_rv_nb_transits","dr2_rv_template_teff","dr2_rv_template_logg",
#                "dr2_rv_template_fe_h","l","b","ecl_lon","ecl_lat")

#path <- "Y:/Gaia/"
#filename <- paste0(path, "GaiaSource_000000-003111.csv.gz")
#data <- read_csv(filename, col_names = var_names, col_types = var_types, skip = 1)
#readed <- nrow(data)
# data2 <- data %>% filter(duplicated==FALSE) %>% filter((params==31) || (params==63) || (params==95)) %>% filter(parallax_over_error>0) %>% filter(parallax_over_error <= 1) %>% filter(gPx>0.2)

  gaia_list <- read_csv(paste0(path, "gaia.list"), col_names = c("filename"), col_types = "c", skip = 0)


  readed_all <- 0;
  filtred_all <-0;
  written_all <- 0;
  vol_count <- 0;
  
  gdr3_data <- data.table()

  t0 <- Sys.time();
  
  for (i in 1:nrow(gaia_list))
  {
    t1 = Sys.time();
    gc()
    filename <- paste0(path, gaia_list$filename[i])
    cat("step:", i, "file:", filename, "\n")
    data <- data.table(read_csv(filename, col_names = var_names, col_types = var_types, skip = 1))
    readed <- nrow(data)
    cat("readed: ", readed, "\n")
    readed_all <- readed_all + readed
    cat("Total read:", readed_all, "\n")
    data <- data %>% filter(duplicated==FALSE) %>% filter((params==31) || (params==63) || (params==95)) %>% filter(parallax_over_error>0) %>% filter(parallax_over_error > 2) %>% filter(gPx>0.1)  
    filtred <- nrow(data)
    cat("filtred:", filtred, "\n")
    filtred_all <- filtred_all +filtred
    cat("Total filtred:", filtred_all, "\n")
    if (filtred>0)
    {
      data <- data %>% mutate(RA = RA*pi/180, DE = DE*pi/180)
      data <- data %>% mutate(l = l*pi/180, b = b*pi/180)
      data$source_id <- as.integer64(data$source_id)
      #saveRDS(data2, file = paste0(path, tools::file_path_sans_ext(tools::file_path_sans_ext(gaia_list$filename[i])), ".rds"))
      #write_csv(x = data2, path = paste0(path, "gaia_dr3_good.csv"), append = TRUE, col_names = FALSE)
      gdr3_data <- rbind(gdr3_data, data)
      cat("Catalogue size (Gb):", object.size(gdr3_data)/1000000000, "\n")
    }
    
    if (object.size(gdr3_data)>10000000000)
    {
      vol_count <- vol_count+1
      filename <- paste0(path, "GaiaSource_vol", as.character(vol_count),".fst")
      cat("Writing volume ", vol_count, " to ", filename,"\n")
      # saveRDS(gdr3_data, file = filename);
      write.fst(x = gdr3_data, filename, 75)
      written_all <- written_all + nrow(gdr3_data)
      gdr3_data <- data.table();
    }
    
    t2 = Sys.time();
    cat("---------------",i/nrow(gaia_list),",",(Sys.time()-t0),"-----------------\n")
    cat(t2-t1, "\n\n")
  }
  
  if (nrow(gdr3_data)>0)
  {
    vol_count <- vol_count+1
    filename <- paste0(path, "GaiaSource_vol", as.character(vol_count),".fst")
    cat("Writing volume ", vol_count, " to ", filename,"\n")
    # saveRDS(gdr3_data, file = filename);
    write.fst(x = gdr3_data, filename, 75)
    gdr3_data <- data.table();
  }
  
  cat("Read total:", readed_all, "\n")
  cat("Write total:", written_all, "\n")
  
}

# складывает маленькие тома в большие, пока не нужна
GDR3_list_vol <- function(path = "/VMstorage/Gaia/", list_name = "gaia_vol_fst.list")
{
  gaia_list <- read_csv(paste0(path, list_name), col_names = c("filename"), col_types = "c", skip = 0)
  
  gdr3_data <- data.table()
  rdata <- data.table()
  
  readed_all <- 0;
  filtred_all <-0;
  written_all <- 0;
  
  vol_count <- 0;
  
  t0 <- Sys.time();
  
  for (i in 1:nrow(gaia_list))
    {
        t1 = Sys.time();
        gc()
        filename <- paste0(path, gaia_list$filename[i])
        cat("step:", i, "file:", filename, "\n")
        
        data <- data.table(readRDS(file = filename))
        
        readed <- nrow(data)
        cat("readed: ", readed, "\n")
        
        readed_all <- readed_all + readed
        cat("Total read:", readed_all, "\n")
        
        gdr3_data <- rbind(gdr3_data, data)
        cat("Catalogue size (Gb):", object.size(gdr3_data)/1000000000, "\n")

        if (object.size(gdr3_data)>10000000000)
        {
          vol_count <- vol_count+1
          filename <- paste0(path, "GaiaSource_vol", as.character(vol_count),".fst")
          cat("Writing volume ", vol_count, " to ", filename,"\n")
          # saveRDS(gdr3_data, file = filename);
          write.fst(x = gdr3_data, filename, 75)
          written_all <- written_all + nrow(gdr3_data)
          gdr3_data <- data.table();
        }

        # #rdata0 <- data.frame(cbind(data$RA, data$DE, data$gPx, data$parallax_error, data$parallax_over_error, data$l, data$b, data$random_index))
        # #rdata <- rbind(rdata, rdata0)
        # #cat("Sample size (Gb):", object.size(rdata)/1000000000, "\n")
      
        t2 = Sys.time();
        cat("---------------",i/nrow(gaia_list),",",(Sys.time()-t0),"-----------------\n")
        cat(t2-t1, "\n\n")
        
  }
  
  if (nrow(gdr3_data)>0)
  {
    vol_count <- vol_count+1
    filename <- paste0(path, "GaiaSource_vol", as.character(vol_count),".fst")
    cat("Writing volume ", vol_count, " to ", filename,"\n")
    # saveRDS(gdr3_data, file = filename);
    write.fst(x = gdr3_data, filename, 75)
    gdr3_data <- data.table();
  }
  
  cat("Read total:", readed_all, "\n")
  cat("Write total:", written_all, "\n")
  
  return(0)
}

# перечитывает каталог в больших томах и отбирает нужные столбцы
GDR3_list_FST <- function(path = "/Catalogues/Gaia/", list_name = "gaia_vol_fst.list")
{
  gaia_list <- read_csv(paste0(path, list_name), col_names = c("filename"), col_types = "c", skip = 0)
  
  rdata <- data.table(source_id=integer64(), RA=numeric(), DE=numeric(), gPx=numeric(), parallax_error=numeric(), parallax_over_error=numeric(), 
                      gpmRA=numeric(), pmra_error=numeric(), gpmDE=numeric(), pmdec_error=numeric(), Vr=numeric(), Vr_error=numeric(),
                      Teff=numeric(), FeH=numeric(), 
                      gl=numeric(), gb=numeric(), 
                      gMag=numeric(), bMag=numeric())
  
  readed_all <- 0;
  filtred_all <-0;
  written_all <- 0;
  
  vol_count <- 0;
  
  t0 <- Sys.time();
  
  for (i in 1:nrow(gaia_list))
  {
    t1 = Sys.time();
    gc()
    filename <- paste0(path, gaia_list$filename[i])
    cat("step:", i, "file:", filename, "\n")
    
    data <- data.table(read.fst(filename))
    colnames(data)[colnames(data) == 'l'] <- "gl"
    colnames(data)[colnames(data) == 'b'] <- "gb"
    
    readed <- nrow(data)
    cat("readed: ", readed, "\n")
    
    readed_all <- readed_all + readed
    cat("Total read:", readed_all, "\n")
    t2 = Sys.time();
    cat(t2-t1, "\n\n")
    #gdr3_data <- rbind(gdr3_data, data)
    #cat("Catalogue size (Gb):", object.size(gdr3_data)/1000000000, "\n")
    
    #rdata0 <- data.table(cbind(data$source_id,  data$RA, data$DE, data$gPx, data$parallax_error, data$parallax_over_error, data$l, data$b, data$gMag, data$bMag))
    #cat("Columns selected...")
    #colnames(rdata0) <- c("source_id", "RA", "DE", "gPx", "parallax_error", "parallax_over_error", "gl", "gb", "gMag", "bMag")
    #rdata0 <- rdata0 %>% mutate(gl = gl*pi/180, gb = gb*pi/180)
    #rdata <- rbind(rdata, rdata0)
    data$source_id <- as.integer64(data$source_id)
    rdata <- rbind(rdata, data[,.(source_id, RA, DE, gPx, parallax_error, parallax_over_error, gpmRA, pmra_error, gpmDE, pmdec_error, Vr, Vr_error, Teff, FeH, gl, gb, gMag, bMag)])    
    rm(data)
    cat("New data binded..")
    cat("Sample size (Gb):", object.size(rdata)/1000000000, "\n")
    
    t3 = Sys.time();
    cat(t3-t2, "\n\n")
    cat("---------------",i/nrow(gaia_list),"  ", (Sys.time()-t0),"-----------------\n")
    #print("Total time:",(Sys.time()-t0))
    cat(t3-t1, "\n\n")
    #print(t2-t1, "\n\n")
    
  }

   
  cat("Read total:", readed_all, "\n")
  cat("Write total:", written_all, "\n")
  
  return(rdata)
}


gaia_diag <- function(data)
{
  ggplot(data, aes(x = gMag)) + stat_bin(binwidth = 0.02) + theme_bw() +  scale_x_continuous(breaks = c(1:25))
  ggplot(data, aes(x = bv)) + stat_bin(binwidth = 0.02) + theme_bw() + scale_x_continuous(breaks = seq(-1, 3, 0.25), limits = c(-2, 4))
  
  ggplot(gaia[sample(.N, 1000000),], aes(x = x, y = y)) + 
    stat_bin2d(binwidth = 0.01, aes(fill = ..count..)) + 
    theme_bw() + 
    scale_x_continuous(breaks = seq(-5, 5, 1), limits = c(-5, 5)) + 
    scale_y_continuous(breaks = seq(-5, 5, 1), limits = c(-5, 5)) + 
    scale_fill_gradient(low = "white", high = "black")
  
  ggplot(gaia[sample(.N, 1000000),]) + 
    stat_density_2d(binwidth = 0.01, geom = "raster", aes(x = x, y = y, fill = after_stat(density)), contour = FALSE, alpha = 0.8) + 
    theme_bw() + 
    scale_x_continuous(breaks = seq(-5, 5, 1), limits = c(-6, 6)) + 
    scale_y_continuous(breaks = seq(-5, 5, 1), limits = c(-6, 6)) + 
    scale_fill_gradient(low = "white", high = "black", trans = "log")
  
}


# ==============================================================================
# --------------------------------- APASS9 -------------------------------------
# 
# Byte-by-byte Description of file: apass9.sam
# --------------------------------------------------------------------------------
#   Bytes Format Units     Label   Explanations
# --------------------------------------------------------------------------------
# recno, I, 
#  1- 10  F10.6 deg       RAdeg   Right ascension in decimal degrees (J2000)
# 12- 21  F10.6 deg       DEdeg   Declination in decimal degrees (J2000)
# 23- 27  F5.3  arcsec  e_RAdeg   [0/2.4] RA uncertainty
# 29- 33  F5.3  arcsec  e_DEdeg   [0/2.4] DEC uncertainty
# 35- 44  I10   ---       Field   [20110001/9999988888] Field name
# 46- 48  I3    ---       nobs    [2/387] Number of observed nights
# 50- 53  I4    ---       mobs    [2/3476] Number of images for this field, usually nobs*5
# 55- 60  F6.3  mag       B-V     [-7.5/13]? B-V color index
# 62- 67  F6.3  mag     e_B-V     [0/10.1]? B-V uncertainty
# 69- 74  F6.3  mag       Vmag    [5.5/27.4]? Johnson V-band magnitude
# 76- 81  F6.3  mag     e_Vmag    [0/7]? Vmag uncertainty
# 83      I1    ---   u_e_Vmag    [0/1]? Uncertainty flag on e_Vmag (1)
# 85- 90  F6.3  mag       Bmag    [5.4/27.3]? Johnson B-band magnitude
# 92- 97  F6.3  mag     e_Bmag    [0/10]? Bmag uncertainty
# 99      I1    ---   u_e_Bmag    [0/1]? Uncertainty flag on e_Bmag (1)
# 101-106  F6.3  mag       g'mag   [5.9/24.2]? g'-band AB magnitude, Sloan filter
# 108-113  F6.3  mag     e_g'mag   [0/9.7]? g'mag uncertainty
# 115      I1    ---   u_e_g'mag   [0/1]? Uncertainty flag on e_g'mag (1)
# 117-122  F6.3  mag       r'mag   [5.1/23.9]? r'-band AB magnitude, Sloan filter
# 124-129  F6.3  mag     e_r'mag   [0/6.5]? r'mag uncertainty
# 131      I1    ---   u_e_r'mag   [0/1]? Uncertainty flag on e_r'mag (1)
# 133-138  F6.3  mag       i'mag   [4.2/29.1]? i'-band AB magnitude, Sloan filter
# 140-145  F6.3  mag     e_i'mag   [0/9.6]? i'mag uncertainty
# 147      I1    ---   u_e_i'mag   [0/1]? Uncertainty flag on e_i'mag (1)
# --------------------------------------------------------------------------------
#   Note (1): Uncertainty flag as follows:
#   0 = Standard deviation of N observations
# 1 = Poissonian error (unique observation)
# --------------------------------------------------------------------------------

APASS9_prepare <- function(path = "/VMstorage/APASS9/")
{
  
  var_types <- "iddddiiiddddlddlddlddlddl";
  var_names <- c("recno", "RA", "DE", "RA_err", "DE_err", "Field", "nobs", "mobs", "BV", "e_BV", "Vmag", "eVmag", "uVmag", "Bmag", "eBmag", "uBmag","Gmag", "eGmag", "uGmag", "Rmag", "eRmag", "uRmag", "Imag", "eImag", "uImag")
  
  readed_all <- 0;
  filtred_all <-0;
  written_all <- 0;
  
  apass9_data <- data.table()
  
  t0 <- Sys.time();
  
  filename <- paste0(path, "apass9.csv.gz")
  cat("file:", filename, "\n")
  apass9_data <- data.table(read_csv(filename, col_names = var_names, col_types = var_types, skip = 1))
  readed_all <- nrow(apass9_data)
  cat("readed: ", readed_all, "\n")
  
  apass9_data <- apass9_data %>% filter(!is.na(apass$BV))
  
  filtred_all <- nrow(apass9_data)
  cat("Filtred:", filtred_all, "\n")
  if (filtred_all>0)
  {
    apass9_data <- apass9_data %>% mutate(RA = RA*pi/180, DE = DE*pi/180)
    
    #saveRDS(data2, file = paste0(path, tools::file_path_sans_ext(tools::file_path_sans_ext(gaia_list$filename[i])), ".rds"))
    #write_csv(x = data2, path = paste0(path, "gaia_dr3_good.csv"), append = TRUE, col_names = FALSE)
    cat("Catalogue size (Gb):", object.size(apass9_data)/1000000000, "\n")
    filename <- paste0(filename, ".fst")
    cat("Writing APASS9 to ", filename,"\n")
    write.fst(x = apass9_data, filename, 75)
    written_all <- nrow(apass9_data)
  }
  
  cat("---------------",(Sys.time()-t0),"-----------------\n")
  
  return(apass9_data)
}

# write.fst(apass, "/VMstorage/APASS9/apass9_good.csv.gz.fst", 75)
# apass <- read.fst("/VMstorage/APASS9/apass9_good.csv.gz.fst")

# ===============================================================================
# --------------------------------   EGDR3 vs APASS9 ----------------------------

APASS9_GAIA_cross_match_prepare <- function(path = "/VMstorage/Gaia/")
{
  
  var_types <- "ciidiii";
  var_names <- c("source_id", "clean_apass9_oid", "recno", "ang_dist", "nn", "nm", "xm_flag")
  
  cm_list <- read_csv(paste0(path, "apassbestneighbour.list"), col_names = c("filename"), col_types = "c", skip = 0)
  
  readed_all <- 0;
  filtred_all <-0;
  written_all <- 0;
  vol_count <- 0;
  
  
  cm_data <- data.table()
  
  t0 <- Sys.time();
  
  for (i in 1:nrow(cm_list))
  {
    t1 = Sys.time();
    gc()
    filename <- paste0(path, cm_list$filename[i])
    cat("step:", i, "file:", filename, "\n")
    data <- data.table(read_csv(filename, col_names = var_names, col_types = var_types, skip = 1))
    readed <- nrow(data)
    cat("readed: ", readed, "\n")
    readed_all <- readed_all + readed
    cat("Total read:", readed_all, "\n")
    data <- data %>% filter(nm==0) %>% filter(nn==1)
    filtred <- nrow(data)
    cat("filtred:", filtred, "\n")
    filtred_all <- filtred_all +filtred
    cat("Total filtred:", filtred_all, "\n")
    if (filtred>0)
    {
      cm_data <- rbind(cm_data, data)
      cat("Catalogue size (Gb):", object.size(cm_data)/1000000000, "\n")
    }
    
    t2 = Sys.time();
    cat("---------------",i/nrow(cm_list),",",(Sys.time()-t0),"-----------------\n")
    cat(t2-t1, "\n\n")
  }
  
  if (nrow(cm_data)>0)
  {
    filename <- paste0(path, "APASSBestNeibour.fst")
    cat("Writing to ", filename,"\n")
    write.fst(x = cm_data, filename, 75)
  }
  
  cat("Read total:", readed_all, "\n")
  cat("Write total:", written_all, "\n")
  
  return(cm_data)
}


merge_gaia_apass <- function()
{
  data <- cross_match %>% left_join(gaia[,.(RA, DE, gPx, parallax_error, gl, gb, gMag, bMag, x, y, z) ], by = "source_id")
}

# write.fst(cross_match, paste0(path, "APASSBestNeibour.fst"))

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
     var_names <- c("HIP","TYC","solution_id","source_id","random_index","ref_epoch","RA","ra_error","DE","dec_error",
                    "gPx","parallax_error","gpmRA","pmra_error","gpmDE","pmdec_error","Gm_flux","phot_g_mean_flux_error",
                    "Gm_mag","phot_variable_flag","l","b","ecl_lon","ecl_lat");
  } else
  {
     types <- "iccccddddddddddddddddddddd______d______i__________idddcdddd";
     var_names <- c("HIP","TYC","solution_id","source_id","random_index","ref_epoch","RA","ra_error","DE","dec_error",
                    "gPx","parallax_error","gpmRA","pmra_error","gpmDE","pmdec_error","ra_dec_corr","ra_parallax_corr","ra_pmra_corr",
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

  #print(paste("Total records in catalogue:"), nrow(tgas_data))

  tgas_data <- tgas_data %>% mutate(RA = RA*pi/180, DE = DE*pi/180)


  return(tgas_data)
}

read_tgas_default <- function(start = 1, n = Inf, is_short = TRUE)
{
  tgas_path <- "X:/Data/Catalogues/TGAS/"
  tgas_data <- read_tgas(tgas_path, start, n, is_short)
  return (tgas_data)
}


# ================================================================================================



# ===============================================================================
# --------------------------------   UCAC4 Routings  ----------------------------
# col byte item   fmt unit       explanation                            notes
# -------------------------------------------------------------------------------
# 1  1- 3 ra     I*4 mas        right ascension at  epoch J2000.0 (ICRS) (1)
# 2  5- 8 spd    I*4 mas        south pole distance epoch J2000.0 (ICRS) (1)
# 3  9-10 magm   I*2 millimag   UCAC fit model magnitude                 (2)
# 4 11-12 maga   I*2 millimag   UCAC aperture  magnitude                 (2)
# 5 13    sigmag I*1 1/100 mag  error of UCAC magnitude                  (3)
# 6 14    objt   I*1            object type                              (4)
# 7 15    cdf    I*1            combined double star flag                (5)
#         15 bytes
# 8 16    sigra  I*1 mas        s.e. at central epoch in RA (*cos Dec)   (6)
# 9 17    sigdc  I*1 mas        s.e. at central epoch in Dec             (6)
# 10 18    na1    I*1            total # of CCD images of this star
# 11 19    nu1    I*1            # of CCD images used for this star       (7)
# 12 20    cu1    I*1            # catalogs (epochs) used for proper motions
#         5 bytes
# 13 21-22 cepra  I*2 0.01 yr    central epoch for mean RA, minus 1900
# 14 23-24 cepdc  I*2 0.01 yr    central epoch for mean Dec,minus 1900
# 15 25-26 pmrac  I*2 0.1 mas/yr proper motion in RA*cos(Dec)             (8)
# 16 27-28 pmdc   I*2 0.1 mas/yr proper motion in Dec
# 17 29    sigpmr I*1 0.1 mas/yr s.e. of pmRA * cos Dec                   (9)
# 18 30    sigpmd I*1 0.1 mas/yr s.e. of pmDec                            (9)
#        10 bytes
# 19 31-34 pts_key I*4           2MASS unique star identifier            (10)
# 20 35-36 j_m    I*2 millimag   2MASS J  magnitude
# 21 37-38 h_m    I*2 millimag   2MASS H  magnitude
# 22 39-40 k_m    I*2 millimag   2MASS K_s magnitude
# 23 41    icqflg I*1            2MASS cc_flg*10 + ph_qual flag for J    (11)
# 24 42     (2)   I*1            2MASS cc_flg*10 + ph_qual flag for H    (11)
# 25 43     (3)   I*1            2MASS cc_flg*10 + ph_qual flag for K_s  (11)
# 26 44    e2mpho I*1 1/100 mag  error 2MASS J   magnitude               (12)
# 27 45     (2)   I*1 1/100 mag  error 2MASS H   magnitude               (12)
# 28 46     (3)   I*1 1/100 mag  error 2MASS K_s magnitude               (12)
#        16 bytes
# 29 47-48 apasm  I*2 millimag   B magnitude from APASS                  (13)
# 30 49-50  (2)   I*2 millimag   V magnitude from APASS                  (13)
# 31 51-52  (3)   I*2 millimag   g magnitude from APASS                  (13)
# 32 53-54  (4)   I*2 millimag   r magnitude from APASS                  (13)
# 33 55-56  (5)   I*2 millimag   i magnitude from APASS                  (13)
# 34 57    apase  I*1 1/100 mag  error of B magnitude from APASS         (14)
# 35 58     (2)   I*1 1/100 mag  error of V magnitude from APASS         (14)
# 36 59     (3)   I*1 1/100 mag  error of g magnitude from APASS         (14)
# 37 60     (4)   I*1 1/100 mag  error of r magnitude from APASS         (14)
# 38 61     (5)   I*1 1/100 mag  error of i magnitude from APASS         (14)
# 39 62    gcflg  I*1            Yale SPM g-flag*10  c-flag              (15)
#        16 bytes
# 40 63-66 icf(1) I*4            FK6-Hipparcos-Tycho source flag         (16)
# 41       icf(2) ..             AC2000       catalog match flag         (17)
# 42       icf(3) ..             AGK2 Bonn    catalog match flag         (17)
# 43       icf(4) ..             AKG2 Hamburg catalog match flag         (17)
# 44       icf(5) ..             Zone Astrog. catalog match flag         (17)
# 45       icf(6) ..             Black Birch  catalog match flag         (17)
# 46       icf(7) ..             Lick Astrog. catalog match flag         (17)
# 47       icf(8) ..             NPM  Lick    catalog match flag         (17)
# 48       icf(9) ..             SPM  YSJ1    catalog match flag         (17)
#          4 bytes
# 49 67    leda   I*1            LEDA galaxy match flag                  (18)
# 50 68    x2m    I*1            2MASS extend.source flag                (19)
# 51 69-72 rnm    I*4            unique star identification number       (20)
# 52 73-74 zn2    I*2            zone number of UCAC2 (0 = no match)     (21)
# 53 75-78 rn2    I*4            running record number along UCAC2 zone  (21)
#           12 bytes
# ---------------------------------------------------------------------------
#          78 = total number of bytes per star record


make_ucac4_df <- function(num = 0)
{
  res <-data.frame(ra = integer(num), spd = integer(num), u_magm = integer(num), u_maga = integer(num), smag = integer(num),
                   objt=integer(num), cdf=integer(num), sigra=integer(num), sigdc=integer(num),
                   na1=integer(num), nu1=integer(num),cu1=integer(num),
                   cepra=integer(num), cepdc=integer(num), pmrac=integer(num), pmdc=integer(num), sigpmra=integer(num), sigpmdc=integer(num),
                   pts_key=integer(num), j_m=integer(num), h_m=integer(num), k_m=integer(num),
                   icqflag_1=integer(num),icqflag_2=integer(num), icqflag_3=integer(num),
                   e2mpho_1=integer(num), e2mpho_2=integer(num),  e2mpho_3=integer(num),
                   apasm_b=integer(num), apasm_v=integer(num), apasm_g=integer(num), apasm_r=integer(num), apasm_i=integer(num),
                   eapasm_b=integer(num), eapasm_v=integer(num), eapasm_g=integer(num), eapasm_r=integer(num), eapasm_i=integer(num),
                   gcflg=integer(num), icf=integer(num), leda=integer(num), x2m=integer(num), uc4_id=integer(num), zn2=integer(num), rn2=integer(num),
                   isHIP = character(num))
  return(res)
}

read_ucac4_bin_rec <- function(con, num = 1)
{
  if (num == 1)
  {
    res <- vector(length = 45)
    res[1:2] <- readBin(con, "integer", n = 2L, size = 4)
    res[3:4] <- readBin(con, "integer", n = 2L, size = 2)
    res[5:12] <- readBin(con, "integer", n = 8L, size = 1)
    res[13:16] <- readBin(con, "integer", n = 4L, size = 2)
    res[17:18] <- readBin(con, "integer", n = 2L, size = 1)
    res[19] <- readBin(con, "integer", n = 1L, size = 4)
    res[20:22] <- readBin(con, "integer", n = 3L, size = 2)
    res[23:28] <- readBin(con, "integer", n = 6L, size = 1)
    res[29:33] <- readBin(con, "integer", n = 5L, size = 2)
    res[34:38] <- readBin(con, "integer", n = 5L, size = 1)
    res[39] <- readBin(con, "integer", n = 1L, size = 1)
    res[40] <- readBin(con, "integer", n = 1L, size = 4)
    res[41:42] <- readBin(con, "integer", n = 2L, size = 1)
    res[43] <- readBin(con, "integer", n = 1L, size = 4)
    res[44] <- readBin(con, "integer", n = 1L, size = 2)
    res[45] <- readBin(con, "integer", n = 1L, size = 4)
    return(res)
  } else if (num>1)
  {
    rec_size <- 78
    raw <- readBin(con, "raw", n = (num * rec_size))
    
    res <- make_ucac4_df(num)
    res[,1:2] <- t(matrix(readBin(raw[rep(0:(num-1), each=8)*rec_size + 1:8], "integer", n = 2L*num, size = 4), nrow = 2))
    res[,3:4] <- t(matrix(readBin(raw[rep(0:(num-1), each=4)*rec_size + 1:4 + 8], "integer", n = 2L*num, size = 2), nrow = 2))
    res[,5:12] <- t(matrix(readBin(raw[rep(0:(num-1), each=8)*rec_size + 1:8 + 12], "integer", n = 8L*num, size = 1), nrow = 8))
    res[,13:16] <- t(matrix(readBin(raw[rep(0:(num-1), each=8)*rec_size + 1:8 + 20], "integer", n = 4L*num, size = 2), nrow = 4))
    res[,17:18] <- t(matrix(readBin(raw[rep(0:(num-1), each=2)*rec_size + 1:2 + 28], "integer", n = 2L*num, size = 1), nrow = 2))
    res[,19] <- readBin(raw[rep(0:(num-1), each=4)*rec_size + 1:4 + 30], "integer", n = num, size = 4)
    res[,20:22] <- t(matrix(readBin(raw[rep(0:(num-1), each=6)*rec_size + 1:6 + 34], "integer", n = 3L*num, size = 2), nrow = 3))   
    res[,23:28] <- t(matrix(readBin(raw[rep(0:(num-1), each=6)*rec_size + 1:6 + 40], "integer", n = 6L*num, size = 1), nrow = 6))   
    res[,29:33] <- t(matrix(readBin(raw[rep(0:(num-1), each=10)*rec_size + 1:10 + 46], "integer", n = 5L*num, size = 2), nrow = 5))   
    res[,34:38] <- t(matrix(readBin(raw[rep(0:(num-1), each=5)*rec_size + 1:5 + 56], "integer", n = 5L*num, size = 1), nrow = 5))
    res[,39] <- readBin(raw[rep(0:(num-1), each=1)*rec_size + 1 + 61], "integer", n = num, size = 1)
    res[,40] <- readBin(raw[rep(0:(num-1), each=4)*rec_size + 1:4 + 62], "integer", n = num, size = 4)
    res[,41:42] <- t(matrix(readBin(raw[rep(0:(num-1), each=2)*rec_size + 1:2 + 66], "integer", n = 2L*num, size = 1), nrow = 2))
    res[,43] <- readBin(raw[rep(0:(num-1), each=4)*rec_size + 1:4 + 68], "integer", n = num, size = 4)
    res[,44] <- readBin(raw[rep(0:(num-1), each=2)*rec_size + 1:2 + 72], "integer", n = num, size = 2)
    res[,45] <- readBin(raw[rep(0:(num-1), each=4)*rec_size + 1:4 + 74], "integer", n = num, size = 4)
    
    return(res);
  } 
  return (0);
  
}

read_ucac4_rec <- function(con, num, id = NULL, index = NULL)
{
  res <- make_ucac4_df()

  cat("reading...\n")
  
  if (is.null(index))
  {
    t0 <- system.time({res <- read_ucac4_bin_rec(con, num)})
    # t0 <- system.time(
    # for (i in 1:num)
    # {
    #   res[i,1:45] <- read_ucac4_bin_rec(con)
    # }
    # )
    
    cat("reading time per record:", t0/num,"\n")
  } else 
  {
    t1 <- system.time(
    for (i in 1:num)
    {
      seek(con, where = (index[i]-1)*78, rw = "read")
      res[i,1:45] <- read_ucac4_bin_rec(con)
    }
    )
    cat("reading time per record:", t1/num,"\n")
  }
  
  cat("filtring...\n")
  
  t2 <- system.time({
  if( !is.null(id))
  {
    res <- res %>% filter(uc4_id %in% id)
  }

  res <- mutate(res, ra = ra / 3600000, spd = (spd/3600000) - 90, u_magm = u_magm/1000, u_maga = u_magm/1000,
                sigra = sigra +128, sigdc = sigdc +128, sigpmra = sigpmra + 128, sigpmdc = sigpmdc + 128,
                j_m = j_m/1000, h_m = h_m/1000, k_m = k_m/1000, cepra = cepra /100, cepdc = cepdc /100)
  #eapasm_b = eapasm_b/100, eapasm_v = eapasm_v /100, eapasm_g = eapasm_g/100, eapasm_r = eapasm_r/100, eapasm_i = eapasm_i/100
  
  res$apasm_b[res$apasm_b!=20000] <- res$apasm_b[res$apasm_b!=20000] / 1000;
  res$apasm_b[res$apasm_b==20000] <- NA;
  res$apasm_v[res$apasm_v!=20000] <- res$apasm_v[res$apasm_v!=20000] / 1000;
  res$apasm_v[res$apasm_v==20000] <- NA;
  res$apasm_g[res$apasm_g!=20000] <- res$apasm_g[res$apasm_g!=20000] / 1000;
  res$apasm_g[res$apasm_g==20000] <- NA;
  res$apasm_r[res$apasm_r!=20000] <- res$apasm_r[res$apasm_r!=20000] / 1000;
  res$apasm_r[res$apasm_r==20000] <- NA;
  res$apasm_i[res$apasm_i!=20000] <- res$apasm_i[res$apasm_i!=20000] / 1000;
  res$apasm_i[res$apasm_i==20000] <- NA;
  
  res$apasm_b[res$eapasm_b!=99] <- res$eapasm_b[res$eapasm_b!=99] / 100;
  res$apasm_b[res$eapasm_b==99] <- NA;
  res$apasm_v[res$eapasm_v!=99] <- res$eapasm_v[res$eapasm_v!=99] / 100;
  res$apasm_v[res$eapasm_v==99] <- NA;
  res$apasm_g[res$eapasm_g!=99] <- res$eapasm_g[res$eapasm_g!=99] / 100;
  res$apasm_g[res$eapasm_g==99] <- NA;
  res$apasm_r[res$eapasm_r!=99] <- res$eapasm_r[res$eapasm_r!=99] / 100;
  res$apasm_r[res$eapasm_r==99] <- NA;
  res$apasm_i[res$eapasm_i!=99] <- res$eapasm_i[res$eapasm_i!=99] / 100;
  res$apasm_i[res$eapasm_i==99] <- NA;

  res <- mutate(res, isHIP = substr(as.character(icf), 1, 1))
  })
  cat("filtering time per record:",t2/num,"\n")
  return(res)
}

read_ucac4 <- function(path, start = 1, n = Inf, is_tyc_only = TRUE, is_hip = TRUE)
{

  if(is_hip)
  {
    # u4supl.dat = data for supplement stars and cross reference to Hipparcos
    # 
    # This ASCII file contains all UCAC4 stars which have a match to the
    # Hipparcos catalog and all stars which are supplemented, i.e. those
    # without UCAC CCD observations (bright stars).  There are 128631 lines
    # (unique stars) in this file.  These are exactly all those stars with
    # unique star number (col.51 of the main data)  rnm < 1,000,000.
    # 
    # 1         0      0  0  7  0
    # 2         0      0  0  7  0
    # 3         0      0  3  7  0
    # 4         0      0  3  7  0
    # ...
    # 7019         0      0  3  7  0
    # 7020         0      0  3  7  0
    # 7021         0      0  3  7  0
    # 200001         0  43636  2  7  0
    # 200002         0 107949  2  7  0
    # 200003         0 115928  2  7  0
    # ...
    # 321607 181870335  11767  7  6  3
    # 321608 181882163 101884  6  3  0
    # 321609 181883754   3128  7  3  0
    # 321610 181888489  47953  7  3  0
    # 
    # col  explanation
    # 1  rnm  = unique star ID number from UCAC4 (col. 51 main data)
    # 2  MPOS = mean position file (CCD data) internal record number,
    # or zero if not observed with UCAC astrograph
    # 3  hipn = Hipparcos Catalogue star number
    # 4  htsf = FK6-Hipparcos-Tycho source flag (see note 15 main file)
    # 5  objt = object type flag (see note 4 main file)
    # 6  stfl = substitute flag, accumulated from following cases:
    #   0 = no substitute data, use UCAC if available, else external data,
    # else replace UCAC position, proper motion by external data because of:
    #   1 = no "good" image from CCD observation
    #   2 = star is flagged as blended image
    #   4 = position difference to external data is too large (> 50 mas)
    
    cat("Reading UCAC<->HIP indexes... \n")
    start_pos <- c(1, 11, 21, 28, 31, 34);
    end_pos <- c(9, 19, 26, 29, 32, 35);
    var_names <- c("uc4_id", "MPOS", "HIP", "htsf", "odjt", "stfl");
    types <- "iiiiii";
    
    filename <- paste0(path, "u4i/u4supl.dat")
    u4hip_index <- read_fwf(filename, col_positions = fwf_positions(start_pos, end_pos, var_names),
                           col_types = types)
    
    u4hip_index <- u4hip_index[unique(u4hip_index$HIP),]
    
  }
  
  if(is_tyc_only)
  {
    cat("Reading UCAC<->Tycho-2 indexes... \n")
    start_pos <- c(1, 14,  24);
    end_pos <- c(12, 22, 32);
    var_names <- c("TYC", "uc4_index", "uc4_id");
    types <- "cii";

    filename <- paste0(path, "u4i/u4xtycho")
    u4xt_index <- read_fwf(filename, col_positions = fwf_positions(start_pos, end_pos, var_names),
                           col_types = types)

 
    u4xt_index$TYC[(12-nchar(u4xt_index$TYC))==1] <- paste0(" ", u4xt_index$TYC[(12-nchar(u4xt_index$TYC))==1])
    u4xt_index$TYC[(12-nchar(u4xt_index$TYC))==2] <- paste0("  ", u4xt_index$TYC[(12-nchar(u4xt_index$TYC))==2])
    u4xt_index$TYC[(12-nchar(u4xt_index$TYC))==3] <- paste0("   ", u4xt_index$TYC[(12-nchar(u4xt_index$TYC))==3])


    TYC1 <- as.integer(substr(u4xt_index$TYC, 1, 4))
    TYC2 <- as.integer(substr(u4xt_index$TYC, 6, 10))
    TYC3 <- as.integer(substr(u4xt_index$TYC, 12, 12))
    u4xt_index <- mutate(u4xt_index, TYC = paste0(TYC1, "-", TYC2, "-", TYC3),
                         nzone = uc4_index %/% 1000000, sindex = uc4_index %% 1000000)
    
    # UCAC4 содержит около 840 источников, идентифицированных с одиними и теми же зведами каталога Tycho-2, 
    # то есть на один TYC может приходится несколько уникальных записей в UCAC4. Вылоняя ниже left_join
    # мы дублируем записи TGAS, дополеннные разными источниками из UCAC4. Все бы ничего, но TYC и TGAS source_id 
    # теряют при этом уникальность и в последующем при дальнейших оперциях сляиния количество дублей будет множиться. 
    # поэтому по хорошему тут нужно все дубли беспощадно выкинуть. Все равно пока другого способа идентификации звезд
    # TGAS и APASS нет кроме как по TYC, а он вот такой ущербный, соответственно даже если в TGAS и APASS 
    # соответствующие звезд разделены, связать мы их пока никак не можем. 
    # u4xt_index <- u4xt_index[!duplicated(u4xt_index$TYC), ] # !! не тестировалось 
    
  }

  ucac4_data <- make_ucac4_df()
  ucac4_data <- mutate(ucac4_data, TYC = NULL)

  to_read <- n;
  readed <- 0;
  gi <- 0;

  for (i in 1:900)
  {
    if (gi >= (start + n - 1))
      break;

    s <- as.character(i);
    if(nchar(s) == 1){
      s <- paste0("00",s)
    } else if(nchar(s) == 2)
      s <- paste0("0",s);
    filename <- paste0(path, "u4b/","z", s)
    cat(paste0("reading: ", filename), "\n")

    if ( !file.exists(filename))
    {
      stop("catalogue file not found!")
    }


    zone_size <- file.size(filename) %/% 78
    if ((gi + zone_size) < start)
    {
      print("skip")
      gi <- gi + zone_size
      next;
    }

    connection <- file(filename, "rb")
    if (is_tyc_only)
    {
      zone <- filter(u4xt_index, nzone == i)
      #print(zone$sindex)
      data <- read_ucac4_rec(connection, nrow(zone), index = zone$sindex)
      #data <- read_ucac4_rec(connection, zone_size, id = zone$uc4_id)
    } else
    {
      data <- read_ucac4_rec(connection, zone_size)
    }

    readed <- zone_size #nrow(data)
    close(connection)

    cat("adding:", nrow(data), "\n")
    #if(((start + n)>gi)&((start + n)<(gi+1+readed)))
    #{
    #  cat(paste("cut after", as.character(gi), as.character(readed)), "\n")
    #  data <- data[-((start + n - gi):readed),]
    #}
    #if ((start>gi+1)&(start<=(gi+readed)))
    #{
    #  cat(paste("cut before", as.character(gi), as.character(readed)), "\n")
    #  data <- data[-(1:(start - 1 - gi)),]
    #}
    gi <- gi + readed


      
    if (is_hip)
    {
      data <- data %>% left_join(u4hip_index[ , names(u4hip_index) %in% c("HIP", "uc4_id", "stfl")], by = "uc4_id")
      data <- data[!(is.na(data$HIP)),]
      cat("filtered HIP records:", nrow(data), "\n") 
    } else if (is_tyc_only)
    {

      data <- data %>% left_join(u4xt_index[ , names(u4xt_index) %in% c("TYC", "uc4_id")], by = "uc4_id")
      data <- data[!(is.na(data$TYC)),]
      cat("filtered (TYC records):", nrow(rbdata), "\n")
    } 
    
    ucac4_data <- rbind(ucac4_data, data)

    cat("Total records in catalogue:", nrow(ucac4_data), "\n")

  }


  return(ucac4_data)
}

  read_ucac4_default <- function(start = 1, n = Inf, is_tyc_only = TRUE, is_hip = FALSE)
  {
    ucac4_path <- "X:/Data/Catalogues/UCAC4/"
    ucac4_data <- read_ucac4(ucac4_path, start, n, is_tyc_only, is_hip)
    return (ucac4_data)
  }


# ===================================================================================================
#                            URAT1
#-------------------------------------------------------------------------------------------
  #   Byte-by-byte Description of output: urat1.sam
  # -------------------------------------------------------------------------------
  #   Bytes Format Units   Label    Explanations
  # -------------------------------------------------------------------------------
  #   1- 10  A10   ---     URAT1    URAT1 recommended identifier (ZZZ-NNNNNN) (13)
  # 12- 22  F11.7 deg     RAdeg    Right ascension on ICRS, at "Epoch" (1)
  # 24- 34  F11.7 deg     DEdeg    Declination on ICRS, at "Epoch" (1)
  # 36- 38  I3    mas     sigs     Position error per coordinate, from scatter (2)
  # 40- 42  I3    mas     sigm     Position error per coordinate, from model (2)
  # 44- 45  I2    ---     Ns       (nst) Total number of sets the star is in (3)
  # 47- 48  I2    ---     Nu       (nsu) Number of sets used for mean position (3)
  # 50- 57  F8.3  yr      Epoch    (epoc) Mean URAT observation epoch (1)
  # 59- 64  F6.3  mag     f.mag    ?(mmag) mean URAT model fit magnitude (4)
  # 66- 70  F5.3  mag   e_f.mag    ?(sigp) URAT photometry error (5)
  # 72- 73  I2    ---     Nm       (nsm) Number of sets used for URAT magnitude (3)
  # 75  I1    ---     r        (ref) largest reference star flag (6)
  # 77- 79  I3    ---     Nit      (nit) Total number of images (observations)
  # 81- 83  I3    ---     Niu      (niu) Number of images used for mean position
  # 85- 87  I3    ---     Ngt      (ngt) Total number of 1st order grating observations
  # 89- 91  I3    ---     Ngu      (ngu) Number of 1st order grating positions used
  # 93- 98  F6.1  mas/yr  pmRA     ?(pmr) Proper motion RA*cosDec (from 2MASS) (7)
  # 100-105  F6.1  mas/yr  pmDE     ?(pmd) Proper motion in Declination (7)
  # 106-109  F4.1  mas/yr  e_pm     ?(pme) Proper motion error per coordinate (8)
  # 112-113  I2    ---     mf2      [1/11] Match flag URAT with 2MASS (9)
  # 115-116  I2    ---     mfa      [1/11] Match flag URAT with APASS (9)
  # 118  A1    ---     G        [-] "-" if there is no match with GSC2.4 (14)
  # --------------------------------------------------------------------------------
  # 120-129  I10   ---     2Mkey    ?(id2) unique 2MASS star identification number
  # 131-136  F6.3  mag     Jmag     ?(jmag) 2MASS J-band magnitude
  # 138-142  F5.3  mag   e_Jmag     ?(ejmag) Error on Jmag
  # 144-145  A2    ---   q_Jmag     [0,58]? J-band quality-confusion flag (10)
  # 147-152  F6.3  mag     Hmag     ?(hmag) 2MASS H-band magnitude
  # 154-158  F5.3  mag   e_Hmag     ?(ehmag) Error on  H-band magnitude (10)
  # 160-161  A2   ---    q_Hmag     [0,58]? H-band quality-confusion flag (10)
  # 163-168  F6.3  mag     Kmag     ?(kmag) 2MASS Ks-band magnitude
  # 170-174  F5.3  mag   e_Kmag     ?(ekmag) Error on Ks-band magnitude (10)
  # 176-177  A2   ---    q_Kmag     [0,58]? Ks-band quality-confusion flag (10)
  # --------------------------------------------------------------------------------
  # 179-181  I3    ---     Nn       (ann) Number of APASS observation nights (12)
  # 183-185  I3    ---     No       (ano) Number of APASS observations (12)
  # 187-192  F6.3  mag     Bmag     ?(abm) APASS B-band magnitude (11)
  # 194-198  F5.3  mag   e_Bmag     ?(ebm) Error on Bmag
  # 200-205  F6.3  mag     Vmag     ?(avm) APASS V-band magnitude
  # 207-211  F5.3  mag   e_Vmag     ?(evm) Error on Vmag
  # 213-218  F6.3  mag     gmag     ?(agm) APASS g-band magnitude
  # 220-224  F5.3  mag   e_gmag     ?(egm) Error on gmag
  # 226-231  F6.3  mag     rmag     ?(arm) APASS r-band magnitude
  # 233-237  F5.3  mag   e_rmag     ?(erm) Error on rmag
  # 239-244  F6.3  mag     imag     ?(aim) APASS i-band magnitude
  # 246-250  F5.3  mag   e_imag     ?(eim) Error on imag
  # --------------------------------------------------------------------------------
    
  make_urat1_df <- function(num = 0)
  {
    res <-data.frame(
      # ra = integer(num), spd = integer(num), u_magm = integer(num), u_maga = integer(num), smag = integer(num),
      #                objt=integer(num), cdf=integer(num), sigra=integer(num), sigdc=integer(num),
      #                na1=integer(num), nu1=integer(num),cu1=integer(num),
      #                cepra=integer(num), cepdc=integer(num), pmrac=integer(num), pmdc=integer(num), sigpmra=integer(num), sigpmdc=integer(num),
      #                pts_key=integer(num), j_m=integer(num), h_m=integer(num), k_m=integer(num),
      #                icqflag_1=integer(num),icqflag_2=integer(num), icqflag_3=integer(num),
      #                e2mpho_1=integer(num), e2mpho_2=integer(num),  e2mpho_3=integer(num),
      #                apasm_b=integer(num), apasm_v=integer(num), apasm_g=integer(num), apasm_r=integer(num), apasm_i=integer(num),
      #                eapasm_b=integer(num), eapasm_v=integer(num), eapasm_g=integer(num), eapasm_r=integer(num), eapasm_i=integer(num),
      #                gcflg=integer(num), icf=integer(num), leda=integer(num), x2m=integer(num), uc4_id=integer(num), zn2=integer(num), rn2=integer(num),
      #                isHIP = character(num)
      )
    return(res)
  }
    
  read_urat1 <- function(path)
  {
    
    types <- "cdd_____________________________dddddddddd";
    var_names <- c("URAT1_id", RA, DE, apass_v, apass_v_err, apass_b, apass_b_err,apass_g, apass_g_err,apass_r, apass_r_err,apass_i, apass_i_err);
    
    
  }
  
  
# ===================================================================================================
#        TGAS Alternative distance
# -----------------------------------------------------------------------------------------------------------------------------------------------------------------------
#  Bytes     Format    Units     Label       Explanations
# -----------------------------------------------------------------------------------------------------------------------------------------------------------------------
#  1-  6     I6        ---     HIPId       Hipparcos identifier (hip)
#  8- 19     A12       ---     Tycho2      Tycho 2 identifier (tycho2_id)
#  21- 39     I19       ---     SourceId    Source ID (source_id) (G2)
#  41- 58     F18.14    deg     LDeg        Galactic longitude at epoch 2015.0 (l)
#  60- 77     F18.14    deg     BGed        Galactic latitude at epoch 2015.0 (b)
#  79- 96     F18.14    mas     varpi       Absolute barycentric stellar parallax of the source at the reference epoch
#  98-113     F16.14    mas     errVarpi    Standard error of parallax (parallax_error)
#  115-120     F6.3      mag     GMag        G-band mean magnitude
#  122-135     F14.9     pc      rMoExp1     Mode distance of the posterior, using the exponentially decreasing space density prior with L = 0.11 kpc
#  137-150     F14.9     pc      r5Exp1      5th percentile of the posterior, using the exponentially decreasing space density prior with L = 0.11 kpc
#  152-165     F14.9     pc      r50Exp1     50th percentile (i.e. the median) of the posterior, using the exponentially decreasing space density prior with L = 0.11 kpc
#  167-180     F14.9     pc      r95Exp1     95th percentile of the posterior, using the exponentially decreasing space density prior with L = 0.11 kpc
#  182-195     F14.9     pc      sigmaRExp1  Distance standard error, using the exponentially decreasing space density prior with L = 0.11 kpc
#  197-210     F14.9     pc      rMoExp2     Mode distance of the posterior, using the exponentially decreasing space density prior with L = 1.35 kpc      
#  212-225     F14.9     pc      r5Exp2      5th percentile of the posterior, using the exponentially decreasing space density prior with L = 1.35 kpc      
#  227-240     F14.9     pc      r50Exp2     50th percentile (i.e. the median) of the posterior, using the exponentially decreasing space density prior with L = 1.35 kpc      
#  242-255     F14.9     pc      r95Exp2     95th percentile of the posterior, using the exponentially decreasing space density prior with L = 1.35 kpc      
#  257-270     F14.9     pc      sigmaRExp2  Distance standard error, using the exponentially decreasing space density prior with L = 1.35 kpc      
#  272-285     F14.9     pc      rMoMW       Mode distance of the posterior, using Milky Way Prior      
#  287-300     F14.9     pc      r5MW        5th percentile of the posterior, using Milky Way Prior      
#  302-315     F14.9     pc      r50MW       50th percentile (i.e. the median) of the posterior, using Milky Way Prior      
#  317-330     F14.9     pc      r95MW       95th percentile of the posterior, using Milky Way Prior      
#  332-345     F14.9     pc      sigmaRMW    Distance standard error, using Milky Way Prior      
# iccdddd____dd____dd____ddd____d  
  

read_tgas_dist <- function(filename)
{
  if (!file.exists(filename))
    return(NA)
  
  types <- "iccdddddd___dd___dd___d";
  var_names <- c("HIP","TYC","source_id","gl","gb","varpi","errVarpi","GMag",
                 "rMoExp1","sigmaRExp1","rMoExp2","sigmaRExp2","rMoMW","sigmaRMW");
  
  data <- read_csv(filename, col_names = var_names, col_types = types, skip = 1)
  readed <- nrow(data)
  
  return(data)
}
  
  
read_tgas_dist_default  <- function()
{
  filename <- "X:/Data/Catalogues/tgas_dist_all_v01.csv"
  return(read_tgas_dist(filename))
}

add_tgas_dist <- function(tgas_data, dist_data)
{
  
  tgas_data <- tgas_data %>% left_join(dist_data[ , names(dist_data) %in% c("source_id", "varpi","errVarpi","rMoExp1","sigmaRExp1","rMoExp2","sigmaRExp2","rMoMW","sigmaRMW")], by = "source_id")
  
  return(tgas_data)
}
  
# ------------------------------------------------------------------------------
#                                TGAS expanded
# ------------------------------------------------------------------------------

make_tgas_exp <- function(tgas_data, tyc2_data, tyc2_sp_data, hip_data)
{
   tyc2_data <- tyc2_data %>% mutate( TYC = paste(TYC1, TYC2, TYC3, sep = "-"))

   tyc2sp_data <- tyc2sp_data %>% mutate(TYC = paste(TYC1, TYC2, TYC3, sep = "-"))

   # ?? Tycho-2: Mag, B-V
   names(tyc2_data)[7] <- "tyc_pmRA"
   names(tyc2_data)[8] <- "tyc_pmDE"

   tgas_data <- tgas_data %>% left_join(tyc2_data[ , names(tyc2_data) %in% c("TYC", "Mag", "B_V", "tyc_pmRA", "tyc_pmDE")], by = "TYC")
   tgas_data <- tgas_data %>% mutate(tyc_m  = Mag, tyc_bv = B_V)


   # ?? Tycho-2 Spectral Type: LClass, TClass, TSubClass
   tgas_data <- tgas_data %>% left_join(tyc2_sp_data[ , names(tyc2_sp_data) %in% c("TYC", "TClass", "SClass", "LClass", "SpType")], by = "TYC")

   # ?? Hipparcos`?: hPx
   tgas_data <- tgas_data %>% left_join(hip_data[ , names(hip_data) %in% c("HIP", "Px", "e_Px", "hMag", "hBV", "e_hBV", "e_hMag")], by = "HIP")
   
   tgas_data <- tgas_data %>% left_join(ucac[ , names(ucac) %in% c("TYC", "u_magm", "u_maga", "smag", "objt", "cdf",
                                                         "j_m", "h_m", "k_m", "pts_key",
                                                         "apasm_b", "apasm_v", "apasm_g", "apasm_r", "apasm_i", "uc4_id")], by = "TYC")

   
   tgas_data <- tgas_apply_APASS(tgas_data)

}

tgas_apply_APASS <- function(tgas_)
{
  tgas_$B_V <- NA
  tgas_$M <- NA
  tgas_$Mag <- NA
  #index <-(tgas_$gPx>0)&(!is.na(tgas_$apasm_b))&(!is.na(tgas_$apasm_v))
  index <- (!is.na(tgas_$apasm_b))&(!is.na(tgas_$apasm_v))
  tgas_$B_V[index] <- (tgas_$apasm_b[index] - tgas_$apasm_v[index])
  tgas_$Mag[index] <- tgas_$apasm_v[index]
  return(tgas_)
}

tgas_apply_HIP_photometry <- function(tgas_, reset = TRUE)
{
  
  if (reset == TRUE)
  {
    tgas_$B_V <- NA
    tgas_$M <- NA
    tgas_$Mag <- NA
  }
  # index <-(tgas_$gPx>0)&(!is.na(tgas_$hBV))&(!is.na(tgas_$hMag))
  index <- (!is.na(tgas_$hBV))&(!is.na(tgas_$hMag)&(!is.na(tgas_$e_hBV)))
  index2 <- (tgas_$e_hBV!=0)
  tgas_$B_V[index & index2] <- (tgas_$hBV[index & index2])
  tgas_$Mag[index & index2] <- tgas_$hMag[index & index2]
  return(tgas_)
}


fix_TYC <- function(cat)
{
  cat[!is.na(cat$TYC.y),]$TYC.x <- cat[!is.na(cat$TYC.y),]$TYC.y
}
