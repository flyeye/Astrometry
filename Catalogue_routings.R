# ===============================================================================
# ---------------------------------  libraries ----------------------------------
tgas_libs_required <- function()
{
  library("readr", lib.loc="~/R/win-library/3.3")
  library("readxl", lib.loc="~/R/win-library/3.3")
  library("tidyverse", lib.loc="~/R/win-library/3.3")
  library("xlsx", lib.loc="~/R/win-library/3.3")
  library("stargazer", lib.loc="~/R/win-library/3.3")
  library("scales", lib.loc="~/R/win-library/3.3")
  library("gdata", lib.loc="~/R/win-library/3.3")
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
                                    B_V = (0.850*(BT-VT)),
                                    Mag = (VT - 0.090*(BT-VT)),
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

# min_px, max_px - mas
filter_tgs_px <- function(tgs, px = c(-Inf,Inf), e_px = Inf, bv = c(-Inf,Inf), Mg = c(-Inf,Inf), z_lim = c(0, Inf))
{
  #print(bv)
  #print(Mg)
  #print(px)
  #print(e_px)
  tgs <- CalcGalXYZ(tgs)
  tgs <- tgs %>% filter(!is.na(B_V) & !is.na(M)) %>%
                 filter( (gPx > px[1]) & (gPx < px[2])) %>%
                 filter( (B_V>bv[1]) & (B_V<bv[2]) ) %>%
                 filter( (M>Mg[1]) & (M<Mg[2])) %>%
                 filter( (parallax_error/gPx) < e_px ) %>% 
                 filter( (abs(z)>=z_lim[1]) & (abs(z)<z_lim[2]))

  print(paste("stars in sample:", nrow(tgs)))
  return(tgs)
}

tgas_get_stars <- function(tgas_data, src = "TGAS")
{
  stars <- matrix(0,nrow(tgas_data), 6)
  stars[,1] <- tgas_data$gl
  stars[,2] <- tgas_data$gb
  stars[,3] <- (1/tgas_data$gPx)    #kPc
  if (src == "TGAS")
  {
    stars[,4] <- tgas_data$gpm_l
    stars[,5] <- tgas_data$gpm_b
  } else if(src == "TYCHO")
  {
    stars[,4] <- tgas_data$tpm_l
    stars[,5] <- tgas_data$tpm_b
  } else
  {
    stars[,4] <- tgas_data$pm_l
    stars[,5] <- tgas_data$pm_b
  }
  stars[,6] <- 0

  return(stars)
}

tgas_calc_gpm <- function(tgas_)
{

  tgas_$pmRA <- tgas_$gpmRA
  tgas_$pmDE <- tgas_$gpmDE
  tgas_ <- cat_eq2gal(tgas_)
  tgas_ <- mutate(tgas_, gpm_l = pm_l, gpm_b = pm_b)

  tgas_ <- filter(tgas_ , !is.na(tyc_pmRA))
  tgas_$pmRA <- tgas_$tyc_pmRA
  tgas_$pmDE <- tgas_$tyc_pmDE
  tgas_ <- cat_eq2gal(tgas_)
  tgas_ <- mutate(tgas_, tpm_l = pm_l, tpm_b = pm_b)

  return(tgas_)
}


tgas_export_solution <- function(res, name)
{
  s_ <- name
  
  output <- cbind(res$X[,1], res$S_X[,1], res$X[,2], res$S_X[,2], res$X[,3], res$S_X[,3], res$X[,4], res$S_X[,4], res$X[,5], res$S_X[,5], res$X[,6], res$S_X[,6],
                  res$X[,7], res$S_X[,7], res$X[,8], res$S_X[,8], res$X[,9], res$S_X[,9], res$X[,10], res$S_X[,10], res$X[,11], res$S_X[,11], 
                  res$Oort[,1], res$s_Oort[,1],res$Oort[,2], res$s_Oort[,2],res$Oort[,3], res$s_Oort[,3],res$Oort[,4], res$s_Oort[,4],
                  res$Parameters[,1:5])
  output <- t(output)
  rownames(output)[1:30]<-c(colnames(res$X)[1], colnames(res$S_X)[1], colnames(res$X)[2], colnames(res$S_X)[2], colnames(res$X)[3], colnames(res$S_X)[3],
                            colnames(res$X)[4], colnames(res$S_X)[4], colnames(res$X)[5], colnames(res$S_X)[5], colnames(res$X)[6], colnames(res$S_X)[6],
                            colnames(res$X)[7], colnames(res$S_X)[7], colnames(res$X)[8], colnames(res$S_X)[8], colnames(res$X)[9], colnames(res$S_X)[9],
                            colnames(res$X)[10], colnames(res$S_X)[10], colnames(res$X)[11], colnames(res$S_X)[11], 
                            colnames(res$Oort)[1], colnames(res$s_Oort)[1], colnames(res$Oort)[2], colnames(res$s_Oort)[2], 
                            colnames(res$Oort)[3], colnames(res$s_Oort)[3], colnames(res$Oort)[4], colnames(res$s_Oort)[4])
  write.xlsx2(x = output, file = paste0(s_,".xls"), sheetName = "Solution")
  
  
  output <- cbind(res$X, res$Oort)
  output_err <- cbind(res$S_X, res$s_Oort) 
  
  output_txt <- paste0("$", sprintf("%.2f", output), "pm", sprintf("%.2f", output_err), "$")
  output_txt <- t(matrix(output_txt, ncol = 15))
  output_txt <- rbind(output_txt, paste0(sprintf("%.1f", res$Parameters[,1]), "-", sprintf("%.1f", res$Parameters[,2]) ))
  output_txt <- rbind(output_txt, sprintf("%d", res$Parameters[,3]) )
  output_txt <- rbind(output_txt, sprintf("%.2f", res$Parameters[,4]) ) 
  output_txt <- rbind(output_txt, sprintf("%.2f", res$Parameters[,5]) ) 
  rownames(output_txt) <- c(colnames(res$X), colnames(res$Oort), c("r", "N", "r_mean", "ePx"))
  
  write_lines(stargazer(output_txt, type = "latex", digits = 2), paste0(s_,".tex"), append = FALSE)
  
  
  output_txt <- paste0(sprintf("%.2f", output), "±", sprintf("%.2f", output_err))
  output_txt <- t(matrix(output_txt, ncol = 15))
  output_txt <- rbind(output_txt, paste0(sprintf("%.2f", res$Parameters[,1]), "-", sprintf("%.2f", res$Parameters[,2]) ))
  output_txt <- rbind(output_txt, sprintf("%d", res$Parameters[,3]) )
  output_txt <- rbind(output_txt, sprintf("%.2f", res$Parameters[,4]) ) 
  output_txt <- rbind(output_txt, sprintf("%.2f", res$Parameters[,5]) ) 
  rownames(output_txt) <- c(colnames(res$X), colnames(res$Oort), c("r", "N", "r_mean", "ePx"))
  
  write_lines(stargazer(output_txt, type = "text"), paste0(s_,".txt"), append = FALSE)
}

tgas_calc_OM_seq <- function(tgas_ = tgas, src_ = "TGAS", start = 1, step = 0.1, q = 2, px_type = "ANGLE", distance = NULL, save = NULL, ...)
{
  if (!is.null(distance))
    q <- nrow(distance)
  res <- matrix(0, q, 11)
  err <- matrix(0, q, 11)
  par <- matrix(0, q, 5)
  sol <- matrix(0, q, 3)
  colnames(sol) <- c("X", "Y", "Z")
  oort <- matrix(0, q, 4)
  oort_err <- matrix(0, q, 4)
  
  if (!is.null(distance))
  {
    ds <- distance[q,2]%/%1 + 1 
  } else 
  {
    ds <- (start+step*(q+1))%/%1 + 1 
  }
  
  solution <- list();
  
  for (i in 1:q)
  {
    if (!is.null(distance))
    {
      par[i, 1] <- distance[i,1]
      par[i, 2] <- distance[i,2]
    } else {
      par[i,1] <- start+step*i
      par[i,2] <- start+step*(i+1)
    }

    if (px_type!="ANGLE")
    {
      px_ = c(1/par[i,2], 1/par[i,1])
      colnames(par) <- c("r_min","r_max","number of stars","r_mean", "px_err")
    } else
    {
      px_ = c(par[i,1], par[i,2])
      colnames(par) <- c("px_min","px_max","number of stars","px_mean", "px_err")
    }

    cat("filtering...", "\n")
    tgas_sample <- filter_tgs_px(tgas_, px = px_, ...);
    if (nrow(tgas_sample)==0)
    {
      par[i,3] <- 0
      par[i,4] <- 0
      par[i,5] <- 0
      res[i, ] <- c(rep(0, 11))
      err[i, ] <- c(rep(0, 11))
      sol[i, ] <- c(rep(0, 3))
      next()
    }
    cat("EQ to GAL converting...", "\n")
    tgas_sample <- tgas_calc_gpm(tgas_sample)

    if(!is.null(save))
    {
      cat("Saving...", "\n")
      tgas_sample <- mutate(tgas_sample, R = 1/gPx)
      s <- paste0(save, "_", src_, "_", par[i,1], "-",par[i,2])
      #write_csv(select(tgas_sample, RA, DE, gpmRA, gpmDE, tyc_pmRA, tyc_pmDE, gl, gb, R, gpm_l, gpm_b, tpm_l, tpm_b, apasm_b, apasm_v), paste0(s, "_sample.csv"), col_names = TRUE)
      x = as.matrix(select(tgas_sample, RA, DE, gpmRA, gpmDE, tyc_pmRA, tyc_pmDE, gl, gb, R, gpm_l, gpm_b, tpm_l, tpm_b, apasm_b, apasm_v, parallax_error, TYC))
      write.fwf(x, file = paste0(s, "_sample.txt"), colnames = TRUE, sep = "   ")
      HRDiagram(tgas_sample, save = s)
      DrawGalaxyPlane(tgas_sample, plane = "XY", save = s, dscale = ds)
      DrawGalaxyPlane(tgas_sample, plane = "XZ", save = s, dscale = ds)
      DrawGalaxyPlane(tgas_sample, plane = "YZ", save = s, dscale = ds)
    }

    cat("Equation of conditions forming...", "\n")
    stars <- tgas_get_stars(tgas_sample, src_)

    par[i,3] <- nrow(stars)
    par[i,4] <- mean(stars[,3])
    par[i,5] <- mean(tgas_sample$parallax_error)
    cat("Solution calculation...", "\n")
    
    res_tgas <- Calc_OM_Model(stars, use_vr = FALSE, mode = 2, scaling = 0, ef = 11)
    
    res[i, ] <- res_tgas$X
    err[i, ] <- res_tgas$s_X
    sol[i, ] <- res_tgas$X[1:3]/par[i,4]
    oort[i, ] <- res_tgas$Oort
    oort_err[i, ] <- res_tgas$s_Oort
    print(par[i,])
    print(res_tgas$X)
    print(res_tgas$s_X)
    solution[[i]] <- res_tgas
  }
  colnames(res) <- names(res_tgas$X)
  colnames(err) <- names(res_tgas$s_X)
  colnames(oort) <- names(res_tgas$Oort)
  colnames(oort_err) <- names(res_tgas$s_Oort)

  res <- list(X = res, S_X = err, Sol = sol, Parameters = par, Oort = oort, s_Oort = oort_err)
  res$SolutionR <- solution;
  
  return(res)
}

tgas_calc_OM_cond <- function(tgas_ = tgas, lclass = 3, population = "ALL", src = "TGAS")
{

  #APASS photometry
  
  conditions <- list();
  
  conditions$Population <- population;
  if (population == "DISK")
  {
    conditions$Z <- c(0, 0.5)
  } else if (population == "GALO")
  {
    conditions$Z <- c(0.25, Inf)
  } else 
  {
    conditions$Z <- c(0, Inf)
  }

  conditions$LClass <- lclass;
  if(lclass == 1)
  {
    tgas_ <- tgas_[tgas_$LClass_apass == 1,]
    conditions$BV <- c(-Inf, Inf)
    conditions$MG <- c(-Inf, -2)
    conditions$e_Px <- Inf
    
    if (population == "DISK")
    {
      distance_ <- matrix(0, nrow = 4, ncol = 2)
      distance_[,1] <- c(0.0, 2.5, 5.0, 7.5)
      distance_[,2] <- c(2.5, 5.0, 7.5, 10.0)
      if (!dir.exists("SG_Disk")) 
        dir.create("SG_Disk")
      SaveTo <- "SG_Disk/SG"
    } else if (population == "GALO")
    {
      distance_ <- matrix(0, nrow = 7, ncol = 2)
      distance_[,1] <- c(0.0, 2.5, 5.0, 7.5, 10.0, 15.0, 20.0)
      distance_[,2] <- c(2.5, 5.0, 7.5, 10.0, 15.0, 20.0, 50.0)
      if (!dir.exists("SG_Galo")) 
        dir.create("SG_Galo")
      SaveTo <- "SG_Galo/SG"
    } else 
    {
      distance_ <- matrix(0, nrow = 7, ncol = 2)
      distance_[,1] <- c(0.0, 2.5, 5.0, 7.5, 10.0, 15.0, 20.0)
      distance_[,2] <- c(2.5, 5.0, 7.5, 10.0, 15.0, 20.0, 50.0)
      if (!dir.exists("SG")) 
        dir.create("SG")
      SaveTo <- "SG/SG"
    }

  } else if(lclass == 3)
  {
    tgas_ <- tgas_[tgas_$LClass_apass == 3,]
    #conditions$BV <- c(0.75, 1.75)
    #conditions$MG <- c(-1, 2)

    conditions$BV <- c(0.8, 2.5)
    conditions$MG <- c(-1.5, 2.5)
    conditions$e_Px <- Inf #1.5
    
    if (population == "DISK")
    {
      distance_ <- matrix(0, nrow = 14, ncol = 2)
      distance_[,1] <- c(0.2, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.5, 1.75, 2.0)
      distance_[,2] <- c(0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.5, 1.75, 2.0, 3)
      if (!dir.exists("RG_Disk")) 
        dir.create("RG_Disk")
      SaveTo <- "RG_Disk/RG"
    } else if (population == "GALO")
    {
      distance_ <- matrix(0, nrow = 16, ncol = 2)
      distance_[,1] <- c(0.5, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.6, 1.8, 2.0, 2.25, 2.5, 2.75, 3)
      distance_[,2] <- c(0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.6, 1.8, 2.0, 2.25, 2.5, 2.75, 3, 4)
      if (!dir.exists("RG_Galo")) 
        dir.create("RG_Galo")
      SaveTo <- "RG_Galo/RG"
    }else 
    {
      distance_ <- matrix(0, nrow = 20, ncol = 2)
      distance_[,1] <- c(0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.6, 1.8, 2.0, 2.25, 2.5, 2.75, 3)
      distance_[,2] <- c(0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.6, 1.8, 2.0, 2.25, 2.5, 2.75, 3, 4)
      if (!dir.exists("RG")) 
        dir.create("RG")
      SaveTo <- "RG/RG"
    }

  } else if(lclass == 5)
  {
    tgas_ <- tgas_[tgas_$LClass_apass == 5,]
    
    conditions$BV <- c(-Inf, Inf)
    conditions$MG <- c(-Inf, Inf)
    conditions$e_Px <- 0.5
    if (!dir.exists("MS")) 
      dir.create("MS")
    SaveTo <- "MS/MS"

    distance_ <- matrix(0, nrow = 14, ncol = 2)
    distance_[,1] <- c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.5)
    distance_[,2] <- c(0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.5, 2.0)
  } else 
  {
    SaveTo <- "NC/NC"
  }
  
  con <- file(paste0(SaveTo,"_description.txt"), "w")
  cat("B-V = ", conditions$BV, "\n", file=con)
  cat("M = ", conditions$MG, "\n", file=con)
  cat("ePX = ", conditions$e_Px, "\n", file=con)
  cat("z = ", conditions$Z, "\n", file=con)
  cat("LClass = ", conditions$LClass , "\n", file=con)
  cat("Distance = ", distance_ , "\n", file=con)
  cat("Population = ", conditions$Population , "\n", file=con)
  
  close(con)

  solution  <- tgas_calc_OM_seq(tgas_, src_ = src, start = start, step = step, q = q, z_lim = conditions$Z, e_px = conditions$e_Px, bv = conditions$BV, Mg = conditions$MG, px_type = "DIST", distance = distance_, save = SaveTo)
  solution$Conditions <- conditions
  
  
  tgas_export_solution(solution, paste0(SaveTo, "_", src))
  tgas_draw_solution(solution, save = paste0(SaveTo, "_", src,"_", ph, "_"))
    

  return(solution)
  
}
  
tgas_draw_solution <- function (solution, save)
{
  draw_OM(solution, title = paste("Ogorodnikov-Miln Model, TGAS proper motions. Photometry:", ph))
  ggsave(paste0(save, "OM-R", ".png"), width = 10, height = 10)
  
  draw_Oort(solution, title = paste("Oort-Lindblad Model, TGAS proper motions. Photometry:", ph))
  ggsave(paste0(save, "OL-R", ".png"), width = 10, height = 10)
  
  draw_OM_Solar(solution, title = paste("Ogorodnikov-Miln Model, TGAS proper motions, Solar motion. Photometr:", ph))
  ggsave(paste0(save, "Solar-R", ".png"), width = 10, height = 10)
}


tgas_make_all_solutions <- function()
{
  solutions <- list();
 
  solutions$SG_ALL <-  tgas_calc_OM_cond(tgas, lclass = 1)
  solutions$SG_Disk <- tgas_calc_OM_cond(tgas, lclass = 1, population = "DISK")
  solutions$SG_Galo <- tgas_calc_OM_cond(tgas, lclass = 1, population = "GALO")
  solutions$RG_All <-  tgas_calc_OM_cond(tgas, lclass = 3, population = "ALL")
  solutions$RG_Disk <- tgas_calc_OM_cond(tgas, lclass = 3, population = "DISK")
  solutions$RG_Galo <- tgas_calc_OM_cond(tgas, lclass = 3, population = "GALO")
  solutions$MS_All <-  tgas_calc_OM_cond(tgas, lclass = 5, population = "DISK")
  
  return(solutions)
}

  #res_tgas_s  <- tgas_calc_OM_seq(tgas_, src_ = src, start = start, step = step, q = q, z_lim = Z, e_px = e_Px, bv = BV, Mg = MG, px_type = "DIST", distance = distance_, save = SaveTo)

  #draw_OM(res_tgas_s, title = paste("Ogorodnikov-Miln Model, TYCHO proper motions. Photometry:", ph))
  #ggsave(paste0(SaveTo, "OM-Px-TYCHO_02-",ph,".png"), width = 10, height = 10)
  
  #draw_Oort(res_tgas_s, title = paste("Oort-Lindblad Model, TYCHO proper motions. Photometry:", ph))
  #ggsave(paste0(SaveTo, "OL-Px-TYCHO_02-",ph,".png"), width = 10, height = 10)

  #draw_OM_Solar(res_tgas_s, paste("Ogorodnikov-Miln Model, TYCHO proper motions, Solar motion. Photometry:",ph))
  #ggsave(paste0(SaveTo, "Solar-Px_TYCHO_02-",ph,".png"), width = 10, height = 10)

  #res2 <- res_tgas
  #res2$X <- res_tgas$X - res_tgas_s$X

  #draw_OM_diff(res2, title = paste("Ogorodnikov-Miln Model, difference TGAS-TYCHO. Photometry:", ph))
  #ggsave(paste0(SaveTo, "OM-Px-TGAS-TYCHO_02-",ph,".png"), width = 10, height = 10)

  #draw_OM_Solar_diff(res2, title = paste("Ogorodnikov-Miln Model, difference TGAS-TYCHO, Solar motions. Photometry:",ph))
  #ggsave(paste0(SaveTo, "Solar-PX_TGAS-TYCHO_02-",ph,".png"), width = 10, height = 10)



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


make_ucac4_df <- function()
{
  res <-data.frame(ra = integer(0), spd = integer(0), u_magm = integer(0), u_maga = integer(0), smag = integer(0),
                   objt=integer(0), cdf=integer(0), sigra=integer(0), sigdc=integer(0),
                   na1=integer(0), nu1=integer(0),cu1=integer(0),
                   cepra=integer(0), cepdc=integer(0), pmrac=integer(0), pmdc=integer(0), sigpmra=integer(0), sigpmdc=integer(0),
                   pts_key=integer(0), j_m=integer(0), h_m=integer(0), k_m=integer(0),
                   icqflag_1=integer(0),icqflag_2=integer(0), icqflag_3=integer(0),
                   e2mpho_1=integer(0), e2mpho_2=integer(0),  e2mpho_3=integer(0),
                   apasm_b=integer(0), apasm_v=integer(0), apasm_g=integer(0), apasm_r=integer(0), apasm_i=integer(0),
                   eapasm_b=integer(0), eapasm_v=integer(0), eapasm_g=integer(0), eapasm_r=integer(0), eapasm_i=integer(0),
                   gcflg=integer(0), icf=integer(0), leda=integer(0), x2m=integer(0), uc4_id=integer(0), zn2=integer(0), rn2=integer(0),
                   isHIP = character(0))
  return(res)
}

read_ucac4_rec <- function(con, num, id = NULL, index = NULL)
{
  res <- make_ucac4_df()

  #t1 <- system.time(
  for (i in 1:num)
  {
    if (!is.null(index))
    {
      seek(con, where = (index[i]-1)*78, rw = "read")
    }

    res[i,1:2] <- readBin(con, "integer", n = 2L, size = 4)
    res[i,3:4] <- readBin(con, "integer", n = 2L, size = 2)
    res[i,5:12] <- readBin(con, "integer", n = 8L, size = 1)
    res[i,13:16] <- readBin(con, "integer", n = 4L, size = 2)
    res[i,17:18] <- readBin(con, "integer", n = 2L, size = 1)
    res[i,19] <- readBin(con, "integer", n = 1L, size = 4)
    res[i,20:22] <- readBin(con, "integer", n = 3L, size = 2)
    res[i,23:28] <- readBin(con, "integer", n = 6L, size = 1)
    res[i,29:33] <- readBin(con, "integer", n = 5L, size = 2)
    res[i,34:38] <- readBin(con, "integer", n = 5L, size = 1)
    res[i,39] <- readBin(con, "integer", n = 1L, size = 1)
    res[i,40] <- readBin(con, "integer", n = 1L, size = 4)
    res[i,41:42] <- readBin(con, "integer", n = 2L, size = 1)
    res[i,43] <- readBin(con, "integer", n = 1L, size = 4)
    res[i,44] <- readBin(con, "integer", n = 1L, size = 2)
    res[i,45] <- readBin(con, "integer", n = 1L, size = 4)
  }
  #)

  #cat(t1,"\n")
  #t2 <- system.time({
  if( !is.null(id))
  {
    res <- res %>% filter(uc4_id %in% id)
  }

  res <- mutate(res, ra = ra / 3600000, spd = (spd/3600000) - 90, u_magm = u_magm/1000, u_maga = u_magm/1000,
                sigra = sigra +128, sigdc = sigdc +128, sigpmra = sigpmra + 128, sigpmdc = sigpmdc + 128,
                j_m = j_m/1000, h_m = h_m/1000, k_m = k_m/1000, cepra = cepra /100, cepdc = cepdc /100,
                eapasm_b = eapasm_b/100, eapasm_v = eapasm_v /100, eapasm_g = eapasm_g/100, eapasm_r = eapasm_r/100, eapasm_i = eapasm_i/100)
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
  res <- mutate(res, isHIP = substr(as.character(icf), 1, 1))
  #})
  #cat(t2,"\n")
  return(res)
}

read_ucac4 <- function(path, start = 1, n = Inf, is_tyc_only = TRUE)
{

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


    if (is_tyc_only)
    {

      data <- data %>% left_join(u4xt_index[ , names(u4xt_index) %in% c("TYC", "uc4_id")], by = "uc4_id")

      data <- data[!(is.na(data$TYC)),]
      cat("filtered (TYC records):", nrow(data), "\n")
    }

    ucac4_data <- rbind(ucac4_data, data)

    cat("Total records in catalogue:", nrow(ucac4_data), "\n")

  }


  return(ucac4_data)
}

  read_ucac4_default <- function(start = 1, n = Inf, is_tyc_only = TRUE)
  {
    ucac4_path <- "X:/Data/Catalogues/UCAC4/"
    ucac4_data <- read_ucac4(ucac4_path, start, n, is_tyc_only)
    return (ucac4_data)
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
   tgas_data <- tgas_data %>% left_join(hip_data[ , names(hip_data) %in% c("HIP", "Px", "e_Px")], by = "HIP")

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
  index <-(tgas_$gPx>0)&(!is.na(tgas_$apasm_b))&(!is.na(tgas_$apasm_v))
  tgas_$B_V[index] <- (tgas_$apasm_b[index] - tgas_$apasm_v[index])
  # M = m + 5 + 5 lg(px).
  tgas_$M[index] <- tgas_$apasm_v[index] + 5 + 5*log10(tgas_$gPx[index]/1000)
  tgas_$Mag[index] <- tgas_$apasm_v[index]
  return(tgas_)
}

tgas_get_APASS <- function(tgas_)
{
  tgas_ <- tgas_apply_APASS(tgas_)
  tgas_ <- tgas_[(!is.na(tgas_$B_V)) & (!is.na(tgas_$M)),]
  
  return (tgas_)
}

tgas_calc_LClass <- function(tgas_)
{
  tgas_a <- tgas_get_APASS(tgas_)
  tgas_a <- mutate(tgas_a, LClass_apass = NA)
  tgas_a$LClass_apass[is_main_sequence(tgas_a$B_V, tgas_a$M)] <- 5
  tgas_a$LClass_apass[tgas_a$M<(-2)] <- 1
  tgas_a$LClass_apass[(tgas_a$M>-1.5)&(tgas_a$M<2.5)&(tgas_a$B_V>0.8)&(tgas_a$B_V<2.5)] <- 3
  
  tgas_ <- tgas_ %>% left_join(tgas_a[ , names(tgas_a) %in% c("source_id", "LClass_apass")], by = "source_id")
  
  return(tgas_)
}


min_M <- function(bv)
{
  m <- numeric(length(bv))
  
  i1 <- bv < 0.8
  m[i1] <- 5.7 * bv[i1] - 1.22222;
  
  i2 <- (bv>=0.8)&(bv<1.0)
  m[i2] <- (13 * bv[i2] - 7)
  
  i3 <- (bv>=1.0)&(bv<1.5)
  m[i3] <- (4 * bv[i3] + 2)
  
  i4 <- (bv>=1.5)
  m[i4] <- (10 * bv[i4] - 7)
  
  return(m)
}

max_M <- function(bv)
{
  m <- numeric(length(bv))
  
  i1 <- bv < 0.5
  m[i1] <- 6.3333333 * bv[i1] + 2.4333333;
  
  i2 <- (bv>=0.5)&(bv<1.3)
  m[i2] <- (4.25 * bv[i2] + 3.475)
  
  i3 <- (bv>=1.3)
  m[i3] <- (20 * bv[i3] - 17)
  
  return (m)
}


is_main_sequence <- function(bv, M)
{
  return( (M<max_M(bv)) & (M>min_M(bv)))
}

# -----------------------------------------------------------------------------
#                                    Diagrams
# -----------------------------------------------------------------------------


HRDiagram <- function(data, photometric = "APASS", title = "Hertzsprung-Russell", save = NULL)
{
  if (photometric == "TYCHO")
  {
    data <- data %>% mutate(M = NA, B_V = NA)
    s <- (data$gPx>0) & (!is.na(data$tyc_m))
    data$M[s] <- data$tyc_m[s] + 5 + 5*log10(data$gPx[s]/1000)
    hrdata <- data.frame(cbind( M = data$M[!is.na(data$M)], B_V = data$tyc_bv[!is.na(data$M)]))
  } else if (photometric == "APASS")
  {
    data <- data %>% mutate(M = NA)
    index <-(data$gPx>0)&(!is.na(data$apasm_v))
    data$M[index] <- data$apasm_v[index] + 5 + 5*log10(data$gPx[index]/1000)
    hrdata <- data.frame(cbind( M = data$M[index], B_V = (data$apasm_b[index]-data$apasm_v[index]), LC = data$LClass[!is.na(data$M)]))
  }
  else
  {
    data <- data %>% mutate(M = NA)
    data$M[data$gPx>0] <- data$Gm_mag[data$gPx>0] + 5 + 5*log10(data$gPx[data$gPx>0]/1000)
    hrdata <- data.frame(cbind( M = data$M[!is.na(data$M)], B_V = data$B_V[!is.na(data$M)], LC = data$LClass[!is.na(data$M)]))
  }

  if(nrow(hrdata)<1000)
  {
    alpha_ <- 1
    size_ <- 1
  } else if(nrow(hrdata)<10000)
  {
    alpha_ <- 0.5
    size_ <- 1
  } else if (nrow(hrdata)<100000)
  {
    alpha_ <- 0.25
    size_ <- 0.5
  } else if (nrow(hrdata)<1000000)
  {
    alpha_ <- 0.1
    size_ <- 0.1
  } else
  {
    alpha_ <- 0.05
    size_ <- 0.1
  }

  #alpha_ <- 1

  #g <- ggplot() + geom_point(data=hrdata, aes(x = hrdata$B_V, y = hrdata$M), alpha_ = 0.05, na.rm = TRUE, size_ = 0.1) + scale_y_reverse()

  g <- ggplot() +
          #scale_colour_manual("L Class",  breaks = colnames(c(1,2,3,4,5)),
           #           values = c("blue", "brown", "red", "yellow", "green")) +
          #scale_color_continuous(high = "red", low = "blue") + 
          geom_point(data=hrdata, aes(x = hrdata$B_V, y = hrdata$M),  alpha = alpha_, na.rm = TRUE, size = size_, shape = ".") +  
          scale_y_reverse(breaks=seq(10,-10,by=-1), minor_breaks=seq(10,-10,by=-0.5), limits = c(10,-10)) +
          scale_x_continuous(breaks=seq(-1,3,by=0.25), minor_breaks=seq(-1,3,by=0.125), limits = c(-1,3)) +
          xlab("B-V") + ylab("M") + ggtitle(title)

  ms_top_limit <- matrix(0, nrow = 5, ncol = 2)
  ms_top_limit[,1] <- c(-0.1, 0.8, 1.0, 1.5, 1.7)
  ms_top_limit[,2] <- c(-1.8, 3.4, 6.0, 8.0, 10.0 )
  ms_bottom_limit <- matrix(0, nrow = 4, ncol = 2)
  ms_bottom_limit[,1] <- c(-0.1, 0.5, 1.3, 1.35)
  ms_bottom_limit[,2] <- c(1.8, 5.6, 9.0, 10.0)

  g <- g + geom_line(aes(x = ms_top_limit[,1], y = ms_top_limit[,2])) +
           geom_line(aes(x = ms_bottom_limit[,1], y = ms_bottom_limit[,2]));

  a <- seq(from = 0, to = 2, by = 0.1)
  g <- g + geom_line(aes(x = a, y = max_M(a))) +
           geom_line(aes(x = a, y = min_M(a)));

  rg_top_limit <- matrix(0, nrow = 4, ncol = 2)
  rg_top_limit[,1] <- c(0.8, 0.8, 2.5, 2.5)
  rg_top_limit[,2] <- c(2.5, -1.5, -1.5, 2.5)
  rg_bottom_limit <- matrix(0, nrow = 2, ncol = 2)
  rg_bottom_limit[,1] <- c(0.8, 2.5)
  rg_bottom_limit[,2] <- c(2.5, 2.5)
  
  #g <- g + geom_path()

  g <- g + geom_line(aes(x = rg_top_limit[,1], y = rg_top_limit[,2])) +
           geom_line(aes(x = rg_bottom_limit[,1], y = rg_bottom_limit[,2]))


  if(!is.null(save))
  {
    ggsave(paste0(save, "-HR.png"), width = 10, height = 10)
  }

  return(g)
}

DrawGalaxyPlane <- function(data, plane = "XZ", title = "Star distribution", save = NULL, dscale = 5)
{
  
  #data <- mutate (data, z = (1/gPx)*sin(b), x = (1/gPx)*cos(b)*cos(l), y = (1/gPx)*cos(b)*sin(l))
  data <- CalcGalXYZ(data)
  
  g <- ggplot()
  
  if ( (plane == "XY") | (plane == "YX"))
  {
    hrdata <- data.frame(cbind(data$x, data$y))
    g <- g + xlab("X") + ylab("Y")
  } else if ((plane == "XZ") | (plane == "ZX"))
  {
    hrdata <- data.frame(cbind(data$x, data$z))
    g <- g + xlab("X") + ylab("Z")
  }
  else
  {
    hrdata <- data.frame(cbind(data$y, data$z))
    g <- g + xlab("Y") + ylab("Z")
  }
  
  #min_x <- (max(min(hrdata[,1]),-9)%/%1)*1 - 1
  #max_x <- (min(max(hrdata[,1]), 9)%/%1)*1 + 1
  #min_y <- (max(min(hrdata[,2]), -9)%/%1)*1 - 1
  #max_y <- (min(max(hrdata[,2]), 9)%/%1)*1 + 1
  
  min_x <- -dscale
  max_x <- dscale
  min_y <- -dscale
  max_y <- dscale
  
  if(nrow(hrdata)<1000)
  {
    alpha_ <- 1
    size_ <- 1
  } else if(nrow(hrdata)<10000)
  {
    alpha_ <- 0.5
    size_ <- 1
  } else if (nrow(hrdata)<100000)
  {
    alpha_ <- 0.25
    size_ <- 0.5
  } else if (nrow(hrdata)<1000000)
  {
    alpha_ <- 0.1
    size_ <- 0.1
  } else
  {
    alpha_ <- 0.05
    size_ <- 0.1
  }
  
  g <- g +
    geom_point(data=hrdata, aes(x = hrdata[,1], y = hrdata[,2]),  alpha = alpha_, na.rm = TRUE, size = size_, shape = ".") + 
    scale_y_continuous(breaks=seq(min_y,max_y,by=1), minor_breaks=seq(min_y,max_y,by=0.5), limits = c(min_y,max_y)) +
    scale_x_continuous(breaks=seq(min_x,max_x,by=1), minor_breaks=seq(min_x,max_x,by=0.5), limits = c(min_x,max_x)) +
    ggtitle(title)
  
 
  if(!is.null(save))
  {
    ggsave(paste0(save, plane, ".png"), width = 10, height = 10)
  }
  
  return(g)
}


draw_HR_TGAS_T2Sp <- function()
{
  HRDiagram(tgas_, title = "Hertzsprung?Russell")
  ggsave("HR-TGAS-T2sp.png", width = 10, height = 10)
  HRDiagram(filter(tgas_, LClass == 1), title = "Hertzsprung?Russell L1 class")
  ggsave("HR-TGAS-T2sp-L1.png", width = 10, height = 10)
  HRDiagram(filter(tgas_, LClass == 2), title = "Hertzsprung?Russell L2 class")
  ggsave("HR-TGAS-T2sp-L2.png", width = 10, height = 10)
  HRDiagram(filter(tgas_, LClass == 3), title = "Hertzsprung?Russell L3 class")
  ggsave("HR-TGAS-T2sp-L3.png", width = 10, height = 10)
  HRDiagram(filter(tgas_, LClass == 4), title = "Hertzsprung?Russell L4 class")
  ggsave("HR-TGAS-T2sp-L4.png", width = 10, height = 10)
  HRDiagram(filter(tgas_, LClass == 5), title = "Hertzsprung?Russell L5 class")
  ggsave("HR-TGAS-T2sp-L5.png", width = 10, height = 10)
}

draw_OM <- function(res, title = "Ogorodnikov-Miln Model")
{
  min_x <- (min(res$Parameters[,4])%/%0.2)*0.2
  max_x <- (max(res$Parameters[,4])%/%0.2)*0.2 + 0.2
  stepx <- 0.2
  #stepx <- 10
  
  #min_y <- (min(res$X[,4:11])%/%1) - 1
  #max_y <- (max(res$X[,4:11])%/%1) + 1
  min_y <- -17
  max_y <- 18
  
  
  g <- ggplot()
  g <- g +
    #scale_y_continuous(breaks=seq(-18,18,by=3), minor_breaks=seq(-18,18,by=1), limits = c(-17,17)) +
    scale_y_continuous(breaks=seq(min_y,max_y,by=1), minor_breaks=seq(min_y,max_y,by=0.5), limits = c(min_y,max_y)) +
    #scale_x_continuous(breaks=seq(0.25,4,by=0.25), minor_breaks=seq(0.25,4,by=0.125), limits = c(0.25,4))
    scale_x_continuous(breaks=seq(min_x,max_x,by=stepx), minor_breaks=seq(min_x,max_x,by=stepx/2), limits = c(min_x,max_x))
    #scale_x_continuous(breaks=seq(1,4.5,by=0.5), minor_breaks=seq(1,4.5,by=0.25), limits = c(1,4.5))

  g <- g + xlab("<px>, kpc") + ylab("O-M, km/s/kpc") +ggtitle(title) +
    scale_colour_manual("Parameters",  breaks = colnames(res$X)[4:11],
                        values = c("#000000", "#E69F00", "#56B4E9", "#009E73", "green", "#0072B2", "#D55E00", "#CC79A7")) +
    #geom_line(aes(x = res$Parameters[,4], y = res$X[,1], colour = colnames(res$X)[1]), size = 1) +
    #geom_line(aes(x = res$Parameters[,4], y = res$X[,2], colour = colnames(res$X)[2]), size = 1) +
    #geom_line(aes(x = res$Parameters[,4], y = res$X[,3], colour = colnames(res$X)[3])) +
    geom_line(aes(x = res$Parameters[,4], y = res$X[,4], colour = colnames(res$X)[4]), size = 1) +
    geom_errorbar(aes(x = res$Parameters[,4], ymin = res$X[,4] - res$S_X[,4], ymax = res$X[,4] + res$S_X[,4], colour = colnames(res$X)[4])) +
    geom_line(aes(x = res$Parameters[,4], y = res$X[,5], colour = colnames(res$X)[5]), size = 1) +
    geom_errorbar(aes(x = res$Parameters[,4], ymin = res$X[,5] - res$S_X[,5], ymax = res$X[,5] + res$S_X[,5], colour = colnames(res$X)[5])) +
    geom_line(aes(x = res$Parameters[,4], y = res$X[,6], colour = colnames(res$X)[6]), size = 1) +
    geom_errorbar(aes(x = res$Parameters[,4], ymin = res$X[,6] - res$S_X[,6], ymax = res$X[,6] + res$S_X[,6], colour = colnames(res$X)[6])) +
    geom_line(aes(x = res$Parameters[,4], y = res$X[,7], colour = colnames(res$X)[7]), size = 1) +
    geom_errorbar(aes(x = res$Parameters[,4], ymin = res$X[,7] - res$S_X[,7], ymax = res$X[,7] + res$S_X[,7], colour = colnames(res$X)[7])) +
    geom_line(aes(x = res$Parameters[,4], y = res$X[,8], colour = colnames(res$X)[8]), size = 1) +
    geom_errorbar(aes(x = res$Parameters[,4], ymin = res$X[,8] - res$S_X[,8], ymax = res$X[,8] + res$S_X[,8], colour = colnames(res$X)[8])) +
    geom_line(aes(x = res$Parameters[,4], y = res$X[,9], colour = colnames(res$X)[9]), size = 1) +
    geom_errorbar(aes(x = res$Parameters[,4], ymin = res$X[,9] - res$S_X[,9], ymax = res$X[,9] + res$S_X[,9], colour = colnames(res$X)[9])) +
    geom_line(aes(x = res$Parameters[,4], y = res$X[,10], colour = colnames(res$X)[10]), size = 1) +
    geom_errorbar(aes(x = res$Parameters[,4], ymin = res$X[,10] - res$S_X[,10], ymax = res$X[,10] + res$S_X[,10], colour = colnames(res$X)[10])) +
    geom_line(aes(x = res$Parameters[,4], y = res$X[,11], colour = colnames(res$X)[11]), size = 1) +
    geom_errorbar(aes(x = res$Parameters[,4], ymin = res$X[,11] - res$S_X[,11], ymax = res$X[,11] + res$S_X[,11], colour = colnames(res$X)[11])) +
    geom_point(aes(x = res$Parameters[,4], y = rep(0,nrow(res$Parameters)))) #+  
    #geom_smooth(aes(x = res$Parameters[,4], y = res$X[,6], colour = colnames(res$X)[6])) +
    #geom_smooth(aes(x = res$Parameters[,4], y = res$X[,9], colour = colnames(res$X)[9])) +
    #geom_smooth(aes(x = res$Parameters[,4], y = res$X[,10], colour = colnames(res$X)[10])) 
  return(g)
}


draw_OM_diff <- function(res, title = "Ogorodnikov-Miln Model")
{
  min_x <- (min(res$Parameters[,4])%/%0.2)*0.2
  max_x <- (max(res$Parameters[,4])%/%0.2)*0.2 + 0.2

  min_y <- (min(res$X[,4:11])%/%1) - 1
  max_y <- (max(res$X[,4:11])%/%1) + 1

  g <- ggplot()

  g <- g +
    #scale_y_continuous(breaks=seq(-6,7,by=0.5), minor_breaks=seq(-6,7,by=0.25), limits = c(-6,7)) +
    scale_y_continuous(breaks=seq(min_y,max_y,by=1), minor_breaks=seq(min_y,max_y,by=0.5), limits = c(min_y,max_y)) +
    #scale_x_continuous(breaks=seq(0,4,by=0.25), minor_breaks=seq(0,4,by=0.125), limits = c(0.2,4))
    #scale_x_continuous(breaks=seq(0,2.5,by=0.25), minor_breaks=seq(0,2.5,by=0.125), limits = c(0.2,2.5))
    scale_x_continuous(breaks=seq(min_x,max_x,by=0.2), minor_breaks=seq(min_x,max_x,by=0.1), limits = c(min_x,max_x))


  g <- g + xlab("<px>, kpc") + ylab("O-M, km/s/kpc") +ggtitle(title) +
    scale_colour_manual("Parameters",  breaks = colnames(res$X)[4:11],
                        values = c("#000000", "#E69F00", "#56B4E9", "#009E73", "green", "#0072B2", "#D55E00", "#CC79A7")) +
    geom_line(aes(x = res$Parameters[,4], y = res$X[,4], colour = colnames(res$X)[4]), size = 1) +
    geom_line(aes(x = res$Parameters[,4], y = res$X[,5], colour = colnames(res$X)[5]), size = 1) +
    geom_line(aes(x = res$Parameters[,4], y = res$X[,6], colour = colnames(res$X)[6]), size = 1) +
    geom_line(aes(x = res$Parameters[,4], y = res$X[,7], colour = colnames(res$X)[7]), size = 1) +
    geom_line(aes(x = res$Parameters[,4], y = res$X[,8], colour = colnames(res$X)[8]), size = 1) +
    geom_line(aes(x = res$Parameters[,4], y = res$X[,9], colour = colnames(res$X)[9]), size = 1) +
    geom_line(aes(x = res$Parameters[,4], y = res$X[,10], colour = colnames(res$X)[10]), size = 1) +
    geom_line(aes(x = res$Parameters[,4], y = res$X[,11], colour = colnames(res$X)[11]), size = 1) +
    geom_point(aes(x = res$Parameters[,4], y = rep(0,nrow(res$Parameters))))
  return(g)
}

draw_OM_Solar <- function(res, title = "Solar motion")
{
  min_x <- (min(res$Parameters[,4])%/%0.2)*0.2
  max_x <- (max(res$Parameters[,4])%/%0.2)*0.2 + 0.2

  min_y <- (min(res$X[,1:3])%/%1) - 1
  max_y <- (max(res$X[,1:3])%/%1) + 1

  g <- ggplot() +
    #scale_y_continuous(breaks=seq(0,20,by=1), minor_breaks=seq(0,20,by=0.5), limits = c(-0.5,20)) +
    scale_y_continuous(breaks=seq(min_y,max_y,by=1), minor_breaks=seq(min_y,max_y,by=0.5), limits = c(min_y,max_y)) +
    #scale_x_continuous(breaks=seq(0,4,by=0.25), minor_breaks=seq(0,4,by=0.125), limits = c(0.25,4)) +
    #scale_x_continuous(breaks=seq(0,2.5,by=0.25), minor_breaks=seq(0,2.5,by=0.125), limits = c(0.25,2.5)) +
    scale_x_continuous(breaks=seq(min_x,max_x,by=0.2), minor_breaks=seq(min_x,max_x,by=0.1), limits = c(min_x,max_x))+
    xlab("<px>, kpc") + ylab("km/s") +ggtitle(title) +
    scale_colour_manual("Parameters",  breaks = colnames(res$X)[1:3],
                        values = c("green", "blue", "brown")) +
    geom_line(aes(x = res$Parameters[,4], y = res$X[,1], colour = colnames(res$X)[1]), size = 1) +
    geom_errorbar(aes(x = res$Parameters[,4], ymin = res$X[,1] - res$S_X[,1], ymax = res$X[,1] + res$S_X[,1], colour = colnames(res$X)[1])) +
    geom_line(aes(x = res$Parameters[,4], y = res$X[,2], colour = colnames(res$X)[2]), size = 1) +
    geom_errorbar(aes(x = res$Parameters[,4], ymin = res$X[,2] - res$S_X[,2], ymax = res$X[,2] + res$S_X[,2], colour = colnames(res$X)[2])) +
    geom_line(aes(x = res$Parameters[,4], y = res$X[,3], colour = colnames(res$X)[3]), size = 1) +
    geom_errorbar(aes(x = res$Parameters[,4], ymin = res$X[,3] - res$S_X[,3], ymax = res$X[,3] + res$S_X[,3], colour = colnames(res$X)[3])) +
    geom_point(aes(x = res$Parameters[,4], y = rep(0,nrow(res$Parameters))))
  return(g)
}


draw_OM_Solar_diff <- function(res, title = "Solar motion")
{
  min_x <- (min(res$Parameters[,4])%/%0.2)*0.2
  max_x <- (max(res$Parameters[,4])%/%0.2)*0.2 + 0.2

  min_y <- (min(res$X[,1:3])%/%1) - 1
  max_y <- (max(res$X[,1:3])%/%1) + 1
  g <- ggplot() +
    #scale_y_continuous(breaks=seq(-2,7,by=0.5), minor_breaks=seq(-2,7,by=0.25), limits = c(-2,7)) +
    scale_y_continuous(breaks=seq(min_y,max_y,by=1), minor_breaks=seq(min_y,max_y,by=0.5), limits = c(min_y,max_y)) +
    #scale_x_continuous(breaks=seq(0,4,by=0.5), minor_breaks=seq(0,4,by=0.1), limits = c(0.25,4)) +
    scale_x_continuous(breaks=seq(min_x,max_x,by=0.2), minor_breaks=seq(min_x,max_x,by=0.1), limits = c(min_x,max_x))+
    xlab("<px>, kpc") + ylab("km/s") +ggtitle(title) +
    scale_colour_manual("Parameters",  breaks = colnames(res$X)[1:3],
                        values = c("green", "blue", "brown")) +
    geom_line(aes(x = res$Parameters[,4], y = res$X[,1], colour = colnames(res$X)[1]), size = 1) +
    geom_line(aes(x = res$Parameters[,4], y = res$X[,2], colour = colnames(res$X)[2]), size = 1) +
    geom_line(aes(x = res$Parameters[,4], y = res$X[,3], colour = colnames(res$X)[3]), size = 1) +
    geom_point(aes(x = res$Parameters[,4], y = rep(0,nrow(res$Parameters))))
  return(g)
}


draw_Oort <- function(res, title = "Oort`s parameters")
{
  min_x <- (min(res$Parameters[,4])%/%0.2)*0.2
  max_x <- (max(res$Parameters[,4])%/%0.2)*0.2 + 0.2
  
  #min_y <- (min(res$Oort)%/%1) - 1
  #max_y <- (max(res$Oort)%/%1) + 1
  min_y <- -17
  max_y <- 18
  
  
  g <- ggplot() +
    #scale_y_continuous(breaks=seq(0,20,by=1), minor_breaks=seq(0,20,by=0.5), limits = c(-0.5,20)) +
    scale_y_continuous(breaks=seq(min_y,max_y,by=1), minor_breaks=seq(min_y,max_y,by=0.5), limits = c(min_y,max_y)) +
    #scale_x_continuous(breaks=seq(0,4,by=0.25), minor_breaks=seq(0,4,by=0.125), limits = c(0.25,4)) +
    #scale_x_continuous(breaks=seq(0,2.5,by=0.25), minor_breaks=seq(0,2.5,by=0.125), limits = c(0.25,2.5)) +
    scale_x_continuous(breaks=seq(min_x,max_x,by=0.2), minor_breaks=seq(min_x,max_x,by=0.1), limits = c(min_x,max_x))+
    xlab("<r>, kpc") + ylab("km/s/kpc") +ggtitle(title) +
    scale_colour_manual("Parameters",  breaks = colnames(res$Oort)[1:4],
                        values = c("green", "blue", "brown", "black")) +
    geom_line(aes(x = res$Parameters[,4], y = res$Oort[,1], colour = colnames(res$Oort)[1]), size = 1) +
    geom_errorbar(aes(x = res$Parameters[,4], ymin = res$Oort[,1] - res$s_Oort[,1], ymax = res$Oort[,1] + res$s_Oort[,1], colour = colnames(res$Oort)[1])) +
    geom_line(aes(x = res$Parameters[,4], y = res$Oort[,2], colour = colnames(res$Oort)[2]), size = 1) +
    geom_errorbar(aes(x = res$Parameters[,4], ymin = res$Oort[,2] - res$s_Oort[,2], ymax = res$Oort[,2] + res$s_Oort[,2], colour = colnames(res$Oort)[2])) +
    geom_line(aes(x = res$Parameters[,4], y = res$Oort[,3], colour = colnames(res$Oort)[3]), size = 1) +
    geom_errorbar(aes(x = res$Parameters[,4], ymin = res$Oort[,3] - res$s_Oort[,3], ymax = res$Oort[,3] + res$s_Oort[,3], colour = colnames(res$Oort)[3])) +
    geom_line(aes(x = res$Parameters[,4], y = res$Oort[,4], colour = colnames(res$Oort)[4]), size = 1) +
    geom_errorbar(aes(x = res$Parameters[,4], ymin = res$Oort[,4] - res$s_Oort[,4], ymax = res$Oort[,4] + res$s_Oort[,4], colour = colnames(res$Oort)[4])) +
    geom_point(aes(x = res$Parameters[,4], y = rep(0,nrow(res$Parameters))))
  return(g)
}

draw_solution_Oort <- function(solution, title = "Oort`s parameters")
{
  min_x <- 0
  max_x <- 4
  step_x <- 0.5
  
  #min_y <- (min(res$Oort)%/%1) - 1
  #max_y <- (max(res$Oort)%/%1) + 1
  min_y <- -17
  max_y <- 18
  
  g <- ggplot() +
    #scale_y_continuous(breaks=seq(0,20,by=1), minor_breaks=seq(0,20,by=0.5), limits = c(-0.5,20)) +
    scale_y_continuous(breaks=seq(min_y,max_y,by=1), minor_breaks=seq(min_y,max_y,by=0.5), limits = c(min_y,max_y)) +
    #scale_x_continuous(breaks=seq(0,4,by=0.25), minor_breaks=seq(0,4,by=0.125), limits = c(0.25,4)) +
    #scale_x_continuous(breaks=seq(0,2.5,by=0.25), minor_breaks=seq(0,2.5,by=0.125), limits = c(0.25,2.5)) +
    scale_x_continuous(breaks=seq(min_x,max_x,by=step_x), minor_breaks=seq(min_x,max_x,by=step_x/2), limits = c(0.2,3.5))+
    xlab("<r>, kpc") + ylab("km/s/kpc") +ggtitle(title) +
    #scale_colour_manual("Parameters",  breaks = colnames(solution$MS_All$Oort)[1:4],
     #                   values = c("green", "blue", "brown", "black")) +
    scale_colour_manual("Parameters",  breaks = c("A Main Sequence", "A Red Giants Disk", "A Red Giants Galo", 
                                                  "B Main Sequence", "B Red Giants Disk", "B Red Giants Galo", 
                                                  "C Main Sequence", "C Red Giants Disk", "C Red Giants Galo", 
                                                  "K Main Sequence", "K Red Giants Disk", "K Red Giants Galo"),
                        values = c("green1", "green3", "green4", "blue1", "blue3", "blue4","brown1", "brown3", "brown4", "gray0","gray20", "gray40")) +
    scale_fill_manual("Parameters",  breaks = c("A Main Sequence", "A Red Giants Disk", "A Red Giants Galo", 
                                                  "B Main Sequence", "B Red Giants Disk", "B Red Giants Galo", 
                                                  "C Main Sequence", "C Red Giants Disk", "C Red Giants Galo", 
                                                  "K Main Sequence", "K Red Giants Disk", "K Red Giants Galo"),
                        values = c("green1", "green3", "green4", "blue1", "blue3", "blue4","brown1", "brown3", "brown4", "gray0","gray20", "gray40")) +    
    #geom_line(aes(x = solution$MS_All$Parameters[,4], y = solution$MS_All$Oort[,1], colour = colnames(solution$MS_All$Oort)[1]), size = 1) +
    #geom_errorbar(aes(x = solution$MS_All$Parameters[,4], ymin = solution$MS_All$Oort[,1] - solution$MS_All$s_Oort[,1], ymax = solution$MS_All$Oort[,1] + solution$MS_All$s_Oort[,1], colour = colnames(solution$MS_All$Oort)[1])) +
    geom_line(  aes(x = solution$MS_All$Parameters[,4], y = solution$MS_All$Oort[,1], colour = "A Main Sequence"), size = 1) +
    geom_point( aes(x = solution$MS_All$Parameters[,4], y = solution$MS_All$Oort[,1], fill = "A Main Sequence"), shape = 21) + 
    geom_errorbar(aes(x = solution$MS_All$Parameters[,4], ymin = solution$MS_All$Oort[,1] - solution$MS_All$s_Oort[,1], ymax = solution$MS_All$Oort[,1] + solution$MS_All$s_Oort[,1], colour = "A Main Sequence")) +
    
    geom_line( aes(x = solution$MS_All$Parameters[,4], y = solution$MS_All$Oort[,2], colour = "B Main Sequence"), size = 1) +
    geom_point(aes(x = solution$MS_All$Parameters[,4], y = solution$MS_All$Oort[,2], fill = "B Main Sequence"), shape = 21) + 
    geom_errorbar(aes(x = solution$MS_All$Parameters[,4], ymin = solution$MS_All$Oort[,2] - solution$MS_All$s_Oort[,2], ymax = solution$MS_All$Oort[,2] + solution$MS_All$s_Oort[,2], colour = "B Main Sequence")) +
    
    
    geom_line(aes(x = solution$MS_All$Parameters[,4], y = solution$MS_All$Oort[,3], colour = "C Main Sequence"), size = 1) +
    geom_errorbar(aes(x = solution$MS_All$Parameters[,4], ymin = solution$MS_All$Oort[,3] - solution$MS_All$s_Oort[,3], ymax = solution$MS_All$Oort[,3] + solution$MS_All$s_Oort[,3], colour = "C Main Sequence")) +
    geom_point(aes(x = solution$MS_All$Parameters[,4], y = solution$MS_All$Oort[,3], fill = "C Main Sequence"), shape = 21) + 
    
    geom_line(aes(x = solution$MS_All$Parameters[,4], y = solution$MS_All$Oort[,4], colour = "K Main Sequence"), size = 1) +
    geom_errorbar(aes(x = solution$MS_All$Parameters[,4], ymin = solution$MS_All$Oort[,4] - solution$MS_All$s_Oort[,4], ymax = solution$MS_All$Oort[,4] + solution$MS_All$s_Oort[,4], colour = "K Main Sequence")) +
    geom_point(aes(x = solution$MS_All$Parameters[,4], y = solution$MS_All$Oort[,4], fill = "K Main Sequence"), shape = 21) + 
    
    geom_point(aes(x = solution$MS_All$Parameters[,4], y = rep(0,nrow(solution$MS_All$Parameters))))
  
  
  g <- g + 
    #geom_line(aes(x = solution$RG_Disk$Parameters[,4], y = solution$RG_Disk$Oort[,1], colour = colnames(solution$RG_Disk$Oort)[1]), size = 1) +
    #geom_errorbar(aes(x = solution$RG_Disk$Parameters[,4], ymin = solution$RG_Disk$Oort[,1] - solution$RG_Disk$s_Oort[,1], ymax = solution$RG_Disk$Oort[,1] + solution$RG_Disk$s_Oort[,1], colour = colnames(solution$RG_Disk$Oort)[1])) +
    geom_line(aes(x = solution$RG_Disk$Parameters[,4], y = solution$RG_Disk$Oort[,1], colour = "A Red Giants Disk"), size = 1) +
    geom_point( aes(x = solution$RG_Disk$Parameters[,4], y = solution$RG_Disk$Oort[,1], fill = "A Red Giants Disk"), shape = 22) + 
    geom_errorbar(aes(x = solution$RG_Disk$Parameters[,4], ymin = solution$RG_Disk$Oort[,1] - solution$RG_Disk$s_Oort[,1], ymax = solution$RG_Disk$Oort[,1] + solution$RG_Disk$s_Oort[,1], colour = "A Red Giants Disk")) +
    
    geom_line(aes(x = solution$RG_Disk$Parameters[,4], y = solution$RG_Disk$Oort[,2], colour = "B Red Giants Disk"), size = 1) +
    geom_point( aes(x = solution$RG_Disk$Parameters[,4], y = solution$RG_Disk$Oort[,2], fill = "B Red Giants Disk"), shape = 22) + 
    geom_errorbar(aes(x = solution$RG_Disk$Parameters[,4], ymin = solution$RG_Disk$Oort[,2] - solution$RG_Disk$s_Oort[,2], ymax = solution$RG_Disk$Oort[,2] + solution$RG_Disk$s_Oort[,2], colour = "B Red Giants Disk")) +
    
    geom_line(aes(x = solution$RG_Disk$Parameters[,4], y = solution$RG_Disk$Oort[,3], colour = "C Red Giants Disk"), size = 1) +
    geom_point( aes(x = solution$RG_Disk$Parameters[,4], y = solution$RG_Disk$Oort[,3], fill = "C Red Giants Disk"), shape = 22) + 
    geom_errorbar(aes(x = solution$RG_Disk$Parameters[,4], ymin = solution$RG_Disk$Oort[,3] - solution$RG_Disk$s_Oort[,3], ymax = solution$RG_Disk$Oort[,3] + solution$RG_Disk$s_Oort[,3], colour = "C Red Giants Disk")) +
    
    geom_line(aes(x = solution$RG_Disk$Parameters[,4], y = solution$RG_Disk$Oort[,4], colour = "K Red Giants Disk"), size = 1) +
    geom_point( aes(x = solution$RG_Disk$Parameters[,4], y = solution$RG_Disk$Oort[,4], fill = "K Red Giants Disk"), shape = 22) + 
    geom_errorbar(aes(x = solution$RG_Disk$Parameters[,4], ymin = solution$RG_Disk$Oort[,4] - solution$RG_Disk$s_Oort[,4], ymax = solution$RG_Disk$Oort[,4] + solution$RG_Disk$s_Oort[,4], colour = "K Red Giants Disk")) +
    
    geom_point(aes(x = solution$RG_Disk$Parameters[,4], y = rep(0,nrow(solution$RG_Disk$Parameters))))
  
  g <- g + 
    #geom_line(aes(x = solution$RG_Galo$Parameters[,4], y = solution$RG_Galo$Oort[,1], colour = colnames(solution$RG_Galo$Oort)[1]), size = 1) +
    #geom_errorbar(aes(x = solution$RG_Galo$Parameters[,4], ymin = solution$RG_Galo$Oort[,1] - solution$RG_Galo$s_Oort[,1], ymax = solution$RG_Galo$Oort[,1] + solution$RG_Galo$s_Oort[,1], colour = colnames(solution$RG_Galo$Oort)[1])) +
    geom_line(  aes(x = solution$RG_Galo$Parameters[,4], y = solution$RG_Galo$Oort[,1], colour = "A Red Giants Galo"), size = 1) +
    geom_point( aes(x = solution$RG_Galo$Parameters[,4], y = solution$RG_Galo$Oort[,1], fill = "A Red Giants Galo"), shape = 23) + 
    geom_errorbar(aes(x = solution$RG_Galo$Parameters[,4], ymin = solution$RG_Galo$Oort[,1] - solution$RG_Galo$s_Oort[,1], ymax = solution$RG_Galo$Oort[,1] + solution$RG_Galo$s_Oort[,1], colour = "A Red Giants Galo")) +
    
    geom_line(aes(x = solution$RG_Galo$Parameters[,4], y = solution$RG_Galo$Oort[,2], colour = "B Red Giants Galo"), size = 1) +
    geom_point( aes(x = solution$RG_Galo$Parameters[,4], y = solution$RG_Galo$Oort[,2], fill = "B Red Giants Galo"), shape = 23) + 
    geom_errorbar(aes(x = solution$RG_Galo$Parameters[,4], ymin = solution$RG_Galo$Oort[,2] - solution$RG_Galo$s_Oort[,2], ymax = solution$RG_Galo$Oort[,2] + solution$RG_Galo$s_Oort[,2], colour = "B Red Giants Galo")) +
    
    geom_line(aes(x = solution$RG_Galo$Parameters[,4], y = solution$RG_Galo$Oort[,3], colour = "C Red Giants Galo"), size = 1) +
    geom_point( aes(x = solution$RG_Galo$Parameters[,4], y = solution$RG_Galo$Oort[,3], fill = "C Red Giants Galo"), shape = 23) + 
    geom_errorbar(aes(x = solution$RG_Galo$Parameters[,4], ymin = solution$RG_Galo$Oort[,3] - solution$RG_Galo$s_Oort[,3], ymax = solution$RG_Galo$Oort[,3] + solution$RG_Galo$s_Oort[,3], colour = "C Red Giants Galo")) +
    
    geom_line(aes(x = solution$RG_Galo$Parameters[,4], y = solution$RG_Galo$Oort[,4], colour = "K Red Giants Galo"), size = 1) +
    geom_point( aes(x = solution$RG_Galo$Parameters[,4], y = solution$RG_Galo$Oort[,4], fill = "K Red Giants Galo"), shape = 23) + 
    geom_errorbar(aes(x = solution$RG_Galo$Parameters[,4], ymin = solution$RG_Galo$Oort[,4] - solution$RG_Galo$s_Oort[,4], ymax = solution$RG_Galo$Oort[,4] + solution$RG_Galo$s_Oort[,4], colour = "K Red Giants Galo")) +
    
    geom_point(aes(x = solution$RG_Galo$Parameters[,4], y = rep(0,nrow(solution$RG_Galo$Parameters))))
  
  # g <- g + geom_line(aes(x = solution$SG_Disk$Parameters[,4], y = solution$SG_Disk$Oort[,1], colour = colnames(solution$SG_Disk$Oort)[1]), size = 1) +
  #   geom_errorbar(aes(x = solution$SG_Disk$Parameters[,4], ymin = solution$SG_Disk$Oort[,1] - solution$SG_Disk$s_Oort[,1], ymax = solution$SG_Disk$Oort[,1] + solution$SG_Disk$s_Oort[,1], colour = colnames(solution$SG_Disk$Oort)[1])) +
  #   
  #   geom_line(aes(x = solution$SG_Disk$Parameters[,4], y = solution$SG_Disk$Oort[,2], colour = colnames(solution$SG_Disk$Oort)[2]), size = 1) +
  #   geom_errorbar(aes(x = solution$SG_Disk$Parameters[,4], ymin = solution$SG_Disk$Oort[,2] - solution$SG_Disk$s_Oort[,2], ymax = solution$SG_Disk$Oort[,2] + solution$SG_Disk$s_Oort[,2], colour = colnames(solution$SG_Disk$Oort)[2])) +
  #   
  #   geom_line(aes(x = solution$SG_Disk$Parameters[,4], y = solution$SG_Disk$Oort[,3], colour = colnames(solution$SG_Disk$Oort)[3]), size = 1) +
  #   geom_errorbar(aes(x = solution$SG_Disk$Parameters[,4], ymin = solution$SG_Disk$Oort[,3] - solution$SG_Disk$s_Oort[,3], ymax = solution$SG_Disk$Oort[,3] + solution$SG_Disk$s_Oort[,3], colour = colnames(solution$SG_Disk$Oort)[3])) +
  #   
  #   geom_line(aes(x = solution$SG_Disk$Parameters[,4], y = solution$SG_Disk$Oort[,4], colour = colnames(solution$SG_Disk$Oort)[4]), size = 1) +
  #   geom_errorbar(aes(x = solution$SG_Disk$Parameters[,4], ymin = solution$SG_Disk$Oort[,4] - solution$SG_Disk$s_Oort[,4], ymax = solution$SG_Disk$Oort[,4] + solution$SG_Disk$s_Oort[,4], colour = colnames(solution$SG_Disk$Oort)[4])) +
  #   
  #   geom_point(aes(x = solution$SG_Disk$Parameters[,4], y = rep(0,nrow(solution$SG_Disk$Parameters))))
  # 
  # g <- g + geom_line(aes(x = solution$SG_Galo$Parameters[,4], y = solution$SG_Galo$Oort[,1], colour = colnames(solution$SG_Galo$Oort)[1]), size = 1) +
  #   geom_errorbar(aes(x = solution$SG_Galo$Parameters[,4], ymin = solution$SG_Galo$Oort[,1] - solution$SG_Galo$s_Oort[,1], ymax = solution$SG_Galo$Oort[,1] + solution$SG_Galo$s_Oort[,1], colour = colnames(solution$SG_Galo$Oort)[1])) +
  #   
  #   geom_line(aes(x = solution$SG_Galo$Parameters[,4], y = solution$SG_Galo$Oort[,2], colour = colnames(solution$SG_Galo$Oort)[2]), size = 1) +
  #   geom_errorbar(aes(x = solution$SG_Galo$Parameters[,4], ymin = solution$SG_Galo$Oort[,2] - solution$SG_Galo$s_Oort[,2], ymax = solution$SG_Galo$Oort[,2] + solution$SG_Galo$s_Oort[,2], colour = colnames(solution$SG_Galo$Oort)[2])) +
  #   
  #   geom_line(aes(x = solution$SG_Galo$Parameters[,4], y = solution$SG_Galo$Oort[,3], colour = colnames(solution$SG_Galo$Oort)[3]), size = 1) +
  #   geom_errorbar(aes(x = solution$SG_Galo$Parameters[,4], ymin = solution$SG_Galo$Oort[,3] - solution$SG_Galo$s_Oort[,3], ymax = solution$SG_Galo$Oort[,3] + solution$SG_Galo$s_Oort[,3], colour = colnames(solution$SG_Galo$Oort)[3])) +
  #   
  #   geom_line(aes(x = solution$SG_Galo$Parameters[,4], y = solution$SG_Galo$Oort[,4], colour = colnames(solution$SG_Galo$Oort)[4]), size = 1) +
  #   geom_errorbar(aes(x = solution$SG_Galo$Parameters[,4], ymin = solution$SG_Galo$Oort[,4] - solution$SG_Galo$s_Oort[,4], ymax = solution$SG_Galo$Oort[,4] + solution$SG_Galo$s_Oort[,4], colour = colnames(solution$SG_Galo$Oort)[4])) +
  #   
  #   geom_point(aes(x = solution$SG_Galo$Parameters[,4], y = rep(0,nrow(solution$SG_Galo$Parameters))))
  # 
  return(g)
}

