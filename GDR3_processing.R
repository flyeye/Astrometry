# GDR3 Processing


GDR3_apass9_save <- function(data, path = "/Catalogues/Gaia/")
{
 write.fst(ga, paste0(path, "gaia_apass9.fst"), 75)
}

# gaia <- Gaia_load_goog()
GDR3_load_good <- function(path =  "/Catalogues/Gaia/")
{
  data <- read.fst(paste0(path, "gaia_good.fst"))
  setDT(data)
  data <- CalcGalXYZ(data)
  return(data)
}

GDR3_apass9_load <- function(path = "/Catalogues/Gaia/")
{
  data <- read.fst(paste0(path, "gaia_apass9.fst"))
  setDT(data)
  data <- filter(data, !is.na(gMag))
  
  data <- CalcGalXYZ(data)
  data <- GDR3_calc_distance(data)
  data <- GDR3_calc_absolute_mag(data)
  data <- GDR3_calc_LClass(data)
  return(data)
}

GDR3_calc_gpm <- function(data)
{
  return(calc_galaxy_mu_DT(data))
}
# get_galaxy_mu(pm = ga[1:100, .(gpmRA, gpmDE)], gal = ga[1:100, .(gl, gb)], dec = ga[1:100, DE])

GDR3_calc_distance <- function(data, dist_type = "GAIA_PX")
{
  data[, R := Inf]
  
  if (dist_type == "GAIA_PX")
  {
    data[(gPx>0.01), R := 1000/gPx]
  }else if (dist_type == "HIP")
  {
    data[(hPx>0.01), R := 1000/hPx]
  }else if (dist_type == "RGEO")
  {
    data[!is.na(Rgeo), R := Rgeo]
  } else if (dist_type == "RPGEO")
  {
    data[!is.na(Rpgeo), R := Rpgeo]
  } 
  return(data)
}

GDR3_calc_absolute_mag <- function(data)
{
  data[, M := NA]
  index <-(data$R>0)&(data$R<Inf)&(!is.na(data$Mag))&(!is.na(data$B_V))
  data[index, M := Mag + 5 - 5*log10(R)]
  return (data) 
}

GDR3_calc_LClass <- function(data)
{
  
  index <- (!is.na(data$B_V)) & (!is.na(data$M))
  
  #tgas_a <- mutate(tgas_a, LClass_apass = 0)
  data[, LClass_apass := 0]
  
  #tgas_a$LClass_apass[is_main_sequence(tgas_a$B_V, tgas_a$M)] <- 5
  data[is_main_sequence(data$B_V, data$M), LClass_apass := 5]
  
  #tgas_a$LClass_apass[tgas_a$M<(-2)& (tgas_a$M>(-Inf))] <- 1
  data[(M<(-2)), LClass_apass := 1]
  
#  tgas_a$LClass_apass[(tgas_a$M>-1.5)&(tgas_a$M<2.5)&(tgas_a$B_V>0.8)&(tgas_a$B_V<2.5)] <- 3
  data[(M>(-1.5)&(M<(2.5))&(B_V>0.8)&(B_V<2.5)), LClass_apass := 3]
  
  # ctgas_ <- within(tgas_, rm("LClass_apass")) # удаляет строчки в каких то случаях...
  #tgas_$LClass_apass <- NULL
  #tgas_ <- tgas_ %>% left_join(tgas_a[ , names(tgas_a) %in% c("source_id", "LClass_apass")], by = "source_id")
  #tgas_$LClass_apass[is.na(tgas_$LClass_apass)] <- 0
  
  return(data)
}

# min_px, max_px - mas
GDR3_filter <- function(data,                     #  catalog (data frame)
                          px = c(-Inf,Inf),        #  limits on parallax
                          e_px = Inf,              #  limits on parallax error
                          bv_lim = c(-Inf,Inf),    #  limits on B-V
                          Mg = c(-Inf,Inf),        #  limits on absolute magnitude applied to M column in data
                          lclass = c(0,1,2,3,4,5),
                          z_lim = c(0, Inf),       #  limits on distance from galactic equator
                          r_lim = c(0, Inf),       #  limits on solar distance applied to R column in data (parsec)
                          g_b = c(-Inf, Inf))      #  limits on galactic lalitude
{
  #print(bv_lim)
  #print(Mg)
  #print(px)
  #print(e_px)
  #cat(px, "\n")
  #cat(r_lim, "\n")
  cat(g_b)
  
  data <- data %>% filter(!is.na(B_V) & !is.na(M)) %>%
   filter( (gPx > px[1]) & (gPx <= px[2])) %>%  #mas
   filter( (R>r_lim[1]) & (R<=r_lim[2]) ) %>%  #pc
   filter( (B_V>bv_lim[1]) & (B_V<=bv_lim[2]) ) %>%
   filter( (M>Mg[1]) & (M<Mg[2])) %>%
   filter( LClass_apass %in% lclass) %>%
   filter( (parallax_error/gPx) < e_px ) %>%
   filter( (abs(z)>=z_lim[1]) & (abs(z)<z_lim[2])) %>% #kpc
   filter( (gb>g_b[1]) & (gb<=g_b[2]))  #grad
  
  # setkey(data, B_V, M, R, LClass_apass)
  # data <- data[!is.na(data$B_V) & !is.na(data$M),][data$LClass_apass %in% LClass,][(data$gPx > px[1]) & (data$gPx <= px[2]),][(data$R>r_lim[1]) & (data$R<=r_lim[2]),][(data$B_V>bv_lim[1]) & (data$B_V<=bv_lim[2]),][(data$M>Mg[1]) & (data$M<Mg[2]),][(data$parallax_error/data$gPx) < e_px,][(abs(data$z)>=z_lim[1]) & (abs(data$z)<z_lim[2]),][(data$gb>=g_b[1]) & (data$gb<g_b[2]),]

  print(paste("stars in sample:", nrow(data)))
  return(data)
}

GDR3_get_stars <- function(data, src = "TGAS")
{
  stars <- matrix(0,nrow(data), 6)
  stars[,1] <- data$gl
  stars[,2] <- data$gb
  stars[,3] <- data$R/1000
  #if (src == "TGAS")
  # {
     stars[,4] <- data$gpm_l
     stars[,5] <- data$gpm_b
  # } else if(src == "TYCHO")
  # {
  #   stars[,4] <- data$tpm_l
  #   stars[,5] <- data$tpm_b
  # } else
  # {
  #   stars[,4] <- data$pm_l
  #   stars[,5] <- data$pm_b
  # }
  stars[,6] <- data$Vr
  colnames(stars) <- c("l", "b", "r", "pml", "pmb", "vr")
  
  return(stars)
}


GDR3_ph_saturation_correction <- function(data)
{
  #2.0 < G < 8
  data[, gMag_corr := gMag]
  data[(gMag>2.0)&(gMag<8.0), gMag_corr := gMag - 0.09892 + 0.059*gMag   - 0.00977*(gMag^2) + 0.0004934*(gMag^3)]
  
  #2.0 < G < 3.94
  data[, bMag_corr := bMag]
  data[(gMag>2.0)&(gMag<3.94), bMag_corr := bMag - 0.9921  - 0.02598*gMag + 0.1833*(gMag^2)  - 0.02862*(gMag^3)]
  
  return (data)
}

### =======================================================
### -------------------------------------------------------

GDR3_calc_OM_seq_2 <- function(data, src_ = "TGAS", px_type = "DIST", distance = c(1, 1000), save = NULL, type = 1, model = 1, bv = c(-Inf, Inf),  
                               dist_type = "GAIA_PX", use = c(TRUE, TRUE, FALSE), lclass = c(0,1,2,3,4,5), g_b = c(-Inf, Inf), z_lim = c(-Inf, Inf), ...)
{
  distance <- matrix(distance, ncol = 2)
  bv <- matrix(bv, ncol = 2)
  q <- nrow(distance) * nrow(bv)
  
  res <- matrix(0, q, 11)
  err <- matrix(0, q, 11)
  par <- matrix(0, q, 10)
  sol <- matrix(0, q, 3)
  colnames(sol) <- c("X", "Y", "Z")
  oort <- matrix(0, q, 6)
  oort_err <- matrix(0, q, 6) 
  
  solution <- list();
  
  for (k in 1:nrow(distance))
  {
    
    for (j in 1:nrow(bv))
    {
      
      gc()
      
      i <- (k-1) * nrow(bv) + j
      
      par[i, 1] <- distance[k,1]
      par[i, 2] <- distance[k,2]
      
      par[i, 7] <- bv[j,1]
      par[i, 8] <- bv[j,2]
      bvl <- c(par[i, 7], par[i, 8])
      
      if (px_type!="ANGLE")
      {
        px_ <- c(-1, Inf)   #c(1/par[i,2], 1/par[i,1])
        dist_ <- c(par[i,1], par[i,2])*1000
        colnames(par) <- c("r_min","r_max","number of stars","r_mean", "px_err", "r_mean_model", "B-V_min", "B-V_max", "B-V_mean", "s0")
      } else
      {
        px_ <- c(par[i,1], par[i,2])
        dist_ <- c(0, Inf)
        colnames(par) <- c("px_min","px_max","number of stars","px_mean", "px_err", "r_mean_model", "B-V_min", "B-V_max", "B-V_mean", "s0")
      }
      
      cat("filtering...", "\n")
      cat("B-V", bvl, "\n")
      cat("dist", dist_, "\n")
      cat("lclass", lclass, "\n")
      cat("px", px_, "\n")
      cat("g_b", g_b, "\n")
      sample <- GDR3_filter(data, px = px_, bv_lim = bvl, r_lim = dist_, lclass = lclass, g_b = g_b, z_lim = z_lim);
      
      if (nrow(sample)<12)
      {
        cat("Not enough samples in data:", nrow(sample), "\n")
        par[i,3:9] <- 0
        res[i, ] <- c(rep(0, ncol(res)))
        err[i, ] <- c(rep(0, ncol(err)))
        sol[i, ] <- c(rep(0, ncol(sol)))
        next()
      } else 
      {
        cat("Stars after filetering:", nrow(sample),"\n")
      }
      
      #cat("EQ to GAL converting...", "\n")
      #tgas_sample <- tgas_calc_gpm(tgas_sample)
      #cat("Stars after EQ to gal converting:", nrow(tgas_sample),"\n")
      
      par[i,4] <- mean(sample$R)/1000
      cat("Mean filtering distance:", par[i,4], "\n")
      
      par[i,9] <- mean(sample$B_V)
      cat("Mean sample B-V:", par[i,9], "\n")
      
      sample <- GDR3_calc_distance(sample, dist_type)
      sample <- sample %>% filter(R<50000)
      
      #tgas_sample$M <- (tgas_sample$apasm_v + 5 - 5*log10(tgas_sample$R))
      cat("Stars after distance re-calc:", nrow(sample),"\n")
      
      # sample_ <<- tgas_sample
      
      if(!is.null(save))
      {
        cat("Saving...", "\n")
        s <- paste0(save, "_", src_, "_", par[i,1], "-",par[i,2], "_", par[i,7], "-",par[i,8])
        #write_csv(select(sample, RA, DE, gpmRA, gpmDE, gl, gb, R, gpm_l, gpm_b, apasm_b, apasm_v), paste0(s, "_sample.csv"), col_names = TRUE)
        #x = as.matrix(select(sample, RA, DE, gpmRA, gpmDE, tyc_pmRA, tyc_pmDE, gl, gb, R, gpm_l, gpm_b, apasm_b, apasm_v, parallax_error))
        #write.fwf(x, file = paste0(s, "_sample.txt"), colnames = TRUE, sep = "   ")
        
        cat("draw HR...", "\n")
        hrd <- HRDiagramGDR3(sample, save = s)
        
        max_dist <- max(sample$R/1000) %/% 1 + 1;
        #cat(max_dist)
        cat("draw XY...", "\n")
        g <- DrawGalaxyPlane(sample, plane = "XY", save = s, dscale = max_dist)
        rm(g)
        
        cat("draw XZ...", "\n")
        g <- DrawGalaxyPlane(sample, plane = "XZ", save = s, dscale = max_dist)
        rm(g)
        
        cat("draw YZ...", "\n")
        g <- DrawGalaxyPlane(sample, plane = "YZ", save = s, dscale = max_dist)
        rm(g)
        
      }
      
      cat("Equation of conditions forming...", "\n")
      stars <- GDR3_get_stars(sample, src_)
      cat("stars N:", nrow(stars), "\n")
      # stars_ <<- stars;
      
      par[i,3] <- nrow(stars)
      par[i,6] <- mean(stars[,3])
      par[i,5] <- mean(sample$parallax_error)
      cat("Solution calculation...", "\n")
      
      res_tgas <- Calc_OM_Model(stars, use = use, mode = 2, model = model, type = type)
      # res <<- res_tgas
      
      if ((i == 1) & (ncol(res) != length(res_tgas$X)))
      {
        res <- matrix(0, q, length(res_tgas$X))
        err <- matrix(0, q, length(res_tgas$X))
      }
      
      res[i, ] <- res_tgas$X
      err[i, ] <- res_tgas$s_X
      sol[i, ] <- res_tgas$X[1:3]/par[i,6]
      oort[i, ] <- res_tgas$Oort
      oort_err[i, ] <- res_tgas$s_Oort
      cat(par[i,], "\n")
      cat(res_tgas$X, "\n")
      cat(res_tgas$s_X, "\n")
      #if(!is.null(save))
      #  res_tgas$HR <- hrd;
      solution[[i]] <- res_tgas
      par[i, 10] <- res_tgas$s0
      
      rm(tgas_res)
      rm(hrd)
      rm(sample)
    }
  }
  colnames(res) <- names(res_tgas$X)
  colnames(err) <- names(res_tgas$s_X)
  colnames(oort) <- names(res_tgas$Oort)
  colnames(oort_err) <- names(res_tgas$s_Oort)
  
  res <- list(X = res, S_X = err, Sol = sol, Parameters = par, Oort = oort, s_Oort = oort_err)
  res$SolutionR <- solution;
  
  return(res)
}



GDR3_ph_correction_LinYang <- function(data)
{
  G <- c(10.0, 10.4, 10.8, 11.2, 11.6, 12.0, 12.4, 12.5, 12.6, 12.7, 12.8, 12.9, 13.0, 13.1, 13.2, 13.3, 13.4, 
         13.5, 13.6, 13.7, 13.8, 13.9, 14.0, 14.1, 14.2, 14.3, 14.4, 14.5, 14.6, 14.7, 14.8, 14.9, 15.0, 15.1,
         15.2, 15.3, 15.4)
  
  dG <- c(???6.7, ???5.9, ???5.0, ???4.1, ???3.2, ???2.3, ???1.0, ???0.6, ???0.4, ???0.4, ???0.2, ???0.9, ???1.4, ???1.9, ???2.7, 
           ???4.0, ???3.5, ???3.5, ???3.6, ???2.6, ???2.7, ???3.4, ???3.5, ???3.7, ???3.7, ???3.5, ???4.2, ???4.1, ???4.4, ???3.9,
           ???3.3, ???3.3, ???3.0, ???3.4, ???3.1, ???3.2, ???3.3, ???3.2, ???3.1, ???3.0, ???2.9, ???2.8, ???2.8, ???2.9, ???3.1, 
           ???3.3, ???3.5, ???3.6, ???3.7, ???3.4, ???3.5, ???3.2, ???3.0, ???2.8, ???2.7, ???2.4, ???2.1, ???1.7, ???1.6, ???1.6, 
           ???1.5, ???1.9, ???2.3, ???2.6, ???2.9, ???2.8, ???2.8, ???3.0, ???3.5, ???4.0)
           
  dGbp <- c(???14.9, ???14.0, ???13.1, ???12.3, ???11.4, ???10.6, ???10.1, ???10.0, ???9.7, ???9.1, ???8.4, ???8.1, ???8.4, ???8.9,
             ???10.0, ???10.6, ???10.8, ???11.0, ???10.8, ???10.5, ???10.3, ???10.0, ???9.7, ???9.7, ???9.8, ???9.9, ???9.9, ???10.0,
             ???9.9, ???9.9, ???9.3, ???9.2, ???9.7, ???9.8, ???9.8, ???9.9, ???10.0, ???10.1, ???10.6, ???10.8, ???10.7, ???10.7, ???10.8, 
             ???11.0, ???11.6, ???12.0, ???12.4, ???12.8, ???13.2, ???13.5, ???13.8, ???14.1, ???14.5, ???14.7, ???15.2, ???15.0, ???16.0, 
             ???16.4, ???16.8, ???17.2, ???17.7, ???18.4, ???19.1, ???20.0, ???20.3, ???20.3, ???20.3, ???20.5, ???20.8, ???21.8) 
             
  dGrp <- c(6.3, 5.7, 5.2, 4.6, 4.1, 3.5, 3.1, 3.1, 3.2, 3.1, 3.0, 2.9, 2.5, 1.9, 1.1, 0.2, 0.1, 0.1, 
            0.1, 0.1, 0.1, 0.0, ???0.1, 0.0, 0.2, 0.3, 0.3, 0.8, 1.3, 1.8, 2.2, 2.5, 2.4, 2.4, 2.3, 2.7,
            3.1, 3.3, 3.3, 3.3, 3.2, 3.2, 3.1, 2.9, 2.6, 2.2, 1.7, 1.2, 0.8, 0.5, 0.1, ???0.2, ???0.5, ???0.7, 
            ???0.8, ???0.9, ???1.0, ???1.3, ???1.8, ???2.3, ???3.1, ???3.8, ???4.4, ???4.7, ???4.9, ???4.9, ???5.1, ???5.2, ???5.9, ???8.6)
}

