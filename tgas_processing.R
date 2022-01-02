

# min_px, max_px - mas
filter_tgs_px <- function(tgs,                     #  catalog (data frame)
                          px = c(-Inf,Inf),        #  limits on parallax
                          e_px = Inf,              #  limits on parallax error
                          bv_lim = c(-Inf,Inf),    #  limits on B-V
                          Mg = c(-Inf,Inf),        #  limits on absolute magnitude applied to M column in data
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
  
  tgs <- CalcGalXYZ(tgs)
  tgs <- tgs %>% filter(!is.na(B_V) & !is.na(M)) %>%
    filter( (gPx > px[1]) & (gPx <= px[2])) %>%  #mas
    filter( (R>r_lim[1]) & (R<=r_lim[2]) ) %>%  #pc
    filter( (B_V>bv_lim[1]) & (B_V<=bv_lim[2]) ) %>%
    filter( (M>Mg[1]) & (M<Mg[2])) %>%
    filter( (parallax_error/gPx) < e_px ) %>% 
    filter( (abs(z)>=z_lim[1]) & (abs(z)<z_lim[2])) %>% #kpc
    filter( (gb>=g_b[1]) & (gb<g_b[2]))  #grad
  
  print(paste("stars in sample:", nrow(tgs)))
  return(tgs)
}

tgas_get_stars <- function(tgas_data, src = "TGAS")
{
  stars <- matrix(0,nrow(tgas_data), 6)
  stars[,1] <- tgas_data$gl
  stars[,2] <- tgas_data$gb
  #stars[,3] <- (1/tgas_data$gPx)    #kPc
  stars[,3] <- tgas_data$R/1000
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
  # GAIA proper motions
  tgas_$pmRA <- tgas_$gpmRA
  tgas_$pmDE <- tgas_$gpmDE
  tgas_ <- cat_eq2gal(tgas_)
  tgas_ <- mutate(tgas_, gpm_l = pm_l, gpm_b = pm_b)
  
  # TYCHO-2 proper motions
  #tgas_ <- filter(tgas_ , !is.na(tyc_pmRA))
  tgas_$pmRA <- tgas_$tyc_pmRA
  tgas_$pmDE <- tgas_$tyc_pmDE
  tgas_ <- cat_eq2gal(tgas_)
  tgas_ <- mutate(tgas_, tpm_l = pm_l, tpm_b = pm_b)
  
  return(tgas_)
}


tgas_calc_distance <- function(tgas_, dist_type = "TGAS_PX")
{
  tgas_$R <- 0;
  
  if (dist_type == "TGAS_PX")
  {
    tgas_$R[tgas_$gPx>0.01] <- 1000/tgas_$gPx[tgas_$gPx>0.01]
    tgas_$R[tgas_$gPx<=0.01] <- Inf
  }else if (dist_type == "HIP")
  {
    tgas_$R[tgas_$Px>0.01] <- 1000/tgas_$Px[tgas_$Px>0.01]
    tgas_$R[tgas_$Px<=0.01] <- Inf
  }else if (dist_type == "rMoMW")
  {
    tgas_ <- mutate(tgas_, R = rMoMW)
  } else if (dist_type == "rMoExp2")
  {
    tgas_ <- mutate(tgas_, R = rMoExp2)
  } else if (dist_type == "rMoExp1")
  {
    tgas_ <- mutate(tgas_, R = rMoExp1)
  }
  return(tgas_)
}

tgas_calc_absolute_mag <- function(tgas_)
{
  index <- (!is.na(tgas_$R)) & (!is.na(tgas_$Mag))   
  tgas_$M[index] <- (tgas_$Mag[index] + 5 - 5*log10(tgas_$R[index]))
  #tgas_$M[index] <- tgas_$apasm_v[index] + 5 + 5*log10(tgas_$gPx[index]/1000)  
  return (tgas_) 
}

tgas_calc_LClass <- function(tgas_, dist_ = "TGAS_PX", ph = "APASS")
{
  if (ph == "APASS")
  {
    tgas_ <- tgas_apply_APASS(tgas_)
  } else if (ph == "HIP")
  {
    tgas_ <- tgas_apply_HIP_photometry(tgas_)  
  } else
  {
    tgas_ <- tgas_apply_APASS(tgas_)
    tgas_ <- tgas_apply_HIP_photometry(tgas_, reset = FALSE)  
  }
  
  tgas_ <- tgas_calc_distance(tgas_, dist_) 
  tgas_ <- tgas_calc_absolute_mag(tgas_)
  
  tgas_a <- tgas_[(!is.na(tgas_$B_V)) & (!is.na(tgas_$M)),]
  
  tgas_a <- mutate(tgas_a, LClass_apass = 0)
  tgas_a$LClass_apass[is_main_sequence(tgas_a$B_V, tgas_a$M)] <- 5
  tgas_a$LClass_apass[tgas_a$M<(-2)& (tgas_a$M>(-Inf))] <- 1
  tgas_a$LClass_apass[(tgas_a$M>-1.5)&(tgas_a$M<2.5)&(tgas_a$B_V>0.8)&(tgas_a$B_V<2.5)] <- 3
  
  # ctgas_ <- within(tgas_, rm("LClass_apass")) # удал€ет строчки в каких то случа€х...
  tgas_$LClass_apass <- NULL
  tgas_ <- tgas_ %>% left_join(tgas_a[ , names(tgas_a) %in% c("source_id", "LClass_apass")], by = "source_id")
  tgas_$LClass_apass[is.na(tgas_$LClass_apass)] <- 0
  
  return(tgas_)
}

min_M <- function(bv)
{
  m <- numeric(length(bv))
  
  # ¬ерхн€€ граница подн€та относительно Ѕови
  i1 <- bv < 0.75
  m[i1] <- 5.73 * bv[i1] - 2.32222;
  
  i2 <- (bv>=0.75)&(bv<1.0)
  m[i2] <- (16 * bv[i2] - 10)
  
  #i1 <- bv < 0.8
  #m[i1] <- 5.7 * bv[i1] - 1.22222;
  
  #i2 <- (bv>=0.8)&(bv<1.0)
  #m[i2] <- (13 * bv[i2] - 7)
  
  
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


### =======================================================
### -------------------------------------------------------

tgas_calc_OM_seq_2 <- function(tgas_ = tgas, src_ = "TGAS", px_type = "ANGLE", distance = c(0.001, 100), save = NULL, type = 0, model = 1, bv = c(-Inf, Inf),  dist_type = "TGAS_PX", use = c(TRUE, TRUE, FALSE),...)
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
      tgas_sample <- filter_tgs_px(tgas_, px = px_, bv_lim = bvl, r_lim = dist_, ...);

      if (nrow(tgas_sample)<12)
      {
        cat("Not enough samples in data:", nrow(tgas_sample), "\n")
        par[i,3:9] <- 0
        res[i, ] <- c(rep(0, 11))
        err[i, ] <- c(rep(0, 11))
        sol[i, ] <- c(rep(0, 3))
        next()
      } else 
      {
        cat("Stars after filetering:", nrow(tgas_sample),"\n")
      }
      
      cat("EQ to GAL converting...", "\n")
      tgas_sample <- tgas_calc_gpm(tgas_sample)
      cat("Stars after EQ to gal converting:", nrow(tgas_sample),"\n")
      
      par[i,4] <- mean(tgas_sample$R)/1000
      cat("Mean filtering distance:", par[i,4], "\n")
      
      par[i,9] <- mean(tgas_sample$B_V)
      cat("Mean sample B-V:", par[i,7], "\n")
      
      tgas_sample <- tgas_calc_distance(tgas_sample, dist_type)
      tgas_sample <- tgas_sample %>% filter(R<50000)
      
      #tgas_sample$M <- (tgas_sample$apasm_v + 5 - 5*log10(tgas_sample$R))
      cat("Stars after distance re-calc:", nrow(tgas_sample),"\n")
      
      # sample_ <<- tgas_sample
      
      if(!is.null(save))
      {
        cat("Saving...", "\n")
        s <- paste0(save, "_", src_, "_", par[i,1], "-",par[i,2], "_", par[i,7], "-",par[i,8])
        # write_csv(select(tgas_sample, RA, DE, gpmRA, gpmDE, tyc_pmRA, tyc_pmDE, gl, gb, R, gpm_l, gpm_b, tpm_l, tpm_b, apasm_b, apasm_v), paste0(s, "_sample.csv"), col_names = TRUE)
        x = as.matrix(select(tgas_sample, RA, DE, gpmRA, gpmDE, tyc_pmRA, tyc_pmDE, gl, gb, R, gpm_l, gpm_b, tpm_l, tpm_b, apasm_b, apasm_v, parallax_error, TYC))
        write.fwf(x, file = paste0(s, "_sample.txt"), colnames = TRUE, sep = "   ")
        
        cat("draw HR...", "\n")
        hrd <- HRDiagram(tgas_sample, save = s, photometric = "none")
        
        max_dist <- max(tgas_sample$R/1000) %/% 1 + 1;
        #cat(max_dist)
        cat("draw XY...", "\n")
        DrawGalaxyPlane(tgas_sample, plane = "XY", save = s, dscale = max_dist)
        cat("draw XZ...", "\n")
        DrawGalaxyPlane(tgas_sample, plane = "XZ", save = s, dscale = max_dist)
        cat("draw YZ...", "\n")
        DrawGalaxyPlane(tgas_sample, plane = "YZ", save = s, dscale = max_dist)
      }
      
      cat("Equation of conditions forming...", "\n")
      stars <- tgas_get_stars(tgas_sample, src_)
      # stars_ <<- stars;
      
      par[i,3] <- nrow(stars)
      par[i,6] <- mean(stars[,3])
      par[i,5] <- mean(tgas_sample$parallax_error)
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

### -------------------------------------------------------
### =======================================================


tgas_calc_OM_seq <- function(tgas_ = tgas, src_ = "TGAS", start = 1, step = 0.1, q = 2, px_type = "DIST", distance = NULL, save = NULL, type = 0, model = 1, dist_type = "TGAS_PX", use = c(TRUE, TRUE, FALSE),...)
{
  if (!is.null(distance))
    q <- nrow(distance)
  
  res <- matrix(0, q, 11)
  err <- matrix(0, q, 11)
  par <- matrix(0, q, 6)
  sol <- matrix(0, q, 3)
  colnames(sol) <- c("X", "Y", "Z")
  oort <- matrix(0, q, 6)
  oort_err <- matrix(0, q, 6)
  
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
      px_ <- c(-1, Inf)   #c(1/par[i,2], 1/par[i,1])
      dist_ <- c(par[i,1], par[i,2])*1000
      colnames(par) <- c("r_min","r_max","number of stars","r_mean", "px_err", "r_mean_model")
    } else
    {
      px_ <- c(par[i,1], par[i,2])
      dist_ <- c(0, Inf)
      colnames(par) <- c("px_min","px_max","number of stars","px_mean", "px_err", "r_mean_model")
    }
    
    
    cat("filtering...", "\n")
    tgas_sample <- filter_tgs_px(tgas_, px = px_, r_lim = dist_, bv_lim = bv, ...);
    if (nrow(tgas_sample)<12)
    {
      cat("Not enough samples in data:", nrow(tgas_sample), "\n")
      par[i,3] <- 0
      par[i,4] <- 0
      par[i,5] <- 0
      par[i,6] <- 0
      res[i, ] <- c(rep(0, 11))
      err[i, ] <- c(rep(0, 11))
      sol[i, ] <- c(rep(0, 3))
      next()
    } else 
    {
      cat("Stars after filetering:", nrow(tgas_sample),"\n")
    }
    
    cat("EQ to GAL converting...", "\n")
    tgas_sample <- tgas_calc_gpm(tgas_sample)
    cat("Stars after EQ to gal converting:", nrow(tgas_sample),"\n")
    
    par[i,4] <- mean(tgas_sample$R)/1000
    cat("Mean filtering distance:", par[i,4], "\n")
    
    tgas_sample <- tgas_calc_distance(tgas_sample, dist_type)
    tgas_sample <- tgas_sample %>% filter(R<50000)
    
    #tgas_sample$M <- (tgas_sample$apasm_v + 5 - 5*log10(tgas_sample$R))
    cat("Stars after distance re-calc:", nrow(tgas_sample),"\n")
    
    #sample_ <<- tgas_sample
    
    if(!is.null(save))
    {
      cat("Saving...", "\n")
      s <- paste0(save, "_", src_, "_", par[i,1], "-",par[i,2])
      # write_csv(select(tgas_sample, RA, DE, gpmRA, gpmDE, tyc_pmRA, tyc_pmDE, gl, gb, R, gpm_l, gpm_b, tpm_l, tpm_b, apasm_b, apasm_v), paste0(s, "_sample.csv"), col_names = TRUE)
      x = as.matrix(select(tgas_sample, RA, DE, gpmRA, gpmDE, tyc_pmRA, tyc_pmDE, gl, gb, R, gpm_l, gpm_b, tpm_l, tpm_b, apasm_b, apasm_v, parallax_error, TYC))
      write.fwf(x, file = paste0(s, "_sample.txt"), colnames = TRUE, sep = "   ")
      
      hrd <- HRDiagram(tgas_sample, save = s, photometric = "none")
      DrawGalaxyPlane(tgas_sample, plane = "XY", save = s, dscale = ds)
      DrawGalaxyPlane(tgas_sample, plane = "XZ", save = s, dscale = ds)
      DrawGalaxyPlane(tgas_sample, plane = "YZ", save = s, dscale = ds)
    }
    
    cat("Equation of conditions forming...", "\n")
    stars <- tgas_get_stars(tgas_sample, src_)
    # stars_ <<- stars;
    
    par[i,3] <- nrow(stars)
    par[i,6] <- mean(stars[,3])
    par[i,5] <- mean(tgas_sample$parallax_error)
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
    if(!is.null(save))
    {
      res_tgas$HR <- hrd;
    }
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


# BV <- matrix(0, nrow = 1, ncol = 2)
# BV[,1] <- c(-Inf, -0.30, 0.00, 0.30, 0.58, 0.85, 1.42)
# BV[,2] <- c(-0.30, 0.00, 0.30, 0.58, 0.85, 1.42, Inf)

tgas_calc_OM_cond <- function(tgas_ = tgas, lclass = 3, population = "ALL", src = "TGAS", 
                              type = 0, model = 1, 
                              use = c(TRUE, TRUE, FALSE), dist_type = "TGAS_PX", filter_dist = "TGAS_PX", 
                              saveto = "", 
                              g_b = c(-Inf, Inf))
{
  
  #APASS photometry
  
  conditions <- list();
  
  conditions$Src <- src;
  conditions$Dist_Type <- dist_type;
  conditions$Filter_Dist <- filter_dist;
  conditions$use <- use;
  conditions$KinModel <- model
  conditions$KinModelType <- type
  conditions$g_B <- g_b
  conditions$SaveTo <- saveto
  
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
      if (conditions$Filter_Dist == "TGAS_PX")
      {
        distance_ <- matrix(0, nrow = 4, ncol = 2)
        distance_[,1] <- c(0.0, 2.5, 5.0, 7.5)
        distance_[,2] <- c(2.5, 5.0, 7.5, 10.0)
      } else 
      {
        distance_ <- matrix(0, nrow = 1, ncol = 2)
        distance_[,1] <- c(0.0 )
        distance_[,2] <- c(10.0)
      }
    } else if (population == "GALO")
    {
      
      if (conditions$Filter_Dist == "TGAS_PX")
      {
        distance_ <- matrix(0, nrow = 7, ncol = 2)
        distance_[,1] <- c(0.0, 2.5, 5.0, 7.5, 10.0, 15.0, 20.0)
        distance_[,2] <- c(2.5, 5.0, 7.5, 10.0, 15.0, 20.0, 50.0)
      } else 
      {
        distance_ <- matrix(0, nrow = 1, ncol = 2)
        distance_[,1] <- c(0.0)
        distance_[,2] <- c(10.0)
      }
      
    } else 
    {
      if (conditions$Filter_Dist == "TGAS_PX")
      {
        distance_ <- matrix(0, nrow = 1, ncol = 2)
        distance_[,1] <- c(0.0, 10.0, 15.0, 20.0)
        distance_[,2] <- c(2.5, 5.0, 7.5, 10.0, 15.0, 20.0, 50.0)
      } else 
      {
        distance_ <- matrix(0, nrow = 1, ncol = 2)
        distance_[,1] <- c(0.0)
        distance_[,2] <- c(10.0)
      }
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
      if (conditions$Filter_Dist == "TGAS_PX")
      {
        distance_ <- matrix(0, nrow = 14, ncol = 2)
        distance_[,1] <- c(0.2, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.5, 1.75, 2.0)
        distance_[,2] <- c(0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.5, 1.75, 2.0, 3)
      } else 
      {
        distance_ <- matrix(0, nrow = 14, ncol = 2)
        distance_[,1] <- c(0.2, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.5, 1.75, 2.0)
        distance_[,2] <- c(0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.5, 1.75, 2.0, 3)
      }
    } else if (population == "GALO")
    {
      
      if (conditions$Filter_Dist == "TGAS_PX")
      {
        distance_ <- matrix(0, nrow = 16, ncol = 2)
        distance_[,1] <- c(0.5, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.6, 1.8, 2.0, 2.25, 2.5, 2.75, 3)
        distance_[,2] <- c(0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.6, 1.8, 2.0, 2.25, 2.5, 2.75, 3, 4)
      } else 
      {
        distance_ <- matrix(0, nrow = 12, ncol = 2)
        distance_[,1] <- c(0.5, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.6, 1.8, 2.0)
        distance_[,2] <- c(0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.6, 1.8, 2.0, 3)
      }
    }else 
    {
      if (conditions$Filter_Dist == "TGAS_PX")
      {
        distance_ <- matrix(0, nrow = 20, ncol = 2)
        distance_[,1] <- c(0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.6, 1.8, 2.0, 2.25, 2.5, 2.75, 3)
        distance_[,2] <- c(0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.6, 1.8, 2.0, 2.25, 2.5, 2.75, 3, 4)
      } else 
      {
        distance_ <- matrix(0, nrow = 16, ncol = 2)
        distance_[,1] <- c(0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.6, 1.8, 2.0 )
        distance_[,2] <- c(0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.6, 1.8, 2.0, 3)
      }
    }
    
  } else if(lclass == 5)
  {
    tgas_ <- tgas_[tgas_$LClass_apass == 5,]
    
    conditions$BV <- c(-Inf, Inf)
    conditions$MG <- c(-Inf, Inf)
    conditions$e_Px <- Inf
    conditions$Z <- c(0, Inf)
    
    distance_ <- matrix(0, nrow = 14, ncol = 2)
    distance_[,1] <- c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.5)
    distance_[,2] <- c(0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.5, 2.0)
  } 
  
  con <- file(paste0(conditions$SaveTo,"description.txt"), "w")
  cat("B-V = ", conditions$BV, "\n", file=con)
  cat("M = ", conditions$MG, "\n", file=con)
  cat("ePX = ", conditions$e_Px, "\n", file=con)
  cat("z = ", conditions$Z, "\n", file=con)
  cat("LClass = ", conditions$LClass , "\n", file=con)
  cat("Distance = ", distance_ , "\n", file=con)
  cat("Population = ", conditions$Population , "\n", file=con)
  cat("Distance type = ", conditions$Dist_Type, "\n", file=con)
  cat("Filtering Distance type = ", conditions$Filter_Dist, "\n", file=con)
  cat("used equations (mu_l, mu_b, v_r) = ", conditions$use, "\n", file=con)
  cat("Kinematic model = ", conditions$KinModel, "\n", file=con)
  cat("Kinematic model type = ", conditions$KinModelType, "\n", file=con)
  cat("gal B = ", conditions$g_B, "\n", file = con)
  
  close(con)
  
  solution  <- tgas_calc_OM_seq(tgas_, src_ = src, 
                                start = start, step = step, q = q, 
                                z_lim = conditions$Z, e_px = conditions$e_Px, bv = conditions$BV, Mg = conditions$MG, 
                                px_type = "DIST", distance = distance_, 
                                save = conditions$SaveTo, 
                                type = conditions$KinModelType, model = conditions$KinModel, 
                                dist_type = dist_type, use = conditions$use, 
                                g_b = conditions$g_B)
  
  solution$Conditions <- conditions
  
  return(solution)
  
}


tgas_make_Oort_solutions <- function(dist_type = "TGAS_PX", filter_dist = "TGAS_PX", src = "TGAS", name = "")
{
  tgas <- tgas_calc_LClass(tgas, dist_ = filter_dist)
  
  if (!dir.exists("solutions")) 
    dir.create("solutions")
  saveto_ <- paste0("solutions/solution_", name, "_", filter_dist, "-", dist_type)
  if (!dir.exists(saveto_)) 
    dir.create(saveto_)
  
  solutions <- list();
  
  # --------------------------------------
  
  saveto_2 <- paste0(saveto_, "/sphere")
  if (!dir.exists(saveto_2)) 
    dir.create(saveto_2)
  
  cat("Red Giants processing", "\n")
  saveto_3 <- paste0(saveto_2, "/RG_ALL")
  if (!dir.exists(saveto_3)) 
    dir.create(saveto_3)
  solutions$RG_All <-  tgas_calc_OM_cond(tgas, lclass = 3, population = "ALL", 
                                         model = 2, type = 2, 
                                         src = src, 
                                         use = c(TRUE, FALSE, FALSE),
                                         dist_type = dist_type, filter_dist = filter_dist, 
                                         saveto = paste0(saveto_3, "/"))
  solutions$RG_All$Name <- "Red Giants"
  
  cat("Main Sequence processing", "\n")
  saveto_3 <- paste0(saveto_2, "/MS_ALL")
  if (!dir.exists(saveto_3)) 
    dir.create(saveto_3)
  solutions$MS_All <-  tgas_calc_OM_cond(tgas, lclass = 5, population = "ALL", 
                                         model = 2, type = 2, 
                                         use = c(TRUE, FALSE, FALSE),
                                         src = src, 
                                         dist_type = dist_type, filter_dist = filter_dist, 
                                         saveto = paste0(saveto_3, "/"))
  solutions$MS_All$Name <- "Main Sequence"
  
  # --------------------------------------
  
  saveto_2 <- paste0(saveto_, "/north_hemisphere")
  if (!dir.exists(saveto_2)) 
    dir.create(saveto_2)
  
  cat("Red Giants processing North", "\n")
  saveto_3 <- paste0(saveto_2, "/RG_North")
  if (!dir.exists(saveto_3)) 
    dir.create(saveto_3)
  solutions$RG_North <-  tgas_calc_OM_cond(tgas, lclass = 3, population = "ALL", 
                                           model = 2, type = 3, 
                                           use = c(TRUE, FALSE, FALSE),
                                           src = src, 
                                           dist_type = dist_type, filter_dist = filter_dist, 
                                           g_b = c(0, Inf),
                                           saveto = paste0(saveto_3, "/"))
  solutions$RG_North$Name <- "Red Giants North"
  
  cat("Main Sequence processing North", "\n")
  saveto_3 <- paste0(saveto_2, "/MS_North")
  if (!dir.exists(saveto_3)) 
    dir.create(saveto_3)
  solutions$MS_North <-  tgas_calc_OM_cond(tgas, lclass = 5, population = "ALL", 
                                           model = 2, type = 3, 
                                           use = c(TRUE, FALSE, FALSE),
                                           src = src, 
                                           dist_type = dist_type, filter_dist = filter_dist, 
                                           g_b = c(0, Inf),
                                           saveto = paste0(saveto_3, "/"))
  
  solutions$MS_North$Name <- "Main Sequence North"
  
  # --------------------------------------
  
  saveto_2 <- paste0(saveto_, "/south_hemisphere")
  if (!dir.exists(saveto_2)) 
    dir.create(saveto_2)
  cat("Red Giants processing", "\n")
  saveto_3 <- paste0(saveto_2, "/RG_South")
  if (!dir.exists(saveto_3)) 
    dir.create(saveto_3)
  solutions$RG_South <-  tgas_calc_OM_cond(tgas, lclass = 3, population = "ALL", 
                                           model = 2, type = 3, 
                                           use = c(TRUE, FALSE, FALSE),
                                           src = src, 
                                           dist_type = dist_type, filter_dist = filter_dist, 
                                           g_b = c(-Inf, 0),
                                           saveto = paste0(saveto_3, "/"))
  solutions$RG_South$Name <- "Red Giants South"
  
  cat("Main Sequence processing", "\n")
  saveto_3 <- paste0(saveto_2, "/MS_South")
  if (!dir.exists(saveto_3)) 
    dir.create(saveto_3)
  solutions$MS_South <-  tgas_calc_OM_cond(tgas, lclass = 5, population = "ALL", 
                                           model = 2, type = 3, 
                                           use = c(TRUE, FALSE, FALSE),
                                           src = src, 
                                           dist_type = dist_type, filter_dist = filter_dist, 
                                           g_b = c(-Inf, 0),
                                           saveto = paste0(saveto_3, "/"))
  solutions$MS_South$Name <- "Main Sequence South"
  
  cat("Calc physical parameters...", "\n")
  solutions <- calc_all_physical_params(solutions)
  cat("Calc weited parameters...", "\n")
  solutions <- tgas_calc_all_weighted(solutions)
  
  return(solutions)
  
}


tgas_make_OM_solutions_dist <- function()
{
  
  tgas_make_OM_solutions(dist_type = "TGAS_PX", filter_dist = "TGAS_PX")
  tgas_make_OM_solutions(dist_type = "rMoMW", filter_dist = "TGAS_PX")
  tgas_make_OM_solutions(dist_type = "rMoExp1", filter_dist = "TGAS_PX")
  tgas_make_OM_solutions(dist_type = "rMoExp2", filter_dist = "TGAS_PX")
  
  tgas_make_OM_solutions(dist_type = "TGAS_PX", filter_dist = "rMoMW")
  tgas_make_OM_solutions(dist_type = "rMoMW", filter_dist = "rMoMW")
  tgas_make_OM_solutions(dist_type = "rMoExp1", filter_dist = "rMoMW")
  tgas_make_OM_solutions(dist_type = "rMoExp2", filter_dist = "rMoMW")
}

tgas_make_OM_solutions <- function(dist_type = "TGAS_PX", filter_dist = "TGAS_PX", src = "TGAS")
{
  tgas <- tgas_calc_LClass(tgas, dist_ = filter_dist)
  
  if (!dir.exists("solutions")) 
    dir.create("solutions")
  saveto_ <- paste0("solutions/solution_", filter_dist, "-", dist_type)
  if (!dir.exists(saveto_)) 
    dir.create(saveto_)
  
  solutions <- list();
  
  cat("Red Giants processing", "\n")
  saveto_2 <- paste0(saveto_, "/RG_ALL")
  if (!dir.exists(saveto_2)) 
    dir.create(saveto_2)
  solutions$RG_All <-  tgas_calc_OM_cond(tgas, lclass = 3, population = "ALL", type = 1, model = 1, src = src, dist_type = dist_type, filter_dist = filter_dist, saveto = paste0(saveto_2, "/"))
  solutions$RG_All$Name <- "Red Giants"
  
  cat("Main Sequence processing", "\n")
  saveto_2 <- paste0(saveto_, "/MS_ALL")
  if (!dir.exists(saveto_2)) 
    dir.create(saveto_2)
  solutions$MS_All <-  tgas_calc_OM_cond(tgas, lclass = 5, population = "ALL", type = 1, model = 1, src = src, dist_type = dist_type, filter_dist = filter_dist, saveto = paste0(saveto_2, "/"))
  solutions$MS_All$Name <- "Main Sequence"
  
  # cat("Red Giants Disk processing", "\n")
  # solutions$RG_Disk <- tgas_calc_OM_cond(tgas, lclass = 3, population = "DISK", type = 1, src = src, dist_type = dist_type, filter_dist = filter_dist, saveto = saveto_)
  # solutions$RG_Disk$Name <- "Red Giants Disk"
  # 
  # cat("Red Giants Galo processing", "\n")
  # solutions$RG_Galo <- tgas_calc_OM_cond(tgas, lclass = 3, population = "GALO", type = 1, src = src, dist_type = dist_type, filter_dist = filter_dist, saveto = saveto_)
  # solutions$RG_Galo$Name <- "Red Giants Galo"
  
  # cat("Super Giants processing", "\n")
  # solutions$SG_ALL <-  tgas_calc_OM_cond(tgas, lclass = 1, type = 1, src = src, dist_type = dist_type, filter_dist = filter_dist, saveto = saveto_)
  # solutions$SG_ALL$Name <- "Super Giants"
  # 
  # cat("Super Giants Disk processing", "\n")
  # solutions$SG_Disk <- tgas_calc_OM_cond(tgas, lclass = 1, population = "DISK", type = 1, src = src, dist_type = dist_type, filter_dist = filter_dist, saveto = saveto_)
  # solutions$SG_Disk$Name <- "Super Giants Disk"
  # 
  # cat("Super Giants Galo processing", "\n")
  # solutions$SG_Galo <- tgas_calc_OM_cond(tgas, lclass = 1, population = "GALO", type = 1, src = src, dist_type = dist_type, filter_dist = filter_dist, saveto = saveto_)
  # solutions$SG_Galo$Name <- "Super Giants Galo"
  
  cat("Calc physical parameters...", "\n")
  solutions <- calc_all_physical_params(solutions)
  cat("Calc weited parameters...", "\n")
  solutions <- tgas_calc_all_weighted(solutions)
  
  cat("Export solutions...", "\n")
  tgas_export_all_solution(solutions)
  cat("Export kinematics...", "\n")
  tgas_draw_all_kinematic(solutions, saveto = saveto_)
  cat("Export physics...", "\n")
  tgas_draw_physics(solutions, src, saveto = saveto_)
  
  #g <- tgas_draw_HR_facet(solutions$SG_ALL, M_lim = c(-2, -10), BV_lim = c(-1, 3))
  #ggsave(file = paste0(saveto_, "SG_HR_all.png"), plot = g, width = 15, height = 12)
  
  #cat("Draw RG HR facet...")
  #g <- tgas_draw_HR_facet(solutions$RG_All, M_lim = c(3, -2), BV_lim = c(0.5, 2.5))
  #ggsave(file = paste0(saveto_, "RG_ALL_HR_all.png"), plot = g, width = 15, height = 12)
  
  #cat("Draw MS HR facet...")
  #g <- tgas_draw_HR_facet(solutions$MS_All, M_lim = c(10, -3), BV_lim = c(-0.5, 2.5))
  #ggsave(file = paste0(saveto_, "MS_HR_all.png"), plot = g, width = 15, height = 12)
  
  # g <- tgas_draw_HR_facet(solutions$RG_Disk, M_lim = c(3, -2), BV_lim = c(0.5, 2.5))
  # ggsave(file = paste0(saveto_, "RG_DISK_HR_all.png"), plot = g, width = 15, height = 12)
  # 
  # g <- tgas_draw_HR_facet(solutions$RG_Galo, M_lim = c(3, -2), BV_lim = c(0.5, 2.5))
  # ggsave(file = paste0(saveto_, "RG_GALO_HR_all.png"), plot = g, width = 15, height = 12)
  
  #save(solutions, file = paste0("solution_",filter_dist, "-",dist_type,".Rdata"))
  
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

