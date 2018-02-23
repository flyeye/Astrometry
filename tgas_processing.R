

# min_px, max_px - mas
filter_tgs_px <- function(tgs,                     #  catalog (data frame)
                          px = c(-Inf,Inf),        #  limits on parallax
                          e_px = Inf,              #  limits on parallax error
                          bv_lim = c(-Inf,Inf),    #  limits on B-V
                          Mg = c(-Inf,Inf),        #  limits on absolute magnitude applied to M column in data
                          z_lim = c(0, Inf),       #  limits on distance from galactic equator
                          r_lim = c(0, Inf),       #  limits on solar distance applied to R column in data
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
  }
  tgas_ <- tgas_calc_distance(tgas_, dist_) 
  tgas_ <- tgas_calc_absolute_mag(tgas_)
  
  tgas_a <- tgas_[(!is.na(tgas_$B_V)) & (!is.na(tgas_$M)),]
  
  tgas_a <- mutate(tgas_a, LClass_apass = 0)
  tgas_a$LClass_apass[is_main_sequence(tgas_a$B_V, tgas_a$M)] <- 5
  tgas_a$LClass_apass[tgas_a$M<(-2)& (tgas_a$M>(-Inf))] <- 1
  tgas_a$LClass_apass[(tgas_a$M>-1.5)&(tgas_a$M<2.5)&(tgas_a$B_V>0.8)&(tgas_a$B_V<2.5)] <- 3
  
  tgas_ <- within(tgas_, rm("LClass_apass"))
  tgas_ <- tgas_ %>% left_join(tgas_a[ , names(tgas_a) %in% c("source_id", "LClass_apass")], by = "source_id")
  tgas_$LClass_apass[is.na(tgas_$LClass_apass)] <- 0
  
  return(tgas_)
}

min_M <- function(bv)
{
  m <- numeric(length(bv))
  
  # Верхняя граница поднята относительно Бови
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
        cat(max_dist)
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
      res_tgas$HR <- hrd;
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


tgas_calc_OM_seq <- function(tgas_ = tgas, src_ = "TGAS", start = 1, step = 0.1, q = 2, px_type = "ANGLE", distance = NULL, save = NULL, type = 0, model = 1, dist_type = "TGAS_PX", use = c(TRUE, TRUE, FALSE),...)
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
    res_tgas$HR <- hrd;
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

tgas_write_conditions <- function(conditions)
{
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
}
  

calc_physical_params <- function(solution, Rs = 8.09)
{
  
  physical <- matrix(0, nrow = nrow(solution$Oort), ncol = 8)
  colnames(physical) <- c("Vg", "P", "S", "F", "M", "apex_L", "apex_B", "Vs")
  
  e_physical <- matrix(0, nrow = nrow(solution$Oort), ncol = 8)
  colnames(e_physical) <- c("eVg", "eP", "eS", "eF", "eM", "apex_eL", "apex_eB", "eVs")
  
  physical[,1] <- Rs * (solution$Oort[,1] - solution$Oort[,2]) # Vs - линейная скорость вращения галактики на расстоянии Rs
  physical[,2] <- (2 * pi * Rs * 3.086e+16 /  physical[,1]) / (86400*365*1000000)   # P - период вращения галактики
  physical[,3] <- -(solution$Oort[,1] + solution$Oort[,2])        # S - наклон кривой вращения галактикки
  physical[,4] <- 2 * sqrt(-solution$Oort[,2]/(solution$Oort[,1] - solution$Oort[,2])) # отношение эпициклической частоты к угловой скорости вращения Галактики в окресностях Солнца
  physical[,5] <- (Rs*3.086e+19) * ((physical[,1]*1000)**2) / 132712438e+12 #6.67408e-11           # масса вещества галактики, сосредоточенная внутри орбиты Солнца
  physical[,6] <- (180/pi) * atan2(solution$X[,2], solution$X[,1])        # L
  physical[,7] <- (180/pi) * atan2(solution$X[,3], sqrt(solution$X[,1]**2 + solution$X[,2]**2)) # B
  physical[,8] <- sqrt(solution$X[,1]**2 + solution$X[,2]**2 + solution$X[,3]**2) # скоросто Солнца
  
  e_physical[, 1] <- Rs * sqrt(solution$s_Oort[,1]**2 + solution$s_Oort[,2]**2)
  e_physical[, 2] <- ((3.086e+16/(86400*365*1000000)) / (solution$Oort[,1] - solution$Oort[,2])**2) * sqrt(solution$s_Oort[,1]**2 + solution$s_Oort[,2]**2)
  e_physical[, 3] <- sqrt(solution$s_Oort[,1]**2 + solution$s_Oort[,2]**2)
  e_physical[, 4] <- sqrt(-solution$Oort[,2]/(solution$Oort[,1]-solution$Oort[,2])**3) * sqrt(solution$s_Oort[,1]**2 + (solution$s_Oort[,2]*solution$Oort[,1]/solution$Oort[,2])**2)
  
  #e_physical[, 5] <- (2 * (Rs*3.086e+19) * (physical[,1]*1000) / 132712438e+12) * e_physical[, 1]  #  6.67408e-11
  e_physical[, 5] <- (2 * (Rs*3.086e+19) * (physical[,1]*1000) / 132712438e+12) * (e_physical[, 1]*1000)  #  6.67408e-11
  
  e_physical[, 6] <- (180/pi) * (1/(solution$X[,1]**2 + solution$X[,2]**2)) * sqrt((solution$X[,2]**2) * (solution$S_X[,1]**2) + (solution$X[,1]**2) * (solution$S_X[,2]**2))
  e_physical[, 7] <- (180/pi) * (1/(solution$X[,1]**2 + solution$X[,2]**2 + solution$X[,3]**2)) * 
    sqrt((solution$X[,1]**2) * (solution$X[,3]**2) * (solution$S_X[,1]**2)/(solution$X[,1]**2 + solution$X[,2]**2) +
           (solution$X[,2]**2) * (solution$X[,3]**2) * (solution$S_X[,2]**2)/(solution$X[,1]**2 + solution$X[,2]**2) + 
           (solution$X[,1]**2 + solution$X[,2]**2)*(solution$S_X[,3]**2)/4)
  e_physical[, 8] <- 2 * sqrt(((solution$X[,1]**2) * (solution$S_X[,1]**2) + (solution$X[,2]**2) * (solution$S_X[,2]**2) + (solution$X[,3]**2) * (solution$S_X[,3]**2))/
                                (solution$X[,1]**2 + solution$X[,2]**2 + solution$X[,3]**2))
  
  
  solution$Physical <- physical
  solution$s_Physical <- e_physical
  
  return (solution)
}

calc_all_physical_params <- function(solutions, Rs = 8.09)
{
  for(i in 1:length(solutions))
  {
    solutions[[i]] <- calc_physical_params(solutions[[i]])
  }
  
  return(solutions) 
}

tgas_export_solution_xls <- function(res)
{
  s_ <- paste0(res$Conditions$SaveTo, "_", res$Conditions$Src)
  
  
  output <- cbind(res$X[,1], res$S_X[,1], res$X[,2], res$S_X[,2], res$X[,3], res$S_X[,3])
  
  for (i in 4:(ncol(res$X)))
  {
    output <- cbind(output, res$X[,i], res$S_X[,i])
  }
  
  for (i in 1:(ncol(res$Physical)))
  {
    output <- cbind(output, res$Physical[,i], res$s_Physical[,i])
  }
  
  #output <- cbind(output, res$Parameters[,1:6])
  output <- cbind(output, res$Parameters)
  
  # output <- cbind(res$X[,1], res$S_X[,1], res$X[,2], res$S_X[,2], res$X[,3], res$S_X[,3], res$X[,4], res$S_X[,4], res$X[,5], res$S_X[,5], res$X[,6], res$S_X[,6],
  #                 res$X[,7], res$S_X[,7], res$X[,8], res$S_X[,8], res$X[,9], res$S_X[,9], res$X[,10], res$S_X[,10], res$X[,11], res$S_X[,11], 
  #                 res$Physical[,1], res$s_Physical[,1],res$Physical[,2], res$s_Physical[,2],res$Physical[,3], res$s_Physical[,3],res$Physical[,4], res$s_Physical[,4],
  #                 res$Physical[,5], res$s_Physical[,5],res$Physical[,6], res$s_Physical[,6],res$Physical[,7], res$s_Physical[,7],res$Physical[,8], res$s_Physical[,8],
  #                 res$Parameters[,1:6])
  output <- t(output)
  
  z <- c(as.vector(rbind(res$wX,res$s_wX)), 
         as.vector(rbind(res$wPhysical, res$s_wPhysical)), 
         rep(x = 0, ncol(res$Parameters)))
  output <- cbind(output, as.matrix(z))
  # output <- cbind(output, as.matrix(c(res$wX[1], res$s_wX[1], res$wX[2], res$s_wX[2], res$wX[3], res$s_wX[3], res$wX[4], res$s_wX[4], 
  #                                     res$wX[5], res$s_wX[5], res$wX[6], res$s_wX[6], res$wX[7], res$s_wX[7], res$wX[8], res$s_wX[8], 
  #                                     res$wX[9], res$s_wX[9], res$wX[10], res$s_wX[10], res$wX[11], res$s_wX[11], 
  #                                     res$wPhysica[1], res$s_wPhysica[1], res$wPhysica[2], res$s_wPhysica[2], 
  #                                     res$wPhysica[3], res$s_wPhysica[3], res$wPhysica[4], res$s_wPhysica[4], 
  #                                     res$wPhysica[5], res$s_wPhysica[5], res$wPhysica[6], res$s_wPhysica[6],
  #                                     res$wPhysica[7], res$s_wPhysica[7], res$wPhysica[8], res$s_wPhysica[8], 
  #                                     rep(x = 0, 6))))
  
  z <- c(as.vector(rbind(colnames(res$X),colnames(res$S_X))), 
         as.vector(rbind(colnames(res$Physical), colnames(res$s_Physical))), 
         colnames(res$Parameters))
  
  # nc <- c(colnames(res$X)[1], colnames(res$S_X)[1], colnames(res$X)[2], colnames(res$S_X)[2], colnames(res$X)[3], colnames(res$S_X)[3],
  #                           colnames(res$X)[4], colnames(res$S_X)[4], colnames(res$X)[5], colnames(res$S_X)[5], colnames(res$X)[6], colnames(res$S_X)[6],
  #                           colnames(res$X)[7], colnames(res$S_X)[7], colnames(res$X)[8], colnames(res$S_X)[8], colnames(res$X)[9], colnames(res$S_X)[9],
  #                           colnames(res$X)[10], colnames(res$S_X)[10], colnames(res$X)[11], colnames(res$S_X)[11], 
  #                           colnames(res$Physical)[1], colnames(res$s_Physical)[1], colnames(res$Physical)[2], colnames(res$s_Physical)[2],
  #                           colnames(res$Physical)[3], colnames(res$s_Physical)[3], colnames(res$Physical)[4], colnames(res$s_Physical)[4], 
  #                           colnames(res$Physical)[5], colnames(res$s_Physical)[5], colnames(res$Physical)[6], colnames(res$s_Physical)[6], 
  #                           colnames(res$Physical)[7], colnames(res$s_Physical)[7], colnames(res$Physical)[8], colnames(res$s_Physical)[8]) 
  
  rownames(output)[1:length(z)]<- z
  
  #colnames(res$Oort)[1], colnames(res$s_Oort)[1], colnames(res$Oort)[2], colnames(res$s_Oort)[2], 
  #colnames(res$Oort)[3], colnames(res$s_Oort)[3], colnames(res$Oort)[4], colnames(res$s_Oort)[4])
  write.xlsx2(x = output, file = paste0(s_,".xls"), sheetName = "Solution")
}

tgas_export_solution_txt <- function(res)
{
  s_ <- paste0(res$Conditions$SaveTo, "_", res$Conditions$Src)
  
  output <- cbind(res$X, res$Physical)
  output_err <- cbind(res$S_X, res$s_Physical)
  
  output_txt <- paste0(sprintf("%.2f", output), "±", sprintf("%.2f", output_err))
  output_txt <- matrix(output_txt, ncol = ncol(output))
  colnames(output_txt) <- c(colnames(res$X), colnames(res$Physical))
  output_txt[,"M"] <- paste0(sprintf("%.1e", output[,"M"]), "±", sprintf("%.2e", output_err[,"eM"]))
  output_txt <- t(output_txt)
  output_txt <- rbind(output_txt, paste0(sprintf("%.2f", res$Parameters[,1]), "-", sprintf("%.2f", res$Parameters[,2]) ))
  output_txt <- rbind(output_txt, sprintf("%d", res$Parameters[,3]) )
  output_txt <- rbind(output_txt, sprintf("%.2f", res$Parameters[,4]) ) 
  output_txt <- rbind(output_txt, sprintf("%.2f", res$Parameters[,5]) ) 
  rownames(output_txt) <- c(colnames(res$X), colnames(res$Physical), c("r", "N", "r_mean", "ePx"))
  
  write_lines(stargazer(output_txt, type = "text"), paste0(s_,".txt"), append = FALSE)
}

my_san <- function (str, type = "latex", ...)
{
  return (str)
}

tgas_export_solution_latex <- function(res)
{
  
  if (res$Conditions$KinModel != 1)
    return()
  
  s_ <- paste0(res$Conditions$SaveTo, "_", res$Conditions$Src)
  
  output <- cbind(res$X[,1:5], res$X[,7:9], res$X[,6:6], res$X[,10:11], res$Physical)
  output <- cbind(res$X[,1:5], res$X[,7:9], res$X[,6:6], res$X[,10:11], res$Physical)
  output_err <- cbind(res$S_X[,1:5], res$S_X[,7:9], res$S_X[,6:6], res$S_X[,10:11], res$s_Physical)
  colnames(output)[4] <- "\\Omega_x"
  colnames(output)[5] <- "\\Omega_y"
  colnames(output)[6] <- "M_{13}"
  colnames(output)[7] <- "M_{23}"
  colnames(output)[8] <- "A"
  colnames(output)[9] <- "B"
  colnames(output)[17] <- "Apex L"
  colnames(output)[18] <- "Apex B"
  colnames(output_err)[8] <- "sA"
  colnames(output_err)[9] <- "sB"
  
  output_txt <- paste0("$", sprintf("%.2f", output), "\\pm", sprintf("%.2f", output_err), "$")
  output_txt <- matrix(output_txt, ncol = ncol(output))
  output_txt[,16] <- paste0("$", sprintf("%.1e", output[,16]), "\\pm", sprintf("%.2e", output_err[,16]), "$")
  output_txt <- t(output_txt)
  output_txt <- rbind(output_txt, paste0(sprintf("%.1f", res$Parameters[,1]), "-", sprintf("%.1f", res$Parameters[,2]) ))
  output_txt <- rbind(output_txt, sprintf("%d", res$Parameters[,3]) )
  output_txt <- rbind(output_txt, sprintf("%.2f", res$Parameters[,4]) ) 
  output_txt <- rbind(output_txt, sprintf("%.2f", res$Parameters[,5]) ) 
  rownames(output_txt) <- c(colnames(output),  c("r", "N", "\\overline{r}", "$\\overline{\\sigma}_{\\pi}$"))
  rownames(output_txt) <- paste0("$", rownames(output_txt), "$")
  
  if (ncol(output_txt)>10)
  {
    i <- ncol(output_txt) %/% 2
    #write_lines(stargazer(output_txt[,1:i], type = "latex", digits = 2), paste0(s_,"_1.tex"), append = FALSE)
    #write_lines(stargazer(output_txt[,(i+1):ncol(output_txt)], type = "latex", digits = 2), paste0(s_,"_2.tex"), append = FALSE)
    print.xtable(xtable(output_txt[,1:i]), sanitize.text.function = my_san, file = paste0(s_,"_1.tex"), append = FALSE)
    print.xtable(xtable(output_txt[,(i+1):ncol(output_txt)]), sanitize.text.function = my_san, file = paste0(s_,"_2.tex"), append = FALSE)
    
  } else 
  {
    #write_lines(stargazer(output_txt, type = "latex", digits = 2), paste0(s_,".tex"), append = FALSE) 
    print.xtable(xtable(output_txt), sanitize.text.function = my_san, file = paste0(s_,".tex"), append = FALSE)
  }
}


tgas_export_physical_latex <- function(res)
{
  s_ <- paste0(res$Conditions$SaveTo, "_", res$Conditions$Src)
  
  output <- cbind(res$Oort, res$Physical)
  output_err <- cbind(res$s_Oort, res$s_Physical)
  
  colnames(output)[1] <- "A"
  colnames(output)[2] <- "B"
  colnames(output)[3] <- "C"
  colnames(output)[4] <- "K"
  colnames(output)[12] <- "Apex L"
  colnames(output)[13] <- "Apex B"
  
  output_txt <- paste0("$", sprintf("%.2f", output), "\\pm", sprintf("%.2f", output_err), "$")
  output_txt <- matrix(output_txt, ncol = 14)
  output_txt[,11] <- paste0("$", sprintf("%.1e", output[,11]), "\\pm", sprintf("%.2e", output_err[,11]), "$")
  #output_txt <- t(matrix(output_txt, ncol = 15))
  output_txt <- t(output_txt)
  output_txt <- rbind(output_txt, paste0(sprintf("%.1f", res$Parameters[,1]), "-", sprintf("%.1f", res$Parameters[,2]) ))
  output_txt <- rbind(output_txt, sprintf("%d", res$Parameters[,3]) )
  output_txt <- rbind(output_txt, sprintf("%.2f", res$Parameters[,4]) ) 
  output_txt <- rbind(output_txt, sprintf("%.2f", res$Parameters[,5]) ) 
  #rownames(output_txt) <- c(colnames(res$X), colnames(res$Oort), c("r", "N", "r_mean", "ePx"))
  rownames(output_txt) <- c( colnames(output), 
                             c("r", "N", "\\overline{r}", "$\\overline{\\sigma}_{\\pi}$"))
  rownames(output_txt) <- paste0("$", rownames(output_txt), "$")
  
  if (ncol(output_txt)>10)
  {
    i <- ncol(output_txt) %/% 2
    #write_lines(stargazer(output_txt[,1:i], type = "latex", digits = 2), paste0(s_,"_physical_1.tex"), append = FALSE)
    #write_lines(stargazer(output_txt[,(i+1):ncol(output_txt)], type = "latex", digits = 2), paste0(s_,"_physical_2.tex"), append = FALSE)
    print.xtable(xtable(output_txt[,1:i]), sanitize.text.function = my_san, file = paste0(s_,"_physical_1.tex"), append = FALSE)
    print.xtable(xtable(output_txt[,(i+1):ncol(output_txt)]), sanitize.text.function = my_san, file = paste0(s_,"_physical_2.tex"), append = FALSE)
  } else 
  {
    #write_lines(stargazer(output_txt, type = "latex", digits = 2), paste0(s_,"_physical.tex"), append = FALSE) 
    print.xtable(xtable(output_txt), sanitize.text.function = my_san, file = paste0(s_,"_physical.tex"), append = FALSE)
  }
}

tgas_export_solution <- function(solution_)
{
  tgas_export_solution_xls(solution_)
  tgas_export_solution_txt(solution_)
  tgas_export_solution_latex(solution_)
  tgas_export_physical_latex(solution_)
}

tgas_export_all_solution <- function(solutions)
{
  for(i in 1:length(solutions))
  {
    tgas_export_solution(solutions[[i]])
  }
  
}


#------------------------------------------------------
# стоит сразу и сводные диаграммы (все параметры на одной диаграмме) для всех решений
# и отдельные диаграммы для каждого параметра (все решения для каждого параметра на одной диаграмме)
tgas_draw_all_kinematic <- function(solutions, src = "TGAS", saveto = "")
{
  for(i in 1:length(solutions))
  {
    tgas_draw_kinematic(solutions[[i]])
  }
  
  tgas_draw_all_kinematic_comp(solutions, src, saveto)
}


# tgas_draw_all_OM_sol_comp(list(solutions_mw$MS_All, solutions_mw_px$MS_All, solutions_Exp1$MS_All, solutions_Exp2$MS_All), 
# ylims  = matrix(data = c(5, 15, 10, 25, 0, 15, -5, 5, -5, 5, -15, -10, -5, 5, -3, 7 , 10, 20, -7, 3, -8, 2), nrow = 2))

# tgas_draw_all_OM_sol_comp(list(solutions_mw$RG_All, solutions_mw_px$RG_All, solutions_Exp1$RG_All, solutions_Exp2$RG_All), 
# ylims  = matrix(data = c(5, 35, 15, 55, 0, 25, -3, 7, -5, 5, -20, -10, -5, 5, -8, 2 , 7, 20, -7, 3, -10, 5), nrow = 2))

# tgas_draw_all_OM_sol_comp(solutions = solutions_bv, 
#                           ylims  = matrix(data = c(5, 20, 0, 25, 0, 15, -2, 10, -5, 2, -18, -8, -5, 5, -8, 3 , 5, 25, -7, 3, -10, 5), nrow = 2),
#                           xlims = c(-0.5, 1.1, 0.1),
#                           xpar = 9, 
#                           xtitle = "B-V", 
#                           saveto = "solutions/")


# строит отдельные диаграмм для каждого параметра в данном решении
tgas_draw_all_OM_sol_comp <- function(solutions, ylims, xlims = c(0, 2.5, 0.5), xpar = 4, xtitle = "<r>, kpc", saveto = "")
{
  
  for (i in 1:ncol(solutions[[1]]$X))
  {
    g <-draw_OMParameter(solutions, 
                         parameter = i, 
                         y_lim = c(ylims[1, i], ylims[2, i], 1), 
                         x_lim = xlims, 
                         x_title = xtitle, 
                         x_par = xpar, 
                         title = colnames(solutions[[1]]$X)[i])
    
    ggsave(paste0(saveto, "OM_", colnames(solutions[[1]]$X)[i],".png"), plot = g, width = 10, height = 5)
    ggsave(paste0(saveto, "OM_", colnames(solutions[[1]]$X)[i],".eps"), plot = g, width = 10, height = 5)
    
  }
}

#------------------------------------------------------
# строит диаграммы для всех параметров решения по двум солюшинам для сравнения

tgas_draw_all_OM_sol <- function(sol1, sol2, sol1_name, sol2_name, saveto = "")
{
  for (i in 1:ncol(sol1$X))
  {
    cname <- colnames(sol1$X)[i]
    g <- draw_OMParComp(parameter = i, sol1 = sol1, sol2 = sol2, 
                        title = cname, 
                        xat = paste(cname, "km/s/kpc,", sol1_name), 
                        yat = paste(cname, "km/s/kpc,", sol2_name))
    ggsave(paste0(saveto, "OM_", cname, "_", sol1_name, "-", sol2_name,  ".png"), plot = g, width = 5, height = 5)
    ggsave(paste0(saveto, "OM_", cname, "_", sol1_name, "-", sol2_name,  ".eps"), plot = g, width = 5, height = 5)
  }
  
}


#----------------------------------------------------

# строит три отдельные диграммы для заданного решения: все параметры Оорта, все Солнечные члены и все параметры модели ОМ

tgas_draw_kinematic <- function (solution)
{
  save <- paste0(solution$Conditions$SaveTo, "_", solution$Conditions$Src,"_")
  
  g <- draw_OM_Solar(solution, title = paste("Solar motion,", solution$Conditions$Src, " proper motions, Solar motion."))
  ggsave(paste0(save, "Solar-R", ".png"), plot = g, width = 10, height = 10)
  ggsave(paste0(save, "Solar-R", ".eps"), plot = g, width = 10, height = 10)
  
  if (solution$Conditions$KinModel <=2 ){
    g <- draw_Oort(solution, title = paste("Oort`s parameters,", solution$Conditions$Src, " proper motions."))
    ggsave(paste0(save, "OL-R", ".png"), plot = g, width = 10, height = 10)
    ggsave(paste0(save, "OL-R", ".eps"), plot = g, width = 10, height = 10)
  }
  
  if (solution$Conditions$KinModel == 1){
    g <- draw_OM(solution, title = paste("Ogorodnikov-Miln Model,", solution$Conditions$Src, " proper motions."))
    ggsave(paste0(save, "OM-R", ".png"), plot = g, width = 10, height = 10)
    ggsave(paste0(save, "OM-R", ".eps"), plot = g, width = 10, height = 10)
  }
  
}

#----------------------------------------------------

## Строит отдельные диаграммы для каждого параметра со всем решениями на каждой диаграмме

## tgas_draw_all_kinematic_comp(solutions_bv, src = "TGAS", saveto = "solutions/bottlinger_kin_", x_par = 9, x_lim = c(-0.5, 1.2,0.2), x_title = "B-V", is_legend = FALSE, width = 4, height = 4)


tgas_draw_all_kinematic_comp <- function(solutions, src = "TGAS", saveto = "", 
                                         x_par = 4, x_lim = c(0, 4, 0.5), x_title = "<r>, kpc",
                                         is_legend = TRUE, 
                                         width = 10, height = 5)
{
  save <- paste0(saveto, src,"_")
  
  g <- draw_OMParameter(solutions,
                        parameter = 1,
                        x_par = x_par,
                        x_lim = x_lim, 
                        x_title = x_title,
                        y_lim = c(5, 20, 5),
                        y_title = "Компонент движения Солнца U, км/с",
                        is_legend = is_legend,
                        title = paste("Solar motion U, ", src, " proper motions."))
  ggsave(paste0(save, "SolarU-R", ".png"), plot = g, width = width, height = height)
  #ggsave(paste0(save, "SolarU-R", ".eps"), plot = g, width = width, height = height)
  
  g <- draw_OMParameter(solutions,
                        parameter = 2,
                        x_par = x_par,
                        x_lim = x_lim, 
                        x_title = x_title,
                        y_lim = c(5, 35, 5),
                        y_title = "Компонент движения Солнца V, км/с",
                        is_legend = is_legend,
                        title = paste("Компонент движения Солнца V, ", src, " proper motions."))
  ggsave(paste0(save, "SolarV-R", ".png"), plot = g, width = width, height = height)
  #ggsave(paste0(save, "SolarV-R", ".eps"), plot = g, width = width, height = height)
  
  g <- draw_OMParameter(solutions,
                        parameter = 3,
                        x_par = x_par,
                        x_lim = x_lim,
                        x_title = x_title,
                        y_lim = c(5, 20, 2.5),
                        y_title = "Компонент движения Солнца W, км/с",
                        is_legend = is_legend,
                        title = paste("Solar motion W, ", src, " proper motions."))
  ggsave(paste0(save, "SolarW-R", ".png"), plot = g, width = width, height = height)
  #ggsave(paste0(save, "SolarW-R", ".eps"), plot = g, width = width, height = height)
  
  g <- draw_OortParameter(solutions, 
                          parameter = 1,
                          x_lim = x_lim, 
                          x_par = x_par,
                          x_title = x_title,
                          #data_x_lim = c(1, 1),
                          is_legend = is_legend,
                          y_lim = c(0, 40, 5), 
                          y_title = "Параметр Оорта A, км/с/кпк",
                          title = paste("Oort`s parameter A, ", src, " proper motions."))
  ggsave(paste0(save, "OortA-R", ".png"), plot = g, width = width, height = height)
  #ggsave(paste0(save, "OortA-R", ".eps"), plot = g, width = width, height = height)
  
  
  g <- draw_OortParameter(solutions, 
                          parameter = 2,
                          x_lim = x_lim, 
                          x_par = x_par,
                          x_title = x_title,
                          #data_x_lim = c(1, 14),
                          is_legend = is_legend,
                          title = paste("Oort`s parameter B, ", src, " proper motions."),
                          y_lim = c(-20, -8, 2), 
                          y_title = "Параметр Оорта B, км/с/кпк")
  ggsave(paste0(save, "OortB-R", ".png"), plot = g, width = width, height = height)
  #ggsave(paste0(save, "OortB-R", ".eps"), plot = g, width = width, height = height)
  
  if ((solutions[[1]]$Conditions$KinModel==1) |
      ((solutions[[1]]$Conditions$KinModel==2) & (solutions[[1]]$Conditions$KinModelType>=1)) |
      (solutions[[1]]$Conditions$KinModel==4)
      )
  {
    g <- draw_OortParameter(solutions, 
                            parameter = 3,
                            x_par = x_par,
                            x_lim = x_lim, 
                            x_title = x_title,
                            #data_x_lim = c(2, 13),
                            is_legend = is_legend,
                            title = paste("Oort`s parameter C, ", src, " proper motions."),
                            y_lim = c(-15, 3, 3), 
                            y_title = "Параметр Оорта С, км/с/кпк")
    ggsave(paste0(save, "OortC-R", ".png"), plot = g, width = width, height = height)
    #ggsave(paste0(save, "OortC-R", ".eps"), plot = g, width = width, height = height)
    
    g <- draw_OortParameter(solutions, 
                            parameter = 4,
                            x_par = x_par,
                            x_lim = x_lim, 
                            x_title = x_title,
                            #data_x_lim = c(3, 13),
                            is_legend = is_legend,
                            title = paste("Oort`s parameter K, ", src, " proper motions."),
                            y_lim = c(-15, 3, 3), 
                            y_title = "Параметр Оорта K, км/с/кпк")
    ggsave(paste0(save, "OortK-R", ".png"), plot = g, width = width, height = height)
    #ggsave(paste0(save, "OortK-R", ".eps"), plot = g, width = width, height = height)  
  }
  
  if (solutions[[1]]$Conditions$KinModel==1)
  {
    g <- draw_OMParameter(solutions,
                          parameter = 4,
                          x_par = x_par,
                          x_lim = x_lim, 
                          x_title = x_title,
                          #data_x_lim = c(3, 13),
                          y_lim = c(-6, 8, 2),
                          y_title = expression(omega[x]*", км/с/кпк"), 
                          is_legend = is_legend,
                          title = paste("Ogorodnikov-Miln Wx, ", src, " proper motions."))
    ggsave(paste0(save, "Wx", ".png"), plot = g, width = width, height = height)
    #ggsave(paste0(save, "Wx", ".eps"), plot = g, width = width, height = height)
    
    g <- draw_OMParameter(solutions,
                          parameter = 5,
                          x_par = x_par,
                          x_lim = x_lim, 
                          x_title = x_title,
                          #data_x_lim = c(3, 13),
                          y_lim = c(-6, 8, 2),
                          y_title = expression(omega[y]*", км/с/кпк"), 
                          is_legend = is_legend,
                          title = paste("Ogorodnikov-Miln Wy, ", src, " proper motions."))
    ggsave(paste0(save, "Wy", ".png"), plot = g, width = width, height = height)
    #ggsave(paste0(save, "Wy", ".eps"), plot = g, width = width, height = height)
    
    g <- draw_OMParameter(solutions,
                          parameter = 7,
                          x_par = x_par,
                          x_lim = x_lim, 
                          x_title = x_title,
                          #data_x_lim = c(3, 13),
                          y_lim = c(-8, 4, 2),
                          y_title = expression(M[13]*", км/с/кпк"), 
                          is_legend = is_legend,
                          title = paste("Ogorodnikov-Miln M13, ", src, " proper motions."))
    ggsave(paste0(save, "M13", ".png"), plot = g, width = width, height = height)
    #ggsave(paste0(save, "M13", ".eps"), plot = g, width = width, height = height)
    
    g <- draw_OMParameter(solutions,
                          parameter = 8,
                          x_par = x_par,
                          x_lim = x_lim, 
                          x_title = x_title,
                          #data_x_lim = c(3, 13),
                          y_lim = c(-8, 4, 2),
                          y_title = expression(M[23]*", км/с/кпк"), 
                          is_legend = is_legend,
                          title = paste("Ogorodnikov-Miln M23, ", src, " proper motions."))
    ggsave(paste0(save, "M23", ".png"), plot = g, width = width, height = height)
    #ggsave(paste0(save, "M23", ".eps"), plot = g, width = width, height = height)
    
  }
  
  if ((solutions[[1]]$Conditions$KinModel==2) & (solutions[[1]]$Conditions$KinModelType==2))
  {
    g <- draw_OortParameter(solutions, 
                            parameter = 5,
                            x_par = x_par,
                            x_lim = x_lim, 
                            x_title = x_title,
                            is_legend = is_legend,
                            title = paste("Gx, ", src, " proper motions."),
                            y_title = "Gx, km/s/kpc", 
                            y_lim = c(-40, 40, 2))
    ggsave(paste0(save, "OortGx-R", ".png"), plot = g, width = width, height = height)
    #ggsave(paste0(save, "OortGx-R", ".eps"), plot = g, width = width, height = height)
    
    g <- draw_OortParameter(solutions, 
                            parameter = 6,
                            x_par = x_par,
                            x_lim = x_lim, 
                            x_title = x_title,
                            is_legend = is_legend,
                            title = paste("Gy, ", src, " proper motions."),
                            y_title = "Gy, km/s/kpc", 
                            y_lim = c(-40, 40, 3))
    ggsave(paste0(save, "OortGy-R", ".png"), plot = g, width = width, height = height)
    #ggsave(paste0(save, "OortGy-R", ".eps"), plot = g, width = width, height = height)  
  }
  
  if (solutions[[1]]$Conditions$KinModel==4) 
  {
    g <- draw_OMParameter(solutions, 
                            parameter = 4,
                            x_par = x_par,
                            x_lim = x_lim, 
                            x_title = x_title,
                            #data_x_lim = c(1, 14),
                            is_legend = is_legend,
                            title = paste("Bottlinger`s parameter W, ", src, " proper motions."),
                            y_title = expression(Omega[0]*", км/с/кпк"), 
                          y_lim = c(17.5, 35, 2.5))
    ggsave(paste0(save, "Bottlinger_W", ".png"), plot = g, width = width, height = height)
    #ggsave(paste0(save, "Bottlinger_W", ".eps"), plot = g, width = width, height = height)
    
    g <- draw_OMParameter(solutions, 
                            parameter = 5,
                            x_par = x_par,
                            x_lim = x_lim, 
                            x_title = x_title,
                            #data_x_lim = c(1, 16),
                            is_legend = is_legend,
                            title = paste("Bottlinger`s parameter W', ", src, " proper motions."),
                            y_title = expression(Omega[0]*"', км/с/кпк"^2), 
                            y_lim = c(-15, 0, 3))
    ggsave(paste0(save, "Bottlinger_W1", ".png"), plot = g, width = width, height = height)
    #ggsave(paste0(save, "Bottlinger_W1", ".eps"), plot = g, width = width, height = height)
    
    g <- draw_OMParameter(solutions, 
                            parameter = 6,
                            x_par = x_par,
                            x_lim = x_lim, 
                            x_title = x_title,
                            #data_x_lim = c(1, 16),
                            is_legend = is_legend,
                            title = paste("Bottlinger`s parameter W\", ", src, " proper motions."),
                            y_title = expression(Omega[0]*"'', км/с/кпк"^3), 
                            y_lim = c(-18, 6, 3))
    ggsave(paste0(save, "Bottlinger_W2", ".png"), plot = g, width = width, height = height)
    #ggsave(paste0(save, "Bottlinger_W2", ".eps"), plot = g, width = width, height = height)  
  }
}


#----------------------------------------------------


tgas_draw_physics <- function (solution, src = "TGAS", saveto = "", is_legend = TRUE, width = 10, height = 5)
{
  
  if (solution[[1]]$Conditions$KinModel == 3)
    return(0);
  
  save <- paste0(saveto, src,"_")
  
  g <- draw_Physical(solution, 
                     x_par = 4,
                     is_legend = is_legend,
                     title = paste("Linear galactic velocity at Solar distance, ", src, " proper motions."))
  ggsave(paste0(save, "V-R", ".png"), plot = g, width = width, height = height)
  ggsave(paste0(save, "V-R", ".eps"), plot = g, width = width, height = height)
  
  g <- draw_Physical(solution, 
                     parameter = 2,
                     title = paste("Galaxy rotation period, ", src, " proper motions."),
                     x_par = 4,
                     y_lim = c(205, 300, 10),
                     is_legend = is_legend,
                     y_title = "million years")
  ggsave(paste0(save, "Period-R", ".png"), plot = g, width = width, height = height)
  ggsave(paste0(save, "Period-R", ".eps"), plot = g, width = width, height = height)
  
  #draw_GalRotationCurveTilt(solution)
  g <- draw_Physical(solution, 
                     parameter = 3,
                     title = paste("Galaxy rotation curve inclination, ", src, " proper motions."),
                     x_par = 4,
                     y_lim = c(-7, 7, 1),
                     is_legend = is_legend,
                     y_title = "km/s/kpc")
  ggsave(paste0(save, "S-R", ".png"), plot = g, width = width, height = height)
  ggsave(paste0(save, "S-R", ".eps"), plot = g, width = width, height = height)
  
  #draw_GalF(solution)
  g <- draw_Physical(solution, 
                     parameter = 4,
                     title = paste("Epicyclic frequency to angular velocity, ", src, " proper motions."),
                     x_par = 4,
                     y_lim = c(1.2, 1.6, 0.1),
                     is_legend = is_legend,
                     y_title = "")
  ggsave(paste0(save, "F-R", ".png"), plot = g, width = width, height = height)
  ggsave(paste0(save, "F-R", ".eps"), plot = g, width = width, height = height)
  
  #draw_GalMass(solution)
  g <- draw_Physical(solution, 
                     parameter = 5,
                     title = paste("Galaxy mass inside Solar orbit, ", src, " proper motions."),
                     x_par = 4,
                     y_lim = c(5.0e10, 11e10, 1e10),
                     is_legend = is_legend,
                     y_title = "Solar mass")
  ggsave(paste0(save, "M-R", ".png"), plot = g, width = width, height = height)
  ggsave(paste0(save, "M-R", ".eps"), plot = g, width = width, height = height)
  
  #draw_ApexL(solution)
  g <- draw_Physical(solution, 
                     parameter = 6,
                     title = paste("Solar motion apex L, ", src, " proper motions."),
                     x_par = 4,
                     y_lim = c(55, 74, 1),
                     is_legend = is_legend,
                     y_title = "degree")
  ggsave(paste0(save, "ApexL-R", ".png"), plot = g, width = width, height = height)
  ggsave(paste0(save, "ApexL-R", ".eps"), plot = g, width = width, height = height)
  
  #draw_ApexB(solution)
  g <- draw_Physical(solution, 
                     parameter = 7,
                     title = paste("Solar motion apex B, ", src, " proper motions."),
                     x_par = 4,
                     y_lim = c(11, 23, 1),
                     is_legend = is_legend,
                     y_title = "degree")
  ggsave(paste0(save, "ApexB-R", ".png"), plot = g, width = width, height = height)
  ggsave(paste0(save, "ApexB-R", ".eps"), plot = g, width = width, height = height)
  
  #draw_SolarV(solution)
  g <- draw_Physical(solution, 
                     parameter = 8,
                     title = paste("Solar velocity, ", src, " proper motions."),
                     x_par = 4,
                     y_lim = c(10, 100, 10),
                     is_legend = is_legend,
                     y_title = "km/s")
  ggsave(paste0(save, "SolarV-R", ".png"), plot = g, width = width, height = height)
  ggsave(paste0(save, "SolarV-R", ".eps"), plot = g, width = width, height = height)
}


## tgas_draw_physics_BV(solutions_bv, src = "TGAS", saveto = "solutions/bottlinger_physycs_", x_lim = c(-0.5, 1.2, 0.2), is_legend = FALSE, width = 4, height = 4)

tgas_draw_physics_BV <- function (solution, src = "TGAS", saveto = "", x_lim = c(-1, 2, 0.1), is_legend = TRUE, width = 10, height = 5)
{
  
  if (solution[[1]]$Conditions$KinModel == 3)
    return(0);
  
  save <- paste0(saveto, src,"_")
  
  g <- draw_Physical(solution, 
                     x_par = 9, 
                     x_title = "B-V",
                     x_lim = x_lim,
                     data_x_lim = c(1, 14),
                     y_lim = c(120, 300, 25),
                     is_legend = is_legend,
                     y_title = "Линейная скорость Солнца V, км/с",
                     title = paste("Linear galactic velocity at Solar distance, ", src, " proper motions."))
  ggsave(paste0(save, "V-R", ".png"), plot = g, width = width, height = height)
  #ggsave(paste0(save, "V-R", ".eps"), plot = g, width = width, height = height)
  
  g <- draw_Physical(solution, 
                     parameter = 2,
                     title = paste("Период вращения Галактики, ", src, " proper motions."),
                     x_par = 9, 
                     x_title = "B-V",
                     x_lim = x_lim,
                     data_x_lim = c(1, 14),
                     is_legend = is_legend,
                     y_lim = c(175, 300, 25),
                     y_title = "Период вращения Галактики, млн.лет")
  ggsave(paste0(save, "Period-R", ".png"), plot = g, width = width, height = height)
  #ggsave(paste0(save, "Period-R", ".eps"), plot = g, width = width, height = height)
  
  #draw_GalRotationCurveTilt(solution)
  g <- draw_Physical(solution, 
                     parameter = 3,
                     title = paste("Galaxy rotation curve inclination, ", src, " proper motions."),
                     x_par = 9, 
                     x_title = "B-V",
                     data_x_lim = c(2, 14),
                     x_lim = x_lim,
                     is_legend = is_legend,
                     y_lim = c(-12, 6, 2),
                     y_title = "Наклон кривой вращения Галактики, км/с/кпк")
  ggsave(paste0(save, "S-R", ".png"), plot = g, width = width, height = height)
  #ggsave(paste0(save, "S-R", ".eps"), plot = g, width = width, height = height)
  
  #draw_GalF(solution)
  g <- draw_Physical(solution, 
                     parameter = 4,
                     title = paste("Epicyclic frequency to angular velocity, ", src, " proper motions."),
                     x_par = 9, 
                     data_x_lim = c(1, 14),
                     x_title = "B-V",
                     x_lim = x_lim,
                     is_legend = is_legend,
                     y_lim = c(1.0, 1.8, 0.2),
                     y_title = "F")
  ggsave(paste0(save, "F-R", ".png"), plot = g, width = width, height = height)
  #ggsave(paste0(save, "F-R", ".eps"), plot = g, width = width, height = height)
  
  #draw_GalMass(solution)
  g <- draw_Physical(solution, 
                     parameter = 5,
                     title = paste("Galaxy mass inside Solar orbit, ", src, " proper motions."),
                     x_par = 9, 
                     x_title = "B-V",
                     data_x_lim = c(2, 13),
                     x_lim = x_lim,
                     is_legend = is_legend,
                     y_lim = c(5.0e10, 1.50e11, 2.5e10),
                     y_title = "масса Галактики, масса Солнца")
  ggsave(paste0(save, "M-R", ".png"), plot = g, width = width, height = height)
  #ggsave(paste0(save, "M-R", ".eps"), plot = g, width = width, height = height)
  
  #draw_ApexL(solution)
  g <- draw_Physical(solution, 
                     parameter = 6,
                     title = paste("Solar motion apex L, ", src, " proper motions."),
                     x_par = 9, 
                     x_title = "B-V",
                     #data_x_lim = c(1, 19),
                     x_lim = x_lim,
                     is_legend = is_legend,
                     y_lim = c(20, 102, 10),
                     y_title = "Долгота апекса Солнца L, градусы")
  ggsave(paste0(save, "ApexL-R", ".png"), plot = g, width = width, height = height)
  #ggsave(paste0(save, "ApexL-R", ".eps"), plot = g, width = width, height = height)
  
  #draw_ApexB(solution)
  g <- draw_Physical(solution, 
                     parameter = 7,
                     title = paste("Solar motion apex B, ", src, " proper motions."),
                     x_par = 9, 
                     x_title = "B-V",
                     x_lim = x_lim,
                     is_legend = is_legend,
                     y_lim = c(10, 40, 5),
                     y_title = "Широта апекса Солнца B, градусы")
  ggsave(paste0(save, "ApexB-R", ".png"), plot = g, width = width, height = height)
  #ggsave(paste0(save, "ApexB-R", ".eps"), plot = g, width = width, height = height)
  
  #draw_SolarV(solution)
  g <- draw_Physical(solution, 
                     parameter = 8,
                     title = paste("Solar velocity, ", src, " proper motions."),
                     x_par = 9, 
                     x_title = "B-V",
                     x_lim = x_lim,
                     is_legend = is_legend,
                     y_lim = c(5, 30, 5),
                     y_title = "скорость Cолнца относительно центроида, км/с")
  ggsave(paste0(save, "SolarV-R", ".png"), plot = g, width = width, height = height)
  #ggsave(paste0(save, "SolarV-R", ".eps"), plot = g, width = width, height = height)
}


#------------------------------------------------------

tgas_draw_HR_facet <- function(solution, M_lim = c(10,-10), BV_lim = c(-1, 3))
{
  gl <- list()
  for (i in 1:length(solution$SolutionR))
  {
    g <- solution$SolutionR[[i]]$HR
    g <- g + ggtitle(paste(solution$Parameters[i,1], "-",paste(solution$Parameters[i,2]," kpc"))) + 
      scale_y_reverse(breaks=seq(M_lim[1],M_lim[2],by=-1), minor_breaks=seq(M_lim[1],M_lim[2],by=-0.5), limits = M_lim) +
      scale_x_continuous(breaks=seq(BV_lim[1],BV_lim[2],by=0.5), minor_breaks=seq(BV_lim[1],BV_lim[2],by=0.25), limits = BV_lim)
    gl <- append(gl, list(g))
  }
  
  g <- grid.arrange(grobs = gl, ncol = 4)
  
  return (g)
}

#------------------------------------------------------

tgas_calc_weighted_mean <- function(X, sX)
{
  res <- list()
  a <- matrix(nrow = nrow(sX), ncol= ncol(sX))
  p <- numeric(0)
  z <- numeric(0)
  s0_1 <- numeric(0)
  s0_2 <- numeric(0)  
  for(i in 1:ncol(sX))
  {
    a[,i] <- min(sX[,i])**2/(sX[,i])**2
    p[i] <- sum(a[,i])
    z[i] <- sum(X[,i]*a[,i])/p[i]
    s0_1[i] <- sqrt(1/sum(1/sX[,i]**2))
    s0_2[i] <- sqrt(sum(a[,i]*(X[,i]-z[i])**2)/(nrow(X)-1))/sqrt(p[i])
  }
  
  res$wX <- z
  res$s_wX <- (s0_1+s0_2)/2
  
  
  return(res)
}


tgas_calc_solution_weighted <- function(solution_)
{
  wr <- tgas_calc_weighted_mean(solution_$X, solution_$S_X)
  solution_$wX <- wr$wX
  solution_$s_wX <- wr$s_wX
  wr <- tgas_calc_weighted_mean(solution_$Physical, solution_$s_Physical)
  solution_$wPhysical <- wr$wX
  solution_$s_wPhysical <- wr$s_wX
  return(solution_)
}


tgas_calc_all_weighted <- function(solutions)
{
  for(i in 1:length(solutions))
  {
    solutions[[i]] <- tgas_calc_solution_weighted(solutions[[i]]) 
    
    # wr <- tgas_calc_weighted_mean(solutions[[i]]$X, solutions[[i]]$S_X)
    # solutions[[i]]$wX <- wr$wX
    # solutions[[i]]$s_wX <- wr$s_wX
    # wr <- tgas_calc_weighted_mean(solutions[[i]]$Physical, solutions[[i]]$s_Physical)
    # solutions[[i]]$wPhysical <- wr$wX
    # solutions[[i]]$s_wPhysical <- wr$s_wX
  }
  return(solutions)
}

tgas_process_solution <- function(solution_)
{
  solution_ <- calc_physical_params(solution_) 
  solution_ <- tgas_calc_solution_weighted(solution_)
  tgas_export_solution_xls(solution_)
  return(solution_)
}


tgas_export_all <- function(solutions, saveto_)
{
  
  cat("Export solutions...", "\n")
  tgas_export_all_solution(solutions)
  cat("Export kinematics...", "\n")
  tgas_draw_all_kinematic(solutions, saveto_)
  cat("Export physics...", "\n")
  tgas_draw_physics(solutions, src, saveto_)
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


# -----------------------------------------------------------------------------
#                                    Diagrams
# -----------------------------------------------------------------------------

draw_HR_TGAS_T2Sp <- function()
{
  HRDiagram(tgas, title = "Hertzsprung-Russell")
  ggsave("HR-TGAS-T2sp.png", width = 10, height = 10)
  HRDiagram(filter(tgas, LClass == 1), title = "Hertzsprung-Russell L1 class")
  ggsave("HR-TGAS-T2sp-L1.png", width = 10, height = 10)
  HRDiagram(filter(tgas, LClass == 2), title = "Hertzsprung-Russell L2 class")
  ggsave("HR-TGAS-T2sp-L2.png", width = 10, height = 10)
  HRDiagram(filter(tgas, LClass == 3), title = "Hertzsprung-Russell L3 class")
  ggsave("HR-TGAS-T2sp-L3.png", width = 10, height = 10)
  HRDiagram(filter(tgas, LClass == 4), title = "Hertzsprung-Russell L4 class")
  ggsave("HR-TGAS-T2sp-L4.png", width = 10, height = 10)
  HRDiagram(filter(tgas, LClass == 5), title = "Hertzsprung-Russell L5 class")
  ggsave("HR-TGAS-T2sp-L5.png", width = 10, height = 10)
}

