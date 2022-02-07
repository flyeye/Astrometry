### =======================================================
### -------------------------------------------------------
###  по шаровым слоям
GDR3_calc_OM_cond <- function(data, lclass = 3, population = "ALL", 
                              type = 0, model = 1, 
                              use = c(TRUE, TRUE, FALSE), dist_type = "GAIA_PX", filter_dist = "GAIA_PX", 
                              saveto = "", 
                              g_b = c(-Inf, Inf))
{
  
  #APASS photometry
  
  conditions <- list();
  conditions$Time <- Sys.time();
  conditions$Dist_Type <- dist_type;
  conditions$Filter_Dist <- filter_dist;
  #data <- GDR3_calc_distance(data, conditions$Filter_Dist)
  conditions$use <- use;
  conditions$KinModel <- model
  conditions$KinModelType <- type
  conditions$g_B <- g_b
  conditions$SaveTo <- saveto
  
  conditions$BV <- c(-Inf, Inf)
  conditions$MG <- c(-Inf, Inf)
  conditions$e_Px <- Inf
  
  conditions$Population <- population;
  if (population == "DISK")
  {
    conditions$Z <- c(0, 0.5)
  } else if (population == "GALO")
  {
    conditions$Z <- c(0.5, Inf)
  } else 
  {
    conditions$Z <- c(0, Inf)
  }
  
  conditions$LClass <- lclass;
  
  if(lclass == 1)
  {
    
    if (population == "DISK")
    {
      distance_ <- cbind(seq(0, 8, by = 2), seq(2, 8, by = 2))
      
    } else if (population == "GALO")
    {
      
      distance_ <- cbind(seq(0, 8, by = 2), seq(2, 8, by = 2))
      
    } else 
    {
      distance_ <- cbind(seq(0, 8, by = 2), seq(2, 8, by = 2))
    }
    
  } else if(lclass == 3)
  {
    
    if (population == "DISK")
    {
      
      distance_ <- cbind(seq(0.1, 2.9, by = 0.1), seq(0.2, 3.0, by = 0.1))
      
    } else if (population == "GALO")
    {
      distance_ <- cbind(seq(0.1, 2.9, by = 0.1), seq(0.2, 3.0, by = 0.1))
    
    }else 
    {
      distance_ <- cbind(seq(0.1, 2.9, by = 0.1), seq(0.2, 3.0, by = 0.1))
    }
    
  } else if(lclass == 5)
  {
    
    conditions$e_Px <- Inf
    
    distance_ <- cbind(seq(0.1, 2.9, by = 0.1), seq(0.2, 3.0, by = 0.1))
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
  
  
  solution  <- GDR3_calc_OM_seq_2(data, px_type = "DIST", distance = distance_, 
                                  save = conditions$SaveTo,
                                  type = conditions$KinModelType, model = conditions$KinModel,
                                  dist_type = dist_type, use = conditions$use,
                                  z_lim = conditions$Z, e_px = conditions$e_Px, bv = conditions$BV, Mg = conditions$MG,
                                  g_b = conditions$g_B, lclass = conditions$LClass)

  solution$Conditions <- conditions
  
  return(solution)
  
}


#============================================================================
# модель ОМ, по шаровым слоям, отдельно для ГП и КГ, и отдельно для разного звездного населения Галактики
GDR3_make_OM_solutions_R <- function(data, filter_dist = "GAIA_PX", name = "R", ph = "APASS", path = "/Catalogues/Gaia/")
{
  
  # тут надо бы сначала применить фотометрию
  #data <- GDR3_calc_LClass(data, dist_ = filter_dist)
  
  if (!dir.exists(paste0(path, "Solutions")))"/Catalogues/Gaia/"
    dir.create(paste0(path, "Solutions"))
  saveto_ <- paste0(path, "Solutions/solution_",name,"_",Sys.Date(), "_", filter_dist, "_", ph)
  if (!dir.exists(saveto_)) 
    dir.create(saveto_)
  
  solutions <- list();
  
  cat("Red Giants processing", "\n")
  saveto_2 <- paste0(saveto_, "/RG_ALL")
  if (!dir.exists(saveto_2))
    dir.create(saveto_2)
  solutions$RG_All <-  GDR3_calc_OM_cond(data, lclass = 3, population = "ALL", type = 1, model = 1, dist_type = "RGEO", filter_dist = filter_dist, saveto = paste0(saveto_2, "/"))
  solutions$RG_All$Name <- "Red Giants"

  cat("Main Sequence processing", "\n")
  saveto_2 <- paste0(saveto_, "/MS_ALL")
  if (!dir.exists(saveto_2))
    dir.create(saveto_2)
  solutions$MS_All <-  GDR3_calc_OM_cond(data, lclass = 5, population = "ALL", type = 1, model = 1, dist_type = "RGEO", filter_dist = filter_dist, saveto = paste0(saveto_2, "/"))
  solutions$MS_All$Name <- "Main Sequence"
  
  cat("Red Giants Disk processing", "\n")
  saveto_2 <- paste0(saveto_, "/RG_DISK")
  if (!dir.exists(saveto_2)) 
    dir.create(saveto_2)
  solutions$RG_Disk <- GDR3_calc_OM_cond(data, lclass = 3, population = "DISK", type = 1, model = 1, dist_type = "RGEO", filter_dist = filter_dist, saveto = paste0(saveto_2, "/"))
  solutions$RG_Disk$Name <- "Red Giants Disk"
   
  cat("Red Giants Galo processing", "\n")
  saveto_2 <- paste0(saveto_, "/RG_GALO")
  if (!dir.exists(saveto_2)) 
    dir.create(saveto_2)
  solutions$RG_Galo <- GDR3_calc_OM_cond(data, lclass = 3, population = "GALO", type = 1, model = 1, dist_type = "RGEO", filter_dist = filter_dist, saveto = paste0(saveto_2, "/"))
  solutions$RG_Galo$Name <- "Red Giants Galo"
  
  cat("Main Sequence Disk processing", "\n")
  saveto_2 <- paste0(saveto_, "/MS_DISK")
  if (!dir.exists(saveto_2)) 
    dir.create(saveto_2)
  solutions$MS_Disk <-  GDR3_calc_OM_cond(data, lclass = 5, population = "DISK", type = 1, model = 1, dist_type = "RGEO", filter_dist = filter_dist, saveto = paste0(saveto_2, "/"))
  solutions$MS_Disk$Name <- "Main Sequence Disk"
  
  cat("Main Sequence Galo processing", "\n")
  saveto_2 <- paste0(saveto_, "/MS_GALO")
  if (!dir.exists(saveto_2)) 
    dir.create(saveto_2)
  solutions$MS_Galo <-  GDR3_calc_OM_cond(data, lclass = 5, population = "Galo", type = 1, model = 1, dist_type = "RGEO", filter_dist = filter_dist, saveto = paste0(saveto_2, "/"))
  solutions$MS_Galo$Name <- "Main Sequence Galo"
  
  #cat("Super Giants processing", "\n")
  #saveto_2 <- paste0(saveto_, "/SG")
  #if (!dir.exists(saveto_2)) 
  #  dir.create(saveto_2)
  #solutions$SG_ALL <-  GDR3_calc_OM_cond(data, lclass = 1, type = 1, model = 1, dist_type = "RGEO", filter_dist = filter_dist, saveto = paste0(saveto_2, "/"))
  #solutions$SG_ALL$Name <- "Super Giants"
  
    # cat("Super Giants Disk processing", "\n")
  # solutions$SG_Disk <- tgas_calc_OM_cond(data, lclass = 1, population = "DISK", type = 1, src = src, dist_type = dist_type, filter_dist = filter_dist, saveto = saveto_)
  # solutions$SG_Disk$Name <- "Super Giants Disk"
  # 
  # cat("Super Giants Galo processing", "\n")
  # solutions$SG_Galo <- tgas_calc_OM_cond(data, lclass = 1, population = "GALO", type = 1, src = src, dist_type = dist_type, filter_dist = filter_dist, saveto = saveto_)
  # solutions$SG_Galo$Name <- "Super Giants Galo"
  
  cat("Calc physical parameters...", "\n")
  solutions <- calc_all_physical_params(solutions)
  cat("Calc weited parameters...", "\n")
  solutions <- calc_all_weighted(solutions)
  
  #s <<-solutions
  
  cat("Export solutions...", "\n")
  export_all_solution(solutions)
  cat("Draw kinematics...", "\n")
  draw_all_kinematic(solutions, saveto = saveto_)
  cat("Draw physics...", "\n")
  draw_physics(solutions, saveto = saveto_)
  
  cat("Draw O-M...", "\n")
  draw_all_OM_sol_comp(solutions = solutions, 
                       ylims  = matrix(data = c(0, 20, 1,
                                                0, 35, 1, 
                                                0, 15, 1, 
                                                -2, 10, 1, 
                                                -5, 2, 1, 
                                                -25, -5, 1, 
                                                -5, 5, 1, 
                                                -8, 5, 1, 
                                                5, 30, 1, 
                                                -10, 3, 1, 
                                                -10, 5, 1), nrow = 3),
                       xlims = c(0, 3, 0.2),
                       xpar = 4, 
                       xtitle = "R", 
                       saveto = paste0(saveto_, "/"))
  
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


#================================================================
# расширенная модель Оорта, отдельно для звезд ГП и КГ, по сфере и по полусферам, шаровыми слоями
GDR3_make_Oort_solutions <- function(data, dist_type = "GAIA_PX", filter_dist = "GAIA_PX", name = "", path = "/Catalogues/Gaia/")
{
  
  if (!dir.exists(paste0(path, "Solutions")))
    dir.create(paste0(path, "Solutions"))
  
  saveto_ <- paste0(path, "Solutions/solution_",name,"_",Sys.Date(), "_", dist_type)
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
  solutions$RG_All <-  GDR3_calc_OM_cond(data, lclass = 3, population = "ALL",
                                         model = 2, type = 2,
                                         use = c(TRUE, TRUE, FALSE),
                                         g_b = c(-Inf, Inf),
                                         dist_type = dist_type, filter_dist = filter_dist,
                                         saveto = paste0(saveto_3, "/"))
  solutions$RG_All$Name <- "Red Giants"
   
  cat("Main Sequence processing", "\n")
  saveto_3 <- paste0(saveto_2, "/MS_ALL")
  if (!dir.exists(saveto_3)) 
    dir.create(saveto_3)
  solutions$MS_All <-  GDR3_calc_OM_cond(data, lclass = 5, population = "ALL", 
                                         model = 2, type = 2, 
                                         g_b = c(-Inf, Inf),
                                         use = c(TRUE, TRUE, FALSE),
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
  solutions$RG_North <-  GDR3_calc_OM_cond(data, lclass = 3, population = "ALL",
                                           model = 2, type = 3,
                                           use = c(TRUE, TRUE, FALSE),
                                           dist_type = dist_type, filter_dist = filter_dist,
                                           g_b = c(0, Inf),
                                           saveto = paste0(saveto_3, "/"))
  solutions$RG_North$Name <- "Red Giants North"
  
  cat("Main Sequence processing North", "\n")
  saveto_3 <- paste0(saveto_2, "/MS_North")
  if (!dir.exists(saveto_3)) 
    dir.create(saveto_3)
  solutions$MS_North <-  GDR3_calc_OM_cond(data, lclass = 5, population = "ALL", 
                                           model = 2, type = 3, 
                                           use = c(TRUE, TRUE, FALSE),
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
  solutions$RG_South <-  GDR3_calc_OM_cond(data, lclass = 3, population = "ALL",
                                           model = 2, type = 3,
                                           use = c(TRUE, TRUE, FALSE),
                                           dist_type = dist_type, filter_dist = filter_dist,
                                           g_b = c(-Inf, 0),
                                           saveto = paste0(saveto_3, "/"))
  solutions$RG_South$Name <- "Red Giants South"
  
  cat("Main Sequence processing", "\n")
  saveto_3 <- paste0(saveto_2, "/MS_South")
  if (!dir.exists(saveto_3)) 
    dir.create(saveto_3)
  solutions$MS_South <-  GDR3_calc_OM_cond(data, lclass = 5, population = "ALL", 
                                           model = 2, type = 3, 
                                           use = c(TRUE, TRUE, FALSE),
                                           dist_type = dist_type, filter_dist = filter_dist, 
                                           g_b = c(-Inf, 0),
                                           saveto = paste0(saveto_3, "/"))
  solutions$MS_South$Name <- "Main Sequence South"
  
  
  cat("Calc physical parameters...", "\n")
  solutions <- calc_all_physical_params(solutions)
  cat("Calc weited parameters...", "\n")
  solutions <- calc_all_weighted(solutions)
  
  cat("Export solutions...", "\n")
  export_all_solution(solutions)
  # cat("Export kinematics...", "\n")
  # draw_all_kinematic(solutions, saveto = saveto_, xlim = c(-0.5, 1.6, 0.1), xpar = 9)
  # cat("Export physics...", "\n")
  
  return(solutions)
  
}

#==============================================================================================

GDR3_make_Oort_solutions_bv <- function(data, filter_dist = "RGEO", name = "Hemisphere", ph = "APASS", path = "/Catalogues/Gaia/")
{
  
  
  conditions <- list();
  conditions$Date <- Sys.time();
  conditions$Src <- "TGAS";
  conditions$Filter_Dist <- filter_dist;
  conditions$use <- c(TRUE, TRUE, FALSE);
  conditions$Photometry <- ph
  
  solutions_bv <- list()
  
  # conditions$BV <- matrix(0, nrow = 7, ncol = 2)
  # conditions$BV[,1] <- c(-Inf, -0.30, 0.00, 0.30, 0.58, 0.85, 1.42)
  # conditions$BV[,2] <- c(-0.30, 0.00, 0.30, 0.58, 0.85, 1.42, Inf)
  
  
  #conditions$BV <- matrix(0, nrow = 19, ncol = 2) 
  #conditions$BV[,1] <- c(-Inf, seq(-0.3, 1.4, 0.1))
  #conditions$BV[,2] <- c(seq(-0.3, 1.4, 0.1), Inf)
  
  conditions$Z <- c(0, Inf)
  conditions$MG <- c(-Inf, Inf)
  conditions$e_Px <- Inf
  conditions$distance_ <- c(0, Inf)
  conditions$Dist_Type <- "RGEO"
  
  if (!dir.exists(paste0(path, "Solutions")))
    dir.create(paste0(path, "Solutions"))
  
  saveto_ <- paste0(path, "Solutions/solution_",name,"_",Sys.Date(), "_Oort_BV", conditions$KinModel)
  if (!dir.exists(saveto_)) 
    dir.create(saveto_)
  
  # ----------------------------  Whole Sphere, Main Sequence ------------------------------  
  
  conditions$LClass <- 5
  conditions$KinModel <- 2
  conditions$KinModelType <- 2
  conditions$g_B <- c(-Inf, Inf)
  conditions$BV <- matrix(0, nrow = 19, ncol = 2)
  conditions$BV[,1] <- c(-Inf, -0.30, 0.00, 0.10, 0.20, 0.30, 0.34, 0.37, 0.42, 0.47, 0.52,  0.58, 0.61, 0.65, 0.69, 0.75,  0.85, 1.16, 1.42)
  conditions$BV[,2] <- c(-0.30, 0.00, 0.10, 0.20, 0.30, 0.34, 0.37, 0.42, 0.47, 0.52, 0.58, 0.61, 0.65, 0.69, 0.75, 0.85, 1.16, 1.42, Inf)
  
  saveto_2 <- paste0(saveto_, "/MS_SPHERE_",filter_dist, "-", conditions$Dist_Type)
  if (!dir.exists(saveto_2)) 
    dir.create(saveto_2)
  
  conditions$SaveTo <- paste0(saveto_2, "/")
  
  solution_bv  <- GDR3_calc_OM_seq_2(data, src_ = conditions$Src,
                                     z_lim = conditions$Z, 
                                     e_px = conditions$e_Px, 
                                     bv = conditions$BV, 
                                     Mg = conditions$MG,
                                     px_type = "DIST",
                                     distance = conditions$distance_,
                                     save = conditions$SaveTo,
                                     type = conditions$KinModelType, model = conditions$KinModel,
                                     dist_type = conditions$Dist_Type, use = conditions$use,
                                     g_b = conditions$g_B, 
                                     lclass = conditions$LClass)
  
  solution_bv$Conditions <- conditions
  
  solution_bv <- process_solution(solution_bv)
  solutions_bv$MS_Sph <- solution_bv
  solutions_bv$MS_Sph$Name <- paste("1. Main sequence, Sphere", solution_bv$Conditions$Filter_Dist, solution_bv$Conditions$Dist_Type)
  
  gc()
  
  # ----------------------------  Whole Sphere, Main Sequence ------------------------------  
  # 
  # conditions$LClass <- 5
  # conditions$KinModel <- 2
  # conditions$KinModelType <- 3
  # conditions$g_B <- c(-Inf, Inf)
  # conditions$BV <- matrix(0, nrow = 19, ncol = 2)
  # conditions$BV[,1] <- c(-Inf, -0.30, 0.00, 0.10, 0.20, 0.30, 0.34, 0.37, 0.42, 0.47, 0.52,  0.58, 0.61, 0.65, 0.69, 0.75,  0.85, 1.16, 1.42)
  # conditions$BV[,2] <- c(-0.30, 0.00, 0.10, 0.20, 0.30, 0.34, 0.37, 0.42, 0.47, 0.52, 0.58, 0.61, 0.65, 0.69, 0.75, 0.85, 1.16, 1.42, Inf)
  # 
  # saveto_2 <- paste0(saveto_, "/MS_SPHERE_TEXT_",filter_dist, "-", conditions$Dist_Type)
  # if (!dir.exists(saveto_2)) 
  #   dir.create(saveto_2)
  # 
  # conditions$SaveTo <- paste0(saveto_2, "/")
  # 
  # solution_bv  <- GDR3_calc_OM_seq_2(data, src_ = conditions$Src,
  #                                    z_lim = conditions$Z, 
  #                                    e_px = conditions$e_Px, 
  #                                    bv = conditions$BV, 
  #                                    Mg = conditions$MG,
  #                                    px_type = "DIST",
  #                                    distance = conditions$distance_,
  #                                    save = conditions$SaveTo,
  #                                    type = conditions$KinModelType, model = conditions$KinModel,
  #                                    dist_type = conditions$Dist_Type, use = conditions$use,
  #                                    g_b = conditions$g_B, 
  #                                    lclass = conditions$LClass)
  # 
  # solution_bv$Conditions <- conditions
  # 
  # solution_bv <- process_solution(solution_bv)
  # solutions_bv$MS_Sph_Test <- solution_bv
  # solutions_bv$MS_Sph_Test$Name <- paste("1. Main sequence, Sphere Test", solution_bv$Conditions$Filter_Dist, solution_bv$Conditions$Dist_Type)
  # 
  # gc()
  
  # -------------------------- North Hemisphere, Main Sequence --------------------------------
  
  conditions$LClass <- 5
  conditions$KinModel <- 2
  conditions$KinModelType <- 3
  conditions$g_B <- c(0, Inf)
  conditions$BV <- matrix(0, nrow = 19, ncol = 2)
  conditions$BV[,1] <- c(-Inf, -0.30, 0.00, 0.10, 0.20, 0.30, 0.34, 0.37, 0.42, 0.47, 0.52,  0.58, 0.61, 0.65, 0.69, 0.75,  0.85, 1.16, 1.42)
  conditions$BV[,2] <- c(-0.30, 0.00, 0.10, 0.20, 0.30, 0.34, 0.37, 0.42, 0.47, 0.52, 0.58, 0.61, 0.65, 0.69, 0.75, 0.85, 1.16, 1.42, Inf)
  
  saveto_2 <- paste0(saveto_, "/MS_NORTH_",filter_dist, "-", conditions$Dist_Type)
  if (!dir.exists(saveto_2)) 
    dir.create(saveto_2)
  
  conditions$SaveTo <- paste0(saveto_2, "/")
  
  solution_bv  <- GDR3_calc_OM_seq_2(data, src_ = conditions$Src,
                                     z_lim = conditions$Z, 
                                     e_px = conditions$e_Px, 
                                     bv = conditions$BV, 
                                     Mg = conditions$MG,
                                     px_type = "DIST",
                                     distance = conditions$distance_,
                                     save = conditions$SaveTo,
                                     type = conditions$KinModelType, model = conditions$KinModel,
                                     dist_type = conditions$Dist_Type, use = conditions$use,
                                     g_b = conditions$g_B, 
                                     lclass = conditions$LClass)
  
  solution_bv$Conditions <- conditions
  
  solution_bv <- process_solution(solution_bv)
  solutions_bv$MS_North <- solution_bv
  solutions_bv$MS_North$Name <- paste("2. Main sequence, North Hemisphere", solution_bv$Conditions$Filter_Dist, solution_bv$Conditions$Dist_Type)
  
  rm(solution_bv)
  
  gc()
  
  # -------------------------- South Hemisphere, Main Sequence --------------------------------
  
  conditions$LClass <- 5
  conditions$KinModel <- 2
  conditions$KinModelType <- 3
  conditions$g_B <- c(-Inf, 0)
  conditions$BV <- matrix(0, nrow = 19, ncol = 2)
  conditions$BV[,1] <- c(-Inf, -0.30, 0.00, 0.10, 0.20, 0.30, 0.34, 0.37, 0.42, 0.47, 0.52,  0.58, 0.61, 0.65, 0.69, 0.75,  0.85, 1.16, 1.42)
  conditions$BV[,2] <- c(-0.30, 0.00, 0.10, 0.20, 0.30, 0.34, 0.37, 0.42, 0.47, 0.52, 0.58, 0.61, 0.65, 0.69, 0.75, 0.85, 1.16, 1.42, Inf)
  
  saveto_2 <- paste0(saveto_, "/MS_SOUTH_",filter_dist, "-", conditions$Dist_Type)
  if (!dir.exists(saveto_2)) 
    dir.create(saveto_2)
  
  conditions$SaveTo <- paste0(saveto_2, "/")
  
  solution_bv  <- GDR3_calc_OM_seq_2(data, src_ = conditions$Src,
                                     z_lim = conditions$Z, 
                                     e_px = conditions$e_Px, 
                                     bv = conditions$BV, 
                                     Mg = conditions$MG,
                                     px_type = "DIST",
                                     distance = conditions$distance_,
                                     save = conditions$SaveTo,
                                     type = conditions$KinModelType, model = conditions$KinModel,
                                     dist_type = conditions$Dist_Type, use = conditions$use,
                                     g_b = conditions$g_B, 
                                     lclass = conditions$LClass)
  
  solution_bv$Conditions <- conditions
  
  solution_bv <- process_solution(solution_bv)
  solutions_bv$MS_South <- solution_bv
  solutions_bv$MS_South$Name <- paste("3. Main sequence, South Hemisphere", solution_bv$Conditions$Filter_Dist, solution_bv$Conditions$Dist_Type)
  
  rm(solution_bv)
  
  gc()
  
  # -------------------------- Sphere, Red Giants --------------------------------

  conditions$LClass <- 3
  conditions$KinModel <- 2
  conditions$KinModelType <- 2
  conditions$g_B <- c(-Inf, Inf)
  conditions$BV <- matrix(0, nrow = 5, ncol = 2)
  conditions$BV[,1] <- c(0.84, 1.03, 1.45, 1.57, 1.8)
  conditions$BV[,2] <- c(1.03, 1.45, 1.57, 1.8, 2.5)

  saveto_2 <- paste0(saveto_, "/RG_SPHERE_",filter_dist, "-", conditions$Dist_Type)
  if (!dir.exists(saveto_2))
    dir.create(saveto_2)

  conditions$SaveTo <- paste0(saveto_2, "/")

  solution_bv  <- GDR3_calc_OM_seq_2(data, src_ = conditions$Src,
                                     z_lim = conditions$Z,
                                     e_px = conditions$e_Px,
                                     bv = conditions$BV,
                                     Mg = conditions$MG,
                                     px_type = "DIST",
                                     distance = conditions$distance_,
                                     save = conditions$SaveTo,
                                     type = conditions$KinModelType, model = conditions$KinModel,
                                     dist_type = conditions$Dist_Type, use = conditions$use,
                                     g_b = conditions$g_B,
                                     lclass = conditions$LClass)

  solution_bv$Conditions <- conditions

  solution_bv <- process_solution(solution_bv)
  solutions_bv$RG_Sphere <- solution_bv
  solutions_bv$RG_Sphere$Name <- paste("4. Reg Giants, Sphere", solution_bv$Conditions$Filter_Dist, solution_bv$Conditions$Dist_Type)

  rm(solution_bv)

  gc()

  # -------------------------- North Hemisphere, Red Giants --------------------------------

  conditions$LClass <- 3
  conditions$KinModel <- 2
  conditions$KinModelType <- 3
  conditions$g_B <- c(0, Inf)
  conditions$BV <- matrix(0, nrow = 5, ncol = 2)
  conditions$BV[,1] <- c(0.84, 1.03, 1.45, 1.57, 1.8)
  conditions$BV[,2] <- c(1.03, 1.45, 1.57, 1.8, 2.5)

  saveto_2 <- paste0(saveto_, "/RG_NORTH_",filter_dist, "-", conditions$Dist_Type)
  if (!dir.exists(saveto_2))
    dir.create(saveto_2)

  conditions$SaveTo <- paste0(saveto_2, "/")

  solution_bv  <- GDR3_calc_OM_seq_2(data, src_ = conditions$Src,
                                     z_lim = conditions$Z,
                                     e_px = conditions$e_Px,
                                     bv = conditions$BV,
                                     Mg = conditions$MG,
                                     px_type = "DIST",
                                     distance = conditions$distance_,
                                     save = conditions$SaveTo,
                                     type = conditions$KinModelType, model = conditions$KinModel,
                                     dist_type = conditions$Dist_Type, use = conditions$use,
                                     g_b = conditions$g_B,
                                     lclass = conditions$LClass)

  solution_bv$Conditions <- conditions

  solution_bv <- process_solution(solution_bv)
  solutions_bv$RG_North <- solution_bv
  solutions_bv$RG_North$Name <- paste("5. Red Giants, North Hemisphere", solution_bv$Conditions$Filter_Dist, solution_bv$Conditions$Dist_Type)

  rm(solution_bv)

  gc()

  # -------------------------- South Hemisphere, Red Giants --------------------------------

  conditions$LClass <- 3
  conditions$KinModel <- 2
  conditions$KinModelType <- 3
  conditions$g_B <- c(-Inf, 0)
  conditions$BV <- matrix(0, nrow = 5, ncol = 2)
  conditions$BV[,1] <- c(0.84, 1.03, 1.45, 1.57, 1.8)
  conditions$BV[,2] <- c(1.03, 1.45, 1.57, 1.8, 2.5)

  saveto_2 <- paste0(saveto_, "/RG_SOUTH_",filter_dist, "-", conditions$Dist_Type)
  if (!dir.exists(saveto_2))
    dir.create(saveto_2)

  conditions$SaveTo <- paste0(saveto_2, "/")

  solution_bv  <- GDR3_calc_OM_seq_2(data, src_ = conditions$Src,
                                     z_lim = conditions$Z,
                                     e_px = conditions$e_Px,
                                     bv = conditions$BV,
                                     Mg = conditions$MG,
                                     px_type = "DIST",
                                     distance = conditions$distance_,
                                     save = conditions$SaveTo,
                                     type = conditions$KinModelType, model = conditions$KinModel,
                                     dist_type = conditions$Dist_Type, use = conditions$use,
                                     g_b = conditions$g_B,
                                     lclass = conditions$LClass)

  solution_bv$Conditions <- conditions

  solution_bv <- process_solution(solution_bv)
  solutions_bv$RG_South <- solution_bv
  solutions_bv$RG_South$Name <- paste("6. Red Giants, South Hemisphere", solution_bv$Conditions$Filter_Dist, solution_bv$Conditions$Dist_Type)

  rm(solution_bv)

  gc()


  cat("Export solutions...", "\n")
  export_all_solution(solutions_bv)

  # ----------------------------------------------------------
  # g <- draw_OortParameter(solutions_bv, parameter = 1,
  #                    title = "Oort`s parameter A",
  #                    x_lim = c(-0.5, 1.6, 0.1),
  #                    y_lim = c(6, 24, 2),
  #                    clr = c("blue", "green4", "brown", "black", "red", "orange"),
  #                    x_par = 9,
  #                    x_title = "B-V")
  # ggsave(paste0("solutions/",filter_dist,"_OL_A.png"), plot = g, width = 10, height = 5)
  # ggsave(paste0("solutions/",filter_dist,"_OL_A.eps"), plot = g, width = 10, height = 5)
  #
  # g <- draw_OortParameter(solutions_bv, parameter = 2,
  #                         title = "Oort`s parameter B",
  #                         x_lim = c(-0.5, 1.6, 0.1), y_lim = c(-22, -5, 2),
  #                         clr = c("blue", "green4", "brown", "black", "red", "orange"),
  #                         x_par = 9,
  #                         x_title = "B-V")
  # ggsave(paste0("solutions/",filter_dist,"_OL_B.png"), plot = g, width = 10, height = 5)
  # ggsave(paste0("solutions/",filter_dist,"_OL_B.eps"), plot = g, width = 10, height = 5)
  #
  # g <- draw_OortParameter(solutions_bv, parameter = 3,
  #                         title = "Oort`s parameter C",
  #                         x_lim = c(-0.5, 1.6, 0.1), y_lim = c(-10, 6, 2),
  #                         clr = c("blue", "green4", "brown", "black", "red", "orange"),
  #                         x_par = 9,
  #                         x_title = "B-V")
  # ggsave(paste0("solutions/",filter_dist,"_OL_C.png"), plot = g, width = 10, height = 5)
  # ggsave(paste0("solutions/",filter_dist,"_OL_C.eps"), plot = g, width = 10, height = 5)
  #
  # g <- draw_OortParameter(solutions_bv, parameter = 4,
  #                         title = "Oort`s parameter K",
  #                         x_lim = c(-0.5, 1.6, 0.1), y_lim = c(-10, 6, 2),
  #                         clr = c("blue", "green4", "brown", "black", "red", "orange"),
  #                         x_par = 9,
  #                         x_title = "B-V")
  # ggsave(paste0("solutions/",filter_dist,"_OL_K.png"), plot = g, width = 10, height = 5)
  # ggsave(paste0("solutions/",filter_dist,"_OL_K.eps"), plot = g, width = 10, height = 5)
  
  draw_all_OM_sol_comp(solutions = solutions_bv, 
                       ylims  = matrix(data = c(0, 20, 2, 
                                                0, 50, 2, 
                                                0, 15, 1,  
                                                -20, -5, 2, 
                                                5, 20, 1, 
                                                -7, 2, 1, 
                                                -20, 20, 2, 
                                                -10, 10, 2, 
                                                -60, 60, 5, 
                                                -10, 10, 2, 
                                                -10, 5, 2), nrow = 3),
                       xlims = c(-0.5, 1.6, 0.1),
                       xpar = 9, 
                       xtitle = "B-V", 
                       saveto = paste0(saveto_, "/"))
  cat("Export kinematics...", "\n")
  draw_all_kinematic(solutions_bv, saveto = saveto_, xlim = c(-0.2, 1.6, 0.2), xpar = 9, x_title = "B-V", src = "Gaia EDR3")
  cat("Export physics...", "\n")
  draw_physics(solutions_bv, saveto = saveto_, x_lim = c(-0.2, 1.6, 0.2), x_par = 9,  x_title = "B-V", src = "Gaia EDR3")
  
  return(solutions_bv)
  
}

#==============================================================================================

GDR3_make_OM_solutions_bv <- function(data, filter_dist = "GAIA_PX", name = "BV", ph = "APASS", path = "/Catalogues/Gaia/")
{

  
  conditions <- list();
  conditions$Src <- "TGAS";
  conditions$Filter_Dist <- filter_dist;
  conditions$use <- c(TRUE, TRUE, FALSE);
  conditions$KinModel <- 1
  conditions$KinModelType <- 1
  conditions$g_B <- c(-Inf, Inf)
  conditions$Photometry <- ph
  
  solutions_bv <- list()
  
  # conditions$BV <- matrix(0, nrow = 7, ncol = 2)
  # conditions$BV[,1] <- c(-Inf, -0.30, 0.00, 0.30, 0.58, 0.85, 1.42)
  # conditions$BV[,2] <- c(-0.30, 0.00, 0.30, 0.58, 0.85, 1.42, Inf)
  
  conditions$BV <- matrix(0, nrow = 19, ncol = 2)
  conditions$BV[,1] <- c(-Inf, -0.30, 0.00, 0.10, 0.20, 0.30, 0.34, 0.37, 0.42, 0.47, 0.52,  0.58, 0.61, 0.65, 0.69, 0.75,  0.85, 1.16, 1.42)
  conditions$BV[,2] <- c(-0.30, 0.00, 0.10, 0.20, 0.30, 0.34, 0.37, 0.42, 0.47, 0.52, 0.58, 0.61, 0.65, 0.69, 0.75, 0.85, 1.16, 1.42, Inf)
  
  #conditions$BV <- matrix(0, nrow = 19, ncol = 2) 
  #conditions$BV[,1] <- c(-Inf, seq(-0.3, 1.4, 0.1))
  #conditions$BV[,2] <- c(seq(-0.3, 1.4, 0.1), Inf)
  
  conditions$Z <- c(0, Inf)
  conditions$MG <- c(-Inf, Inf)
  conditions$e_Px <- Inf
  conditions$distance_ <- c(0, Inf)
  conditions$LClass <- 5
  
  
  if (!dir.exists(paste0(path, "Solutions")))
    dir.create(paste0(path, "Solutions"))
  
  saveto_ <- paste0(path, "Solutions/solution_",name,"_",Sys.Date(), "_OM", conditions$KinModel)
  if (!dir.exists(saveto_)) 
    dir.create(saveto_)
  
  # ----------------------------------------------------------  
  conditions$Dist_Type <- "GAIA_PX"
  
  saveto_2 <- paste0(saveto_, "/MS_ALL_",filter_dist, "-", conditions$Dist_Type)
  if (!dir.exists(saveto_2)) 
    dir.create(saveto_2)
  
  conditions$SaveTo <- paste0(saveto_2, "/")
  
  solution_bv  <- GDR3_calc_OM_seq_2(data, src_ = conditions$Src,
                                     z_lim = conditions$Z, 
                                     e_px = conditions$e_Px, 
                                     bv = conditions$BV, 
                                     Mg = conditions$MG,
                                     px_type = "DIST",
                                     distance = conditions$distance_,
                                     save = conditions$SaveTo,
                                     type = conditions$KinModelType, model = conditions$KinModel,
                                     dist_type = conditions$Dist_Type, use = conditions$use,
                                     g_b = conditions$g_B, 
                                     lclass = conditions$LClass)
  
  solution_bv$Conditions <- conditions
  
  solution_bv <- process_solution(solution_bv)
  solutions_bv$MS_ALL_Rpi <- solution_bv
  solutions_bv$MS_ALL_Rpi$Name <- paste("1. Main sequence", solution_bv$Conditions$Filter_Dist, solution_bv$Conditions$Dist_Type)
  
  rm(solution_bv)
  gc()
  
  # ----------------------------------------------------------
  conditions$Dist_Type <- "RGEO"
  
  
  saveto_2 <- paste0(saveto_, "/MS_ALL_",filter_dist, "-", conditions$Dist_Type)
  if (!dir.exists(saveto_2))
    dir.create(saveto_2)
  
  conditions$SaveTo <- paste0(saveto_2, "/")
  
  solution_bv  <- GDR3_calc_OM_seq_2(data, src_ = conditions$Src,
                                     z_lim = conditions$Z, e_px = conditions$e_Px, bv = conditions$BV, Mg = conditions$MG,
                                     px_type = "DIST", 
                                     distance = conditions$distance_,
                                     save = conditions$SaveTo,
                                     type = conditions$KinModelType, model = conditions$KinModel,
                                     dist_type = conditions$Dist_Type, use = conditions$use,
                                     g_b = conditions$g_B, 
                                     lclass = conditions$LClass)
  
  solution_bv$Conditions <- conditions
  
  solution_bv <- process_solution(solution_bv)
  solutions_bv$MS_ALL_Rgeo <- solution_bv
  solutions_bv$MS_ALL_Rgeo$Name <- paste("2. Main sequence", solution_bv$Conditions$Filter_Dist, solution_bv$Conditions$Dist_Type)
  
  rm(solution_bv)
  gc()
  
  # ----------------------------------------------------------
  
  conditions$Dist_Type <- "RPGEO"


  saveto_2 <- paste0(saveto_, "/MS_ALL_",filter_dist, "-", conditions$Dist_Type)
  if (!dir.exists(saveto_2))
    dir.create(saveto_2)

  conditions$SaveTo <- paste0(saveto_2, "/")

  solution_bv  <- GDR3_calc_OM_seq_2(data, src_ = conditions$Src,
                                     z_lim = conditions$Z, e_px = conditions$e_Px, bv = conditions$BV, Mg = conditions$MG,
                                     px_type = "DIST",
                                     distance = conditions$distance_,
                                     save = conditions$SaveTo,
                                     type = conditions$KinModelType, model = conditions$KinModel,
                                     dist_type = conditions$Dist_Type, use = conditions$use,
                                     g_b = conditions$g_B,
                                     lclass = conditions$LClass)
  solution_bv$Conditions <- conditions

  solution_bv <- process_solution(solution_bv)
  solutions_bv$MS_ALL_Rpgeo <- solution_bv
  solutions_bv$MS_ALL_Rpgeo$Name <- paste("3. Main sequence", solution_bv$Conditions$Filter_Dist, solution_bv$Conditions$Dist_Type)
  
  rm(solution_bv)
  gc()
  
  # ----------------------------------------------------------
  
  
  saveRDS(solutions_bv, file = paste0(saveto_, name, ".RDS"))
  
  g <- draw_OortParameter(solutions_bv, parameter = 1,
                     title = "Oort`s parameter A",
                     x_lim = c(-0.5, 1.6, 0.1),
                     y_lim = c(6, 24, 2),
                     clr = c("blue", "green4", "brown", "black", "red", "orange"),
                     x_par = 9,
                     x_title = "B-V")
  ggsave(paste0(saveto_,"/","_OL_A.png"), plot = g, width = 10, height = 5)
  ggsave(paste0(saveto_,"/","_OL_A.eps"), plot = g, width = 10, height = 5)

  g <- draw_OortParameter(solutions_bv, parameter = 2,
                          title = "Oort`s parameter B",
                          x_lim = c(-0.5, 1.6, 0.1), y_lim = c(-22, -5, 2),
                          clr = c("blue", "green4", "brown", "black", "red", "orange"),
                          x_par = 9,
                          x_title = "B-V")
  ggsave(paste0(saveto_,"/","_OL_B.png"), plot = g, width = 10, height = 5)
  ggsave(paste0(saveto_,"/","_OL_B.eps"), plot = g, width = 10, height = 5)

  g <- draw_OortParameter(solutions_bv, parameter = 3,
                          title = "Oort`s parameter C",
                          x_lim = c(-0.5, 1.6, 0.1), y_lim = c(-10, 6, 2),
                          clr = c("blue", "green4", "brown", "black", "red", "orange"),
                          x_par = 9,
                          x_title = "B-V")
  ggsave(paste0(saveto_, "/","_OL_C.png"), plot = g, width = 10, height = 5)
  ggsave(paste0(saveto_, "/","_OL_C.eps"), plot = g, width = 10, height = 5)

  g <- draw_OortParameter(solutions_bv, parameter = 4,
                          title = "Oort`s parameter K",
                          x_lim = c(-0.5, 1.6, 0.1), y_lim = c(-10, 6, 2),
                          clr = c("blue", "green4", "brown", "black", "red", "orange"),
                          x_par = 9,
                          x_title = "B-V")
  ggsave(paste0(saveto_,"/","_OL_K.png"), plot = g, width = 10, height = 5)
  ggsave(paste0(saveto_,"/","_OL_K.eps"), plot = g, width = 10, height = 5)
  
  draw_all_OM_sol_comp(solutions = solutions_bv, 
                            ylims  = matrix(data = c(0, 20, 1, 
                                                     0, 35, 1, 
                                                     0, 15, 1, 
                                                     -2, 10, 1, 
                                                     -5, 2, 1, 
                                                     -25, -5, 1, 
                                                     -5, 5, 1, 
                                                     -8, 5, 1, 
                                                     5, 30, 1, 
                                                     -10, 3, 1, 
                                                     -10, 5, 1), nrow =3),
                            xlims = c(-0.5, 1.6, 0.1),
                            xpar = 9, 
                            xtitle = "B-V", 
                            saveto = paste0(saveto_, "/"))
  
  cat("Export solutions...", "\n")
  export_all_solution(solutions_bv)
  cat("Export kinematics...", "\n")
  draw_all_kinematic(solutions_bv, saveto = saveto_, xlim = c(-0.5, 1.6, 0.1), xpar = 9)
  cat("Export physics...", "\n")
  draw_physics(solutions_bv, saveto = saveto_, x_lim = c(-0.5, 1.6, 0.1), x_par = 9)
  
  return(solutions_bv)
  
}

#==============================================================================================

GDR3_make_bottlinger_solutions_bv <- function(data, filter_dist = "GAIA_PX", name = "BV", ph = "APASS", path = "/Catalogues/Gaia/")
{
  
  solutions_bv <- list()

  conditions <- list();
  conditions$Src <- "GAIA";
  conditions$Filter_Dist <- filter_dist;
  conditions$use <- c(TRUE, TRUE, FALSE);
  conditions$KinModel <- 4
  conditions$KinModelType <- 1
  conditions$g_B <- c(-Inf, Inf)
  
  # conditions$BV <- matrix(0, nrow = 7, ncol = 2)
  # conditions$BV[,1] <- c(-Inf, -0.30, 0.00, 0.30, 0.58, 0.85, 1.42)
  # conditions$BV[,2] <- c(-0.30, 0.00, 0.30, 0.58, 0.85, 1.42, Inf)
  
  conditions$BV <- matrix(0, nrow = 19, ncol = 2) 
  conditions$BV[,1] <- c(-Inf, -0.30, 0.00, 0.10, 0.20, 0.30, 0.34, 0.37, 0.42, 0.47, 0.52,  0.58, 0.61, 0.65, 0.69, 0.75,  0.85, 1.16, 1.42)
  conditions$BV[,2] <- c(-0.30, 0.00, 0.10, 0.20, 0.30, 0.34, 0.37, 0.42, 0.47, 0.52, 0.58, 0.61, 0.65, 0.69, 0.75, 0.85, 1.16, 1.42, Inf)
  
  # conditions$BV <- matrix(0, nrow = 19, ncol = 2) 
  # conditions$BV[,1] <- c(-Inf, seq(-0.3, 1.4, 0.1))
  # conditions$BV[,2] <- c(seq(-0.3, 1.4, 0.1), Inf)
  
  conditions$Z <- c(0, Inf)
  conditions$MG <- c(-Inf, Inf)
  conditions$e_Px <- Inf
  conditions$distance_ <- c(0, Inf)

  
  # ----------------------------------------------------------  
  
  conditions$Dist_Type <- "RGEO"
  conditions$LClass <- 5
  
  if (!dir.exists(paste0(path, "Solutions")))
    dir.create(paste0(path, "Solutions"))
  
  saveto_ <- paste0(path, "Solutions/solution_",name,"_",Sys.Date(), "_", filter_dist, "-", conditions$Dist_Type, "_Bottlinger", conditions$KinModel)
  if (!dir.exists(saveto_)) 
    dir.create(saveto_)
  
  saveto_2 <- paste0(saveto_, "/MS_ALL")
  if (!dir.exists(saveto_2)) 
    dir.create(saveto_2)
  
  conditions$SaveTo <- paste0(saveto_2, "/")
  
  solution_bv  <- GDR3_calc_OM_seq_2(data, src_ = conditions$Src,
                                     z_lim = conditions$Z, 
                                     e_px = conditions$e_Px, 
                                     bv = conditions$BV, 
                                     Mg = conditions$MG,
                                     lclass = conditions$LClass,
                                     px_type = "DIST",
                                     distance = distance_,
                                     save = conditions$SaveTo,
                                     type = conditions$KinModelType, model = conditions$KinModel,
                                     dist_type = conditions$Dist_Type, use = conditions$use,
                                     g_b = conditions$g_B)
  solution_bv$Conditions <- conditions
  
  solution_bv <- process_solution(solution_bv)
  solutions_bv$MS_ALL_Rpi <- solution_bv
  solutions_bv$MS_ALL_Rpi$Name <- paste("1. Main sequence", solution_bv$Conditions$Filter_Dist, solution_bv$Conditions$Dist_Type)
  
  rm(solution_bv)
  gc()
  
  # ----------------------------------------------------------  
  # conditions$Dist_Type <- "RGEO"
  # conditions$LClass <- 3
  # 
  # saveto_ <- paste0(path, "Solutions/solution_",name,"_",Sys.Date(), "_", filter_dist, "-", conditions$Dist_Type, "_Bottlinger", conditions$KinModel)
  # if (!dir.exists(saveto_)) 
  #   dir.create(saveto_)
  # 
  # saveto_2 <- paste0(saveto_, "/RG_ALL")
  # if (!dir.exists(saveto_2)) 
  #   dir.create(saveto_2)
  # 
  # conditions$SaveTo <- paste0(saveto_2, "/")
  # 
  # solution_bv  <- GDR3_calc_OM_seq_2(data, src_ = conditions$Src,
  #                                    z_lim = conditions$Z, e_px = conditions$e_Px, bv = conditions$BV, Mg = conditions$MG,
  #                                    px_type = "DIST", distance = distance_,
  #                                    lclass = conditions$LClass,
  #                                    save = conditions$SaveTo,
  #                                    type = conditions$KinModelType, model = conditions$KinModel,
  #                                    dist_type = conditions$Dist_Type, use = conditions$use,
  #                                    g_b = conditions$g_B)
  # solution_bv$Conditions <- conditions
  # 
  # solution_bv <- process_solution(solution_bv)
  # solutions_bv$MS_ALL_Rgeo <- solution_bv
  # solutions_bv$MS_ALL_Rgeo$Name <- paste("2. Red Giants", solution_bv$Conditions$Filter_Dist, solution_bv$Conditions$Dist_Type)
  # #solutions_bv$MS_ALL_Rgeo$Name <- "Main sequence Rpi-rMoMW"
  # 
  # rm(solution_bv)
  # gc()
  # ----------------------------------------------------------    
  
  # conditions$Dist_Type <- "RPGEO"
  # 
  # saveto_ <- paste0(path, "Solutions/solution_",name,"_",Sys.Date(), "_", filter_dist, "-", conditions$Dist_Type, "_Bottlinger", conditions$KinModel)
  # if (!dir.exists(saveto_))
  #   dir.create(saveto_)
  # 
  # saveto_2 <- paste0(saveto_, "/MS_ALL")
  # if (!dir.exists(saveto_2))
  #   dir.create(saveto_2)
  # 
  # conditions$SaveTo <- paste0(saveto_2, "/")
  # 
  # solution_bv  <- GDR3_calc_OM_seq_2(data, src_ = conditions$Src,
  #                                    z_lim = conditions$Z, e_px = conditions$e_Px, bv = conditions$BV, Mg = conditions$MG,
  #                                    px_type = "DIST", distance = distance_,
  #                                    save = conditions$SaveTo,
  #                                    type = conditions$KinModelType, model = conditions$KinModel,
  #                                    dist_type = conditions$Dist_Type, use = conditions$use,
  #                                    g_b = conditions$g_B)
  # solution_bv$Conditions <- conditions
  # 
  # solution_bv <- process_solution(solution_bv)
  # solutions_bv$MS_ALL_RPgeo <- solution_bv
  # solutions_bv$MS_ALL_RPgeo$Name <- paste("3. Main sequence", solution_bv$Conditions$Filter_Dist, solution_bv$Conditions$Dist_Type)
  # 
  # gc()
  # ----------------------------------------------------------    

  saveRDS(solutions_bv, file = paste0(saveto_, name, ".RDS"))
    
  g <- draw_OortParameter(solutions_bv, parameter = 1,
                          title = "Oort`s parameter A", 
                          x_lim = c(-0.5, 1.6, 0.1), y_lim = c(6, 24, 2), 
                          clr = c("blue", "green4", "brown", "black", "red", "orange"),
                          x_par = 9, 
                          x_title = "B-V")
  ggsave(paste0(saveto_, "/Bottlinger_OL_A.png"), plot = g, width = 10, height = 5)
  ggsave(paste0(saveto_, "/Bottlinger_OL_A.eps"), plot = g, width = 10, height = 5)
  
  g <- draw_OortParameter(solutions_bv, parameter = 2,
                          title = "Oort`s parameter B", 
                          x_lim = c(-0.5, 1.6, 0.1), y_lim = c(-18, -8, 2), 
                          clr = c("blue", "green4", "brown", "black", "red", "orange"),
                          x_par = 9, 
                          x_title = "B-V")
  ggsave(paste0(saveto_, "/Bottlinger_OL_B.png"), plot = g, width = 10, height = 5)
  ggsave(paste0(saveto_, "/Bottlinger_OL_B.eps"), plot = g, width = 10, height = 5)
  
  g <- draw_OortParameter(solutions_bv, parameter = 4,
                          title = "Oort`s parameter K", 
                          x_lim = c(-0.5, 1.6, 0.1), y_lim = c(-10, 6, 2), 
                          clr = c("blue", "green4", "brown", "black", "red", "orange"),
                          x_par = 9, 
                          x_title = "B-V")
  ggsave(paste0(saveto_, "/Bottlinger_OL_K.png"), plot = g, width = 10, height = 5)
  ggsave(paste0(saveto_, "/Bottlinger_OL_K.eps"), plot = g, width = 10, height = 5)
  
  tgas_draw_all_OM_sol_comp(solutions = solutions_bv, 
                            ylims  = matrix(data = c(-5, 25, 1,
                                                     0, 25, 1, 
                                                     0, 15, 1, 
                                                     15, 35, 1, 
                                                     -10, 0, 1, 
                                                     -7, 10, 1, 
                                                     -10, 5, 1, 
                                                     5, 25, 1, 
                                                     -20, -5, 1), nrow = 3),
                            xlims = c(-0.5, 1.6, 0.1),
                            xpar = 9, 
                            xtitle = "B-V", 
                            saveto = paste0(saveto_,"/Bottlinger_", filter_dist,"_"))
  
  return(solutions_bv)
  
}



