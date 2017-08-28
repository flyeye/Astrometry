
tgas_make_OM_solutions_bv <- function(filter_dist = "TGAS_PX", src = "TGAS", name = "BV")
{

  solutions_bv <- list()
  
  if (!dir.exists("solutions")) 
    dir.create("solutions")
  
  tgas_ <- tgas_calc_LClass(tgas, dist_ = filter_dist)
  tgas_ <- tgas_[tgas_$LClass_apass == 5,]
  
  conditions <- list();
  conditions$Src <- "TGAS";
  conditions$Filter_Dist <- filter_dist;
  conditions$use <- c(TRUE, TRUE, FALSE);
  conditions$KinModel <- 1
  conditions$KinModelType <- 1
  conditions$g_B <- c(-Inf, Inf)
  conditions$BV <- matrix(0, nrow = 7, ncol = 2)
  conditions$BV[,1] <- c(-Inf, -0.30, 0.00, 0.30, 0.58, 0.85, 1.42)
  conditions$BV[,2] <- c(-0.30, 0.00, 0.30, 0.58, 0.85, 1.42, Inf)
  conditions$Z <- c(0, Inf)
  conditions$MG <- c(-Inf, Inf)
  conditions$e_Px <- Inf
  conditions$distance_ <- c(0, Inf)
  conditions$LClass <- 5
  
# ----------------------------------------------------------  
  conditions$Dist_Type <- "TGAS_PX"
  
  saveto_ <- paste0("solutions/solution_", name, "_", filter_dist, "-", conditions$Dist_Type)
  if (!dir.exists(saveto_)) 
    dir.create(saveto_)
  
  saveto_2 <- paste0(saveto_, "/MS_ALL")
  if (!dir.exists(saveto_2)) 
    dir.create(saveto_2)
  
  conditions$SaveTo <- paste0(saveto_2, "/")
  
  solution_bv  <- tgas_calc_OM_seq_2(tgas_, src_ = conditions$Src,
                                     z_lim = conditions$Z, 
                                     e_px = conditions$e_Px, 
                                     bv = conditions$BV, 
                                     Mg = conditions$MG,
                                     px_type = "DIST",
                                     distance = distance_,
                                     save = conditions$SaveTo,
                                     type = conditions$KinModelType, model = conditions$KinModel,
                                     dist_type = conditions$Dist_Type, use = conditions$use,
                                     g_b = conditions$g_B)
  solution_bv$Conditions <- conditions
  tgas_write_conditions(conditions)
  
  solution_bv <- tgas_process_solution(solution_bv)
  solutions_bv$MS_ALL_Rpi <- solution_bv
  solutions_bv$MS_ALL_Rpi$Name <- paste("Main sequence", solution_bv$Conditions$Filter_Dist, solution_bv$Conditions$Dist_Type)
  
  gc()
  
  # ----------------------------------------------------------  
  conditions$Dist_Type <- "rMoMW"
  
  saveto_ <- paste0("solutions/solution_", name, "_", filter_dist, "-", conditions$Dist_Type)
  if (!dir.exists(saveto_)) 
    dir.create(saveto_)
  
  saveto_2 <- paste0(saveto_, "/MS_ALL")
  if (!dir.exists(saveto_2)) 
    dir.create(saveto_2)
  
  conditions$SaveTo <- paste0(saveto_2, "/")
  
  solution_bv  <- tgas_calc_OM_seq_2(tgas_, src_ = conditions$Src,
                                     z_lim = conditions$Z, e_px = conditions$e_Px, bv = conditions$BV, Mg = conditions$MG,
                                     px_type = "DIST", distance = distance_,
                                     save = conditions$SaveTo,
                                     type = conditions$KinModelType, model = conditions$KinModel,
                                     dist_type = conditions$Dist_Type, use = conditions$use,
                                     g_b = conditions$g_B)
  solution_bv$Conditions <- conditions
  tgas_write_conditions(conditions)
  
  solution_bv <- tgas_process_solution(solution_bv)
  solutions_bv$MS_ALL_rMoMW <- solution_bv
  solutions_bv$MS_ALL_rMoMW$Name <- paste("Main sequence", solution_bv$Conditions$Filter_Dist, solution_bv$Conditions$Dist_Type)
  #solutions_bv$MS_ALL_rMoMW$Name <- "Main sequence Rpi-rMoMW"

  gc()
  # ----------------------------------------------------------    
  
  conditions$Dist_Type <- "rMoExp1"
  
  saveto_ <- paste0("solutions/solution_", name, "_", filter_dist, "-", conditions$Dist_Type)
  if (!dir.exists(saveto_)) 
    dir.create(saveto_)
  
  saveto_2 <- paste0(saveto_, "/MS_ALL")
  if (!dir.exists(saveto_2)) 
    dir.create(saveto_2)
  
  conditions$SaveTo <- paste0(saveto_2, "/")
  
  solution_bv  <- tgas_calc_OM_seq_2(tgas_, src_ = conditions$Src,
                                     z_lim = conditions$Z, e_px = conditions$e_Px, bv = conditions$BV, Mg = conditions$MG,
                                     px_type = "DIST", distance = distance_,
                                     save = conditions$SaveTo,
                                     type = conditions$KinModelType, model = conditions$KinModel,
                                     dist_type = conditions$Dist_Type, use = conditions$use,
                                     g_b = conditions$g_B)
  solution_bv$Conditions <- conditions
  tgas_write_conditions(conditions)
  
  solution_bv <- tgas_process_solution(solution_bv)
  solutions_bv$MS_ALL_rMoExp1 <- solution_bv
  #solutions_bv$MS_ALL_rMoExp1$Name <- "Main sequence Rpi-rMoExp1"
  solutions_bv$MS_ALL_rMoExp1$Name <- paste("Main sequence", solution_bv$Conditions$Filter_Dist, solution_bv$Conditions$Dist_Type)
  
  gc()
  # ----------------------------------------------------------    
  
  conditions$Dist_Type <- "rMoExp2"
  
  saveto_ <- paste0("solutions/solution_", name, "_", filter_dist, "-", conditions$Dist_Type)
  if (!dir.exists(saveto_)) 
    dir.create(saveto_)
  
  saveto_2 <- paste0(saveto_, "/MS_ALL")
  if (!dir.exists(saveto_2)) 
    dir.create(saveto_2)
  
  conditions$SaveTo <- paste0(saveto_2, "/")
  
  solution_bv  <- tgas_calc_OM_seq_2(tgas_, src_ = conditions$Src,
                                     z_lim = conditions$Z, e_px = conditions$e_Px, bv = conditions$BV, Mg = conditions$MG,
                                     px_type = "DIST", distance = distance_,
                                     save = conditions$SaveTo,
                                     type = conditions$KinModelType, model = conditions$KinModel,
                                     dist_type = conditions$Dist_Type, use = conditions$use,
                                     g_b = conditions$g_B)
  solution_bv$Conditions <- conditions
  tgas_write_conditions(conditions)
  
  solution_bv <- tgas_process_solution(solution_bv)
  solutions_bv$MS_ALL_rMoExp2 <- solution_bv
  #solutions_bv$MS_ALL_rMoExp2$Name <- "Main sequence Rpi-rMoExp2"
  solutions_bv$MS_ALL_rMoExp2$Name <- paste("Main sequence", solution_bv$Conditions$Filter_Dist, solution_bv$Conditions$Dist_Type)
  
  gc()
  # ----------------------------------------------------------      
  g <- draw_OortParameter(solutions_bv, parameter = 1,
                     title = "Oort`s parameter A", 
                     x_lim = c(-0.5, 1.1, 0.1), y_lim = c(6, 24, 2), 
                     clr = c("blue", "green4", "brown", "black", "red", "orange"),
                     x_par = 9, 
                     x_title = "B-V")
  ggsave(paste0("solutions/",filter_dist,"_OL_A.png"), plot = g, width = 10, height = 5)
  ggsave(paste0("solutions/",filter_dist,"_OL_A.eps"), plot = g, width = 10, height = 5)
  
  g <- draw_OortParameter(solutions_bv, parameter = 2,
                          title = "Oort`s parameter B", 
                          x_lim = c(-0.5, 1.1, 0.1), y_lim = c(-18, -8, 2), 
                          clr = c("blue", "green4", "brown", "black", "red", "orange"),
                          x_par = 9, 
                          x_title = "B-V")
  ggsave(paste0("solutions/",filter_dist,"_OL_B.png"), plot = g, width = 10, height = 5)
  ggsave(paste0("solutions/",filter_dist,"_OL_B.eps"), plot = g, width = 10, height = 5)
  
  g <- draw_OortParameter(solutions_bv, parameter = 3,
                          title = "Oort`s parameter C", 
                          x_lim = c(-0.5, 1.1, 0.1), y_lim = c(-10, 6, 2), 
                          clr = c("blue", "green4", "brown", "black", "red", "orange"),
                          x_par = 9, 
                          x_title = "B-V")
  ggsave(paste0("solutions/",filter_dist,"_OL_C.png"), plot = g, width = 10, height = 5)
  ggsave(paste0("solutions/",filter_dist,"_OL_C.eps"), plot = g, width = 10, height = 5)
  
  g <- draw_OortParameter(solutions_bv, parameter = 4,
                          title = "Oort`s parameter K", 
                          x_lim = c(-0.5, 1.1, 0.1), y_lim = c(-10, 6, 2), 
                          clr = c("blue", "green4", "brown", "black", "red", "orange"),
                          x_par = 9, 
                          x_title = "B-V")
  ggsave(paste0("solutions/",filter_dist,"_OL_K.png"), plot = g, width = 10, height = 5)
  ggsave(paste0("solutions/",filter_dist,"_OL_K.eps"), plot = g, width = 10, height = 5)
  
  tgas_draw_all_OM_sol_comp(solutions = solutions_bv, 
                            ylims  = matrix(data = c(5, 20, 0, 25, 0, 15, -2, 10, -5, 2, -18, -8, -5, 5, -8, 3 , 5, 25, -7, 3, -10, 5), nrow = 2),
                            xlims = c(-0.5, 1.1, 0.1),
                            xpar = 9, 
                            xtitle = "B-V", 
                            saveto = paste0("solutions/", filter_dist,"_"))
  
  
  return(solutions_bv)
  
}















