
hip_make_OM_solutions_bv <- function(data, filter_dist = "HIP", src = "TGAS", name = "BV_APASS_C", ph = "HIP")
{

  #index <- (((data$e_Px/data$Px)<0.1) & ((data$parallax_error/data$gPx)<0.1) & (!is.na(data$apasm_v)) & (!is.na(data$apasm_b)) & ((!is.na(data$e_hBV)) & (data$e_hBV !=0 )) )
  index <- ( (!is.na(data$apasm_v)) & (!is.na(data$apasm_b)) & ((!is.na(data$e_hBV)) & (data$e_hBV !=0 )) )
  tgas_ <- data[index,]
  
  #tgas_ <- data
  
  #tgas_ <- tgas_calc_LClass(tgas_, dist_ = filter_dist, ph = ph)
  tgas_ <- tgas_[tgas_$LClass_apass==5,]
  
  cat("input sample size:", nrow(tgas_), "\n")
  
  conditions <- list();
  conditions$Src <- "TGAS";
  conditions$Filter_Dist <- filter_dist;
  conditions$use <- c(TRUE, TRUE, FALSE);
  conditions$KinModel <- 1
  conditions$KinModelType <- 1
  conditions$g_B <- c(-Inf, Inf)
  conditions$Photometry <- ph
  
  solutions_bv <- list()
  
  if (!dir.exists("solutions")) 
    dir.create("solutions")
  
  #conditions$BV <- matrix(0, nrow = 7, ncol = 2)
  #conditions$BV[,1] <- c(-Inf, -0.30, 0.00, 0.30, 0.58, 0.85, 1.42)
  #conditions$BV[,2] <- c(-0.30, 0.00, 0.30, 0.58, 0.85, 1.42, Inf)
  
  #conditions$BV <- matrix(0, nrow = 19, ncol = 2) 
  #conditions$BV[,1] <- c(-Inf, -0.30, 0.00, 0.10, 0.20, 0.30, 0.34, 0.37, 0.42, 0.47, 0.52,  0.58, 0.61, 0.65, 0.69, 0.75,  0.85, 1.16, 1.42)
  #conditions$BV[,2] <- c(-0.30, 0.00, 0.10, 0.20, 0.30, 0.34, 0.37, 0.42, 0.47, 0.52, 0.58, 0.61, 0.65, 0.69, 0.75, 0.85, 1.16, 1.42, Inf)
  
  conditions$BV <- matrix(0, nrow = 19, ncol = 2) 
  conditions$BV[,1] <- c(-Inf, seq(-0.3, 1.4, 0.1))
  conditions$BV[,2] <- c(seq(-0.3, 1.4, 0.1), Inf)
  
  conditions$Z <- c(0, Inf)
  conditions$MG <- c(-Inf, Inf)
  conditions$e_Px <- Inf
  conditions$distance_ <- c(0, Inf)
  conditions$LClass <- 5
  # conditions$is_draw_HR <- false
  # conditions$is_draw_XYZ <- false
  
# ----------------------------------------------------------  
  cat("stage 1 \n\n")
  
  conditions$Dist_Type <- "HIP"
  #tgas_ <- tgas_calc_LClass(tgas_, dist_ = filter_dist, ph = "HIP")
  tgas_$B_V <- tgas_$hBV
  
  saveto_ <- paste0("solutions/solution_", name, "_", filter_dist,"_PH_", ph, "_OM", conditions$KinModel)
  if (!dir.exists(saveto_)) 
    dir.create(saveto_)
  
  saveto_2 <- paste0(saveto_, "/HIP")
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
                                     #save = conditions$SaveTo,
                                     type = conditions$KinModelType, model = conditions$KinModel,
                                     dist_type = conditions$Dist_Type, use = conditions$use,
                                     g_b = conditions$g_B)
  solution_bv$Conditions <- conditions
  tgas_write_conditions(conditions)
  
  solution_bv <- tgas_process_solution(solution_bv)
  solutions_bv$MS_ALL_Rpi <- solution_bv
  solutions_bv$MS_ALL_Rpi$Name <- paste("Hip")
  
  gc()
  
  # ----------------------------------------------------------
  cat("stage 2 \n\n")
    
  conditions$Dist_Type <- "HIP"
  #tgas_ <- tgas_calc_LClass(tgas_, dist_ = filter_dist, ph = "APASS")
  tgas_$B_V <- (tgas_$apasm_b - tgas_$apasm_v)

  saveto_ <- paste0("solutions/solution_", name, "_", filter_dist,"_PH_", ph, "_OM", conditions$KinModel)
  if (!dir.exists(saveto_))
    dir.create(saveto_)

  saveto_2 <- paste0(saveto_, "/APASS")
  if (!dir.exists(saveto_2))
    dir.create(saveto_2)

  conditions$SaveTo <- paste0(saveto_2, "/")

  solution_bv  <- tgas_calc_OM_seq_2(tgas_, src_ = conditions$Src,
                                     z_lim = conditions$Z, e_px = conditions$e_Px, bv = conditions$BV, Mg = conditions$MG,
                                     px_type = "DIST", distance = distance_,
                                     #save = conditions$SaveTo,
                                     type = conditions$KinModelType, model = conditions$KinModel,
                                     dist_type = conditions$Dist_Type, use = conditions$use,
                                     g_b = conditions$g_B)
  solution_bv$Conditions <- conditions
  tgas_write_conditions(conditions)

  solution_bv <- tgas_process_solution(solution_bv)
  solutions_bv$MS_ALL_rMoMW <- solution_bv
  solutions_bv$MS_ALL_rMoMW$Name <- paste("APASS")

  gc()
  # ----------------------------------------------------------
  # cat("stage 3 \n\n")
  # 
  # conditions$Dist_Type <- "HIP"
  # #tgas_ <- tgas_calc_LClass(tgas_, dist_ = filter_dist, ph = "APASS")
  # tgas_$B_V <- (tgas_$apasm_b - tgas_$apasm_v)*0.7673513 + 0.2591304
  # 
  # saveto_ <- paste0("solutions/solution_", name, "_", filter_dist,"_PH_", ph, "_OM", conditions$KinModel)
  # if (!dir.exists(saveto_))
  #   dir.create(saveto_)
  # 
  # saveto_2 <- paste0(saveto_, "/APASS-c")
  # if (!dir.exists(saveto_2))
  #   dir.create(saveto_2)
  # 
  # conditions$SaveTo <- paste0(saveto_2, "/")
  # 
  # solution_bv  <- tgas_calc_OM_seq_2(tgas_, src_ = conditions$Src,
  #                                    z_lim = conditions$Z, e_px = conditions$e_Px, bv = conditions$BV, Mg = conditions$MG,
  #                                    px_type = "DIST", distance = distance_,
  #                                    #save = conditions$SaveTo,
  #                                    type = conditions$KinModelType, model = conditions$KinModel,
  #                                    dist_type = conditions$Dist_Type, use = conditions$use,
  #                                    g_b = conditions$g_B)
  # solution_bv$Conditions <- conditions
  # tgas_write_conditions(conditions)
  # 
  # solution_bv <- tgas_process_solution(solution_bv)
  # solutions_bv$MS_ALL_rMoExp1 <- solution_bv
  # solutions_bv$MS_ALL_rMoExp1$Name <- paste("APASS-c")
  # 
  # gc()
  # ----------------------------------------------------------
#  cat("stage 4 \n\n")
# 
#   conditions$Dist_Type <- "TGAS_PX"
#   tgas_ <- tgas_calc_LClass(tgas_, dist_ = filter_dist, ph = "APASS")
# 
#   saveto_ <- paste0("solutions/solution_", name, "_", filter_dist, "_OM", conditions$KinModel)
#   if (!dir.exists(saveto_))
#     dir.create(saveto_)
# 
#   saveto_2 <- paste0(saveto_, "/Rpi-APASS")
#   if (!dir.exists(saveto_2))
#     dir.create(saveto_2)
# 
#   conditions$SaveTo <- paste0(saveto_2, "/")
# 
#   solution_bv  <- tgas_calc_OM_seq_2(tgas_, src_ = conditions$Src,
#                                      z_lim = conditions$Z, e_px = conditions$e_Px, bv = conditions$BV, Mg = conditions$MG,
#                                      px_type = "DIST", distance = distance_,
#                                      #save = conditions$SaveTo,
#                                      type = conditions$KinModelType, model = conditions$KinModel,
#                                      dist_type = conditions$Dist_Type, use = conditions$use,
#                                      g_b = conditions$g_B)
#   solution_bv$Conditions <- conditions
#   tgas_write_conditions(conditions)
# 
#   solution_bv <- tgas_process_solution(solution_bv)
#   solutions_bv$MS_ALL_rMoExp2 <- solution_bv
#   solutions_bv$MS_ALL_rMoExp2$Name <- "Rpi-APASS"
# 
# 
#   gc()
  # # ----------------------------------------------------------      
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
  
  tgas_draw_all_OM_sol_comp(solutions = solutions_bv, 
                            ylims  = matrix(data = c(0, 20, 0, 35, 0, 15, -2, 10, -5, 2, -25, -5, -5, 5, -8, 5 , 5, 30, -10, 3, -10, 5), nrow = 2),
                            xlims = c(-0.5, 1.6, 0.1),
                            xpar = 9, 
                            xtitle = "B-V", 
                            saveto = paste0("solutions/", filter_dist,"_"))
  
  # tgas_draw_all_kinematic_comp(solutions_bv, src = "TGAS", saveto = "solutions/om_kin", x_par = 9, x_lim = c(-0.6, 1.3,0.3), x_title = "B-V", is_legend = FALSE, width = 3.7, height = 3.7)
  
  
  return(solutions_bv)
  
}















