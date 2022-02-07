write_conditions <- function(conditions)
{
  con <- file(paste0(conditions$SaveTo,"description.txt"), "w")
  cat("B-V = ", conditions$Date, "\n", file=con)
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
  cat("Photometry = ", conditions$Photometry, "\n", file = con )
  
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

export_solution_xls <- function(res)
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

export_solution_txt <- function(res)
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

export_solution_latex <- function(res)
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


export_physical_latex <- function(res)
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

export_solution <- function(solution_)
{
  write_conditions(solution_$Conditions)
  export_solution_xls(solution_)
  export_solution_txt(solution_)
  export_solution_latex(solution_)
  export_physical_latex(solution_)
}

export_all_solution <- function(solutions)
{
  for(i in 1:length(solutions))
  {
    export_solution(solutions[[i]])
  }
  
}


#------------------------------------------------------
# стоит сразу и сводные диаграммы (все параметры на одной диаграмме) для всех решений
# и отдельные диаграммы для каждого параметра (все решения для каждого параметра на одной диаграмме)
draw_all_kinematic <- function(solutions, src = "TGAS", saveto = "", xpar = 4, xlim = c(0, 4, 0.5), x_title = "")
{
  for(i in 1:length(solutions))
  {
    draw_kinematic(solutions[[i]], x_par = xpar)
  }
  
  draw_all_kinematic_comp(solutions, src, x_par = xpar, x_lim = xlim, saveto, x_title = x_title)
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
draw_all_OM_sol_comp <- function(solutions, ylims, xlims = c(0, 2.5, 0.5), xpar = 4, xtitle = "<r>, kpc", saveto = "")
{
  
  for (i in 1:ncol(solutions[[1]]$X))
  {
    g <-draw_OMParameter(solutions, 
                         parameter = i, 
                         y_lim = c(ylims[1, i], ylims[2, i], ylims[3, i]), 
                         x_lim = xlims, 
                         x_title = xtitle, 
                         x_par = xpar, 
                         title = colnames(solutions[[1]]$X)[i])
    
    ggsave(paste0(saveto, "OM_", colnames(solutions[[1]]$X)[i],".png"), plot = g, width = 10, height = 5)
    #ggsave(paste0(saveto, "OM_", colnames(solutions[[1]]$X)[i],".eps"), plot = g, width = 10, height = 5)
    
  }
}

#------------------------------------------------------
# строит диаграммы для всех параметров решения по двум солюшинам для сравнения

draw_all_OM_sol <- function(sol1, sol2, sol1_name, sol2_name, saveto = "")
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

draw_kinematic <- function (solution, x_par = 4)
{
  save <- paste0(solution$Conditions$SaveTo, "_", solution$Conditions$Src,"_")
  
  g <- draw_OM_Solar(solution, title = paste("Solar motion,", solution$Conditions$Src, " proper motions, Solar motion."), xpar = x_par, xtitle = "B-V, mag")
  ggsave(paste0(save, "Solar-R", ".png"), plot = g, width = 10, height = 10)
  ggsave(paste0(save, "Solar-R", ".eps"), plot = g, width = 10, height = 10)
  
  if (solution$Conditions$KinModel <=2 ){
    g <- draw_Oort(solution, title = paste("Oort`s parameters,", solution$Conditions$Src, " proper motions."), xpar = x_par, xtitle = "B-V, mag")
    ggsave(paste0(save, "OL-R", ".png"), plot = g, width = 10, height = 10)
    ggsave(paste0(save, "OL-R", ".eps"), plot = g, width = 10, height = 10)
  }
  
  if (solution$Conditions$KinModel == 1){
    g <- draw_OM(solution, title = paste("Ogorodnikov-Miln Model,", solution$Conditions$Src, " proper motions."), xpar = x_par, xtitle = "B-V, mag")
    ggsave(paste0(save, "OM-R", ".png"), plot = g, width = 10, height = 10)
    ggsave(paste0(save, "OM-R", ".eps"), plot = g, width = 10, height = 10)
  }
  
}

#----------------------------------------------------

## Строит отдельные диаграммы для каждого параметра со всем решениями на каждой диаграмме

## tgas_draw_all_kinematic_comp(solutions_bv, src = "TGAS", saveto = "solutions/bottlinger_kin_", x_par = 9, x_lim = c(-0.5, 1.2,0.2), x_title = "B-V", is_legend = FALSE, width = 4, height = 4)

## tgas_draw_all_kinematic_comp(solutions, src = "TGAS", saveto = "solutions/om_kin", x_par = 9, x_lim = c(-0.6, 1.5,0.3), x_title = "B-V", is_legend = FALSE, width = 3.7, height = 3.7)


draw_all_kinematic_comp <- function(solutions, src = "TGAS", saveto = "", 
                                         x_par = 4, x_lim = c(0, 4, 0.5), x_title = "<r>, kpc",
                                         is_legend = TRUE, 
                                         width = 10, height = 5)
{
  cat(saveto)
  save <- paste0(saveto, "/_")
  
  g <- draw_OMParameter(solutions,
                        parameter = 1,
                        x_par = x_par,
                        x_lim = x_lim, 
                        x_title = x_title,
                        y_lim = c(0, 25, 2.5),
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
                        y_lim = c(0, 30, 2.5),
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
                        y_lim = c(0, 20, 2.5),
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
                          y_lim = c(5, 20, 2.5), 
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
                          y_lim = c(-18, -12, 1), 
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
                            y_lim = c(-6, 2, 1), 
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
                            y_lim = c(-10, 8, 2), 
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
                            y_lim = c(-15, 10, 2.5))
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
                            y_lim = c(-60, 60, 10))
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


draw_physics <- function (solution, saveto = "", is_legend = TRUE, width = 10, height = 5, x_par = 4, x_lim = c(0, 4, 0.5), src = "", x_title = "")
{
  
  if (solution[[1]]$Conditions$KinModel == 3)
    return(0);
  
  save <- paste0(saveto, "/_")
  
  g <- draw_Physical(solution, 
                     x_par = x_par,
                     x_lim = x_lim, 
                     x_title = x_title,
                     y_lim = c(170, 280, 10),
                     is_legend = is_legend,
                     title = paste("Linear galactic velocity at Solar distance, ", src, " proper motions."))
  ggsave(paste0(save, "V-R", ".png"), plot = g, width = width, height = height)
  #ggsave(paste0(save, "V-R", ".eps"), plot = g, width = width, height = height)
  
  g <- draw_Physical(solution, 
                     parameter = 2,
                     title = paste("Galaxy rotation period, ", src, " proper motions."),
                     x_par = x_par,
                     x_lim = x_lim, 
                     x_title = x_title,
                     y_lim = c(170, 280, 10),
                     is_legend = is_legend,
                     y_title = "million years")
  ggsave(paste0(save, "Period-R", ".png"), plot = g, width = width, height = height)
  #ggsave(paste0(save, "Period-R", ".eps"), plot = g, width = width, height = height)
  
  #draw_GalRotationCurveTilt(solution)
  g <- draw_Physical(solution, 
                     parameter = 3,
                     title = paste("Galaxy rotation curve inclination, ", src, " proper motions."),
                     x_par = x_par,
                     x_lim = x_lim, 
                     x_title = x_title,
                     y_lim = c(-10, 15, 2),
                     is_legend = is_legend,
                     y_title = "km/s/kpc")
  ggsave(paste0(save, "S-R", ".png"), plot = g, width = width, height = height)
  #ggsave(paste0(save, "S-R", ".eps"), plot = g, width = width, height = height)
  
  #draw_GalF(solution)
  g <- draw_Physical(solution, 
                     parameter = 4,
                     title = paste("Epicyclic frequency to angular velocity, ", src, " proper motions."),
                     x_par = x_par,
                     x_lim = x_lim,
                     x_title = x_title,
                     y_lim = c(1.2, 2.0, 0.1),
                     is_legend = is_legend,
                     y_title = "")
  ggsave(paste0(save, "F-R", ".png"), plot = g, width = width, height = height)
  #ggsave(paste0(save, "F-R", ".eps"), plot = g, width = width, height = height)
  
  #draw_GalMass(solution)
  g <- draw_Physical(solution, 
                     parameter = 5,
                     title = paste("Galaxy mass inside Solar orbit, ", src, " proper motions."),
                     x_par = x_par,
                     x_lim = x_lim,
                     x_title = x_title,
                     y_lim = c(5.0e10, 12e10, 1e10),
                     is_legend = is_legend,
                     y_title = "Solar mass")
  ggsave(paste0(save, "M-R", ".png"), plot = g, width = width, height = height)
  #ggsave(paste0(save, "M-R", ".eps"), plot = g, width = width, height = height)
  
  #draw_ApexL(solution)
  g <- draw_Physical(solution, 
                     parameter = 6,
                     title = paste("Solar motion apex L, ", src, " proper motions."),
                     x_par = x_par,
                     x_lim = x_lim, 
                     x_title = x_title,
                     y_lim = c(30, 90, 5),
                     is_legend = is_legend,
                     y_title = "degree")
  ggsave(paste0(save, "ApexL-R", ".png"), plot = g, width = width, height = height)
  #ggsave(paste0(save, "ApexL-R", ".eps"), plot = g, width = width, height = height)
  
  #draw_ApexB(solution)
  g <- draw_Physical(solution, 
                     parameter = 7,
                     title = paste("Solar motion apex B, ", src, " proper motions."),
                     x_par = x_par,
                     x_lim = x_lim, 
                     x_title = x_title,
                     y_lim = c(10, 35, 2),
                     is_legend = is_legend,
                     y_title = "degree")
  ggsave(paste0(save, "ApexB-R", ".png"), plot = g, width = width, height = height)
  #ggsave(paste0(save, "ApexB-R", ".eps"), plot = g, width = width, height = height)
  
  #draw_SolarV(solution)
  g <- draw_Physical(solution, 
                     parameter = 8,
                     title = paste("Solar velocity, ", src, " proper motions."),
                     x_par = x_par,
                     x_lim = x_lim, 
                     x_title = x_title,
                     y_lim = c(10, 40, 5),
                     is_legend = is_legend,
                     y_title = "km/s")
  ggsave(paste0(save, "SolarVr-R", ".png"), plot = g, width = width, height = height)
  #ggsave(paste0(save, "SolarVr-R", ".eps"), plot = g, width = width, height = height)
}


## tgas_draw_physics_BV(solutions_bv, src = "TGAS", saveto = "solutions/bottlinger_physycs_", x_lim = c(-0.5, 1.2, 0.2), is_legend = FALSE, width = 4, height = 4)

draw_physics_BV <- function (solution, src = "TGAS", saveto = "", x_lim = c(-1, 2, 0.1), is_legend = TRUE, width = 10, height = 5)
{
  
  if (solution[[1]]$Conditions$KinModel == 3)
    return(0);
  
  save <- paste0(saveto, src,"_")
  
  g <- draw_Physical(solution, 
                     x_par = 9, 
                     x_title = "B-V",
                     x_lim = x_lim,
                     data_x_lim = c(1, 13),
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
                     data_x_lim = c(2, 13),
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
                     y_lim = c(-5, 40, 5),
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

draw_HR_facet <- function(solution, M_lim = c(10,-10), BV_lim = c(-1, 3))
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

calc_weighted_mean <- function(X, sX)
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


calc_solution_weighted <- function(solution_)
{
  wr <- calc_weighted_mean(solution_$X, solution_$S_X)
  solution_$wX <- wr$wX
  solution_$s_wX <- wr$s_wX
  wr <- calc_weighted_mean(solution_$Physical, solution_$s_Physical)
  solution_$wPhysical <- wr$wX
  solution_$s_wPhysical <- wr$s_wX
  return(solution_)
}


calc_all_weighted <- function(solutions)
{
  for(i in 1:length(solutions))
  {
    solutions[[i]] <- calc_solution_weighted(solutions[[i]]) 
    
    # wr <- tgas_calc_weighted_mean(solutions[[i]]$X, solutions[[i]]$S_X)
    # solutions[[i]]$wX <- wr$wX
    # solutions[[i]]$s_wX <- wr$s_wX
    # wr <- tgas_calc_weighted_mean(solutions[[i]]$Physical, solutions[[i]]$s_Physical)
    # solutions[[i]]$wPhysical <- wr$wX
    # solutions[[i]]$s_wPhysical <- wr$s_wX
  }
  return(solutions)
}

process_solution <- function(solution_)
{
  solution_ <- calc_physical_params(solution_) 
  solution_ <- calc_solution_weighted(solution_)
  export_solution_xls(solution_)
  return(solution_)
}


export_all <- function(solutions, saveto_)
{
  
  cat("Export solutions...", "\n")
  export_all_solution(solutions)
  cat("Export kinematics...", "\n")
  draw_all_kinematic(solutions, saveto_)
  cat("Export physics...", "\n")
  draw_physics(solutions, src, saveto_)
}

