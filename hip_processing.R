

index <- (((tgas_hip$e_Px/tgas_hip$Px)<0.1) & ((tgas_hip$parallax_error/tgas_hip$gPx)<0.1) & (!is.na(tgas_hip$apasm_v)) & (!is.na(tgas_hip$apasm_b)) & ((!is.na(tgas_hip$e_hBV)) & (tgas_hip$e_hBV !=0 )) )
tgas_hip_px <- tgas_hip[index,]

tgas_hip_apass <- tgas_calc_LClass(tgas_hip_px, dist_ = "HIP", ph = "APASS")
tgas_hip_apass_5 <- tgas_hip_apass[tgas_hip_apass$LClass_apass==5,]


tgas_hip_hip <- tgas_calc_LClass(tgas_hip_px, dist_ = "HIP", ph = "HIP")
tgas_hip_5 <- tgas_hip_hip[tgas_hip_hip$LClass_apass==5,]


index <- (((tgas_hip_5$e_Px/tgas_hip_5$Px)<0.1) & ((tgas_hip_5$parallax_error/tgas_hip_5$gPx)<0.1))
solution_hip_01 <- tgas_make_OM_solutions_bv(tgas_hip_5[index,], filter_dist = "HIP", ph = "APASS", name = "HIP_HIP_APASS_01")