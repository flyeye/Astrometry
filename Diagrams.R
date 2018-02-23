
# построение диаграммы сравнения B-V по фотометрии Hipparcos`а и APASS

tgas_hip_5 <- tgas_hip_hip[tgas_hip_hip$LClass_apass==5,]
index <- ((tgas_hip_5$e_Px/tgas_hip_5$Px) < 0.1) & ((tgas_hip_5$parallax_error/tgas_hip_5$gPx)<0.1)
index <- ((tgas_hip_5$e_Px/tgas_hip_5$Px) > 0.1) & ((tgas_hip_5$parallax_error/tgas_hip_5$gPx)>0.1)

BV_ <- data.frame(cbind(tgas_hip_5$hBV, (tgas_hip_5$apasm_b-tgas_hip_5$apasm_v)))

BV_ <- data.frame(cbind(tgas_hip_5$hBV[index], (tgas_hip_5$apasm_b[index]-tgas_hip_5$apasm_v[index])))

colnames(BV_) <- c("HIP", "APASS")
ggplot() + geom_point(data = BV_, aes(x = BV_$APASS, y = BV_$HIP), alpha = 0.05, size = 0.5) + scale_x_continuous(limits = c(-0.5, 2.5)) + scale_y_continuous(limits = c(-0.5,2))


# построение диаграммы сравнения параллаксов по  Hipparcos`у и TGAS

PX_ <- data.frame(cbind(tgas_hip_5$Px, tgas_hip_5$gPx))

index <- ((tgas_hip_5$e_Px/tgas_hip_5$Px) < 0.1) & ((tgas_hip_5$parallax_error/tgas_hip_5$gPx)<0.1)
PX_ <- data.frame(cbind(tgas_hip_5$Px[index], tgas_hip_5$gPx[index]))

PX_ <- data.frame(cbind(tgas_hip$Px, tgas_hip$gPx))

index <- ((tgas_hip$e_Px/tgas_hip$Px) < 0.1) & ((tgas_hip$parallax_error/tgas_hip$gPx)<0.1)
PX_ <- data.frame(cbind(tgas_hip$Px[index], tgas_hip$gPx[index]))

colnames(PX_) <- c("HIP", "TGAS")
ggplot() + geom_point(data = PX_, aes(x = PX_$TGAS, y = PX_$HIP), alpha = 0.05, size = 0.5) + scale_x_continuous(limits = c(-0.5, 100)) + scale_y_continuous(limits = c(-0.5,100))