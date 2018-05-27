

index <- (((tgas_hip$e_Px/tgas_hip$Px)<0.1) & ((tgas_hip$parallax_error/tgas_hip$gPx)<0.1) & (!is.na(tgas_hip$apasm_v)) & (!is.na(tgas_hip$apasm_b)) & ((!is.na(tgas_hip$e_hBV)) & (tgas_hip$e_hBV !=0 )) )
tgas_hip_px <- tgas_hip[index,]

index <- (tgas_hip_px$eapasm_b<0.25) & (tgas_hip_px$eapasm_v<0.25)
index <- ((tgas_hip_px$eapasm_b/tgas_hip_px$apasm_b)<0.1) & ((tgas_hip_px$eapasm_v/tgas_hip_px$apasm_v)<0.1)

tgas_hip_apass <- tgas_calc_LClass(tgas_hip_px[index,], dist_ = "HIP", ph = "APASS")

ggplot() + geom_histogram(data = tgas_hip_px[index,], aes(x = eapasm_v))

tgas_hip_apass <- tgas_calc_LClass(tgas_hip_px, dist_ = "HIP", ph = "APASS")
tgas_hip_apass_5 <- tgas_hip_apass[tgas_hip_apass$LClass_apass==5,]


tgas_hip_hip <- tgas_calc_LClass(tgas_hip, dist_ = "TGAS", ph = "HIP")
tgas_hip_5 <- tgas_hip_hip[tgas_hip_hip$LClass_apass==5,]


index <- (((tgas_hip_5$e_Px/tgas_hip_5$Px)<0.1) & ((tgas_hip_5$parallax_error/tgas_hip_5$gPx)<0.1))
solution_hip_01 <- tgas_make_OM_solutions_bv(tgas_hip_5[index,], filter_dist = "HIP", ph = "APASS", name = "HIP_HIP_APASS_01")


tgas_hip_hip <- tgas_hip
tgas_hip_hip$B_V <- NA
tgas_hip_hip$M <- NA
tgas_hip_hip$Mag <- NA
tgas_hip_hip$B_V <- tgas_hip_hip$hBV
tgas_hip_hip$Mag <- tgas_hip_hip$hMag
tgas_hip_hip <- tgas_calc_distance(tgas_hip_hip)
tgas_hip_hip <- tgas_calc_absolute_mag(tgas_hip_hip)
tgas_a <- tgas_hip_hip[(!is.na(tgas_hip_hip$B_V)) & (!is.na(tgas_hip_hip$M)),]
tgas_a <- mutate(tgas_a, LClass_apass = 0)
tgas_a$LClass_apass[is_main_sequence(tgas_a$B_V, tgas_a$M)] <- 5
tgas_a$LClass_apass[tgas_a$M<(-2)& (tgas_a$M>(-Inf))] <- 1
tgas_a$LClass_apass[(tgas_a$M>-1.5)&(tgas_a$M<2.5)&(tgas_a$B_V>0.8)&(tgas_a$B_V<2.5)] <- 3
nrow(tgas_hip_hip)
tgas_hip_hip$LClass_apass <- NULL
tgas_hip_hip <- tgas_hip_hip %>% left_join(tgas_a[ , names(tgas_a) %in% c("source_id", "LClass_apass")], by = "source_id")
tgas_hip_hip$LClass_apass[is.na(tgas_hip_hip$LClass_apass)] <- 0


#--------------------
  


tgas_ <- tgas %>% filter(!is.na(B_V) & !is.na(M) & (LClass_apass == 5) & ((B_V>1.0) & (B_V<Inf)) & (is.na(HIP)))

tgas_ <- tgas %>% filter(!is.na(B_V) & !is.na(M) & (LClass_apass == 5) & ((B_V>1.0) & (B_V<Inf)))


draw_tgas_hist_Z(tgas_, name_ = "NEW_MS_Z_M-stars")
ggplot(data = tgas_) + geom_histogram(aes(x = Mag), fill = "gray70", colour = "gray10") +
  scale_x_continuous(breaks=seq(7,15,by=1), minor_breaks=seq(7,15,by=0.25), limits = c(7,15)) +
  theme_bw()
ggsave(filename = "NEW_M_APASS_V.jpeg", width = 5, height = 5)


