# OM model through VSF 

tgas_sample <- filter_tgs_px(tgas[tgas$LClass_apass==5,], r_lim = c(450, 550))
tgas_sample <- tgas_calc_gpm(tgas_sample)
stars <- tgas_get_stars(tgas_sample)
res <- Calc_VSF_Coef(stars, J = 80)

oort <- matrix(0, nrow = 2, ncol = 12)
oort[1,] <- c(GetOM_UR_SphK(res$X, 1), GetOM_VR_SphK(res$X, 1), GetOM_WR_SphK(res$X, 1), GetOM_B_SphK(res$X, 1), GetOM_A_SphK(res$X, 1),GetOM_C_SphK(res$X, 1), GetOM_K_SphK(res$X, 1), GetOM_Gx_SphK(res$X, 1), GetOM_Gy_SphK(res$X, 1), 0, 0, 0)
oort[2,] <- c(GetOM_UR_SphK(res$X, 2), GetOM_VR_SphK(res$X, 2), GetOM_WR_SphK(res$X, 2), GetOM_B_SphK(res$X, 2), GetOM_A_SphK(res$X, 2), GetOM_C_SphK(res$X, 2), GetOM_K_SphK(res$X, 2), GetOM_Gx_SphK(res$X, 2), GetOM_Gy_SphK(res$X, 2), 0, 0, 0)
rownames(oort) <- c("1", "2")
colnames(oort) <- c("UR", "VR", "WR", "B", "A",  "C", "K", "Gx", "Gy", "", "", "")

stargazer(oort, type = "text")

oort_err <- matrix(0, nrow = 2, ncol = 12)
oort_err[1,] <- c(GetOM_UR_SphK(res$s_X, 1), GetOM_VR_SphK(res$s_X, 1), GetOM_WR_SphK(res$s_X, 1), GetOM_B_SphK(res$s_X, 1), GetOM_A_SphK(res$s_X, 1),GetOM_C_SphK(res$s_X, 1), GetOM_K_SphK(res$s_X, 1), GetOM_Gx_SphK(res$s_X, 1), GetOM_Gy_SphK(res$s_X, 1), 0, 0, 0)
oort_err[2,] <- c(GetOM_UR_SphK(res$s_X, 2), GetOM_VR_SphK(res$s_X, 2), GetOM_WR_SphK(res$s_X, 2), GetOM_B_SphK(res$s_X, 2), GetOM_A_SphK(res$s_X, 2), GetOM_C_SphK(res$s_X, 2), GetOM_K_SphK(res$s_X, 2), GetOM_Gx_SphK(res$s_X, 2), GetOM_Gy_SphK(res$s_X, 2), 0, 0, 0)
rownames(oort_err) <- c("1", "2")
colnames(oort_err) <- c("eU", "eV", "eW", "eB", "eA",  "eC", "eK", "eGx", "eGy", "", "", "")
stargazer(oort_err, type = "text")

stargazer(res_tgas$X, type = "text")


#s_oort <- c(res$s_X["eA"], res$s_X["eB"], res$s_X["eC"], res$s_X["eK"], res$s_X["eGx"], res$s_X["eGy"])

# -----------------------

stars <- matrix(0, ncol = 6, nrow = 100000)
stars[,1] <- runif(100000, 0, 2*pi)
stars[,2] <- runif(100000, -pi/2, pi/2)
stars[,3] <- 1
stars[,4] <- 0
stars[,5] <- 0
stars[,6] <- 0

a0 <- GetSphFuncK_matrix(80, stars[,1], stars[,2])
Sph_0 <- GetSphCoefDefault(81, r = 1)
b0 <- rowSums(t(t(a0)*Sph_0))
# b0 <- a0 %*% matrix(Sph_0, ncol = 1)
stars[,4] <- b0 / 4.74064
res <- TLS_Gen(a0, b0, mode = 2, ef = ncol(a0))
Sph_0 - res$X

a0 <- MakeOMCoef(stars, use = c(TRUE, FALSE, FALSE), model = 2, type = 1)
#  calculate B0
#OM_0 <- c(10.3, 15.2, -15,  15, -5, -10, 30) # type = 3, mu_l
#OM_0 <- c(10.3, 15.2, 10, -15,  15, -5, -3) # type = 1, mu_l + mu_b
OM_0 <- c(10.3, 15.2, -15,  15, -5) # type = 1, mu_l
b0 <- rowSums(t(t(a0)*OM_0)) 
stars[,4] <- b0[1:nrow(stars)] / 4.74064
#stars[,5] <- b0[(nrow(stars)+1):length(b0)]
res_tgas <- Calc_OM_Model(stars, use = c(TRUE, FALSE, FALSE), mode = 2, model = 2, type = 1)


res_tgas <- Calc_OM_Model(stars, use = c(TRUE, FALSE, FALSE), mode = 2, model = 2, type = 3)

