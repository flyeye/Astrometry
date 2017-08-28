

bvl <- c(bv[2, 1], bv[2, 2])
tgas_ <- tgas_calc_LClass(tgas, dist_ = "TGAS_PX")
tgas_ <- tgas_[tgas_$LClass_apass == 5,]
tgas_sample <- filter_tgs_px(tgas_, bv_lim = bvl);
tgas_sample <- tgas_calc_gpm(tgas_sample)
tgas_sample <- tgas_calc_distance(tgas_sample, dist_type)
tgas_sample <- tgas_sample %>% filter(R<50000)
stars <- tgas_get_stars(tgas_sample, src_)
a <- MakeOMCoef(stars = stars, use = use, model = model, type = type)
b <- PrepareOMRightSide(stars = stars, use = use)
names
colnames(a) <- names
colnames(b) <- c("y")


ol_mod_string = " model {
    for (i in 1:length(y)) {
      y[i] ~ dnorm(mu[i], prec)
      mu[i] = x[1]*U[i] + x[2]*V[i] + x[3]*W[i]  + x[4]*B[i] + x[5]*A[i]
    }

    for (i in 1:3) {
      x[i] ~ dnorm(10.0, 1.0/1.0e3)
    }

    x[4] ~ dnorm(-13.0, 1.0/1.0e1)
    x[5] ~ dnorm(14.0, 1.0/1.0e1)
    
    prec ~ dgamma(1.0/2.0, 1.0*1500.0/2.0)

   sig2 = 1.0 / prec
   sig = sqrt(sig2)
} "

data_jags = as.list(cbind(as.data.frame(a), b))
#data_jags$y = b


params = c("x", "sig")

inits_ol = function() {
  inits = list("x"=c(rnorm(3,0.0,1000.0), 
                     rnorm(1,-13,10.0), 
                     rnorm(1,14.0,10.0)), 
               "prec"=rgamma(1,1.0,1.0))
}

mod_ol = jags.model(textConnection(ol_mod_string), data=data_jags, inits=inits_ol, n.chains=3)

mod_ol = jags.model(textConnection(ol_mod_string), data=data_jags, n.chains=3)
update(mod_ol, 1e3) # burn-in

mod_ol_sim = coda.samples(model=mod_ol,
                        variable.names=params,
                        n.iter=5e3)

mod_ol_csim = as.mcmc(do.call(rbind, mod_ol_sim))


## convergence diagnostics
plot(mod_ol_sim)

gelman.diag(mod_ol_sim)
autocorr.diag(mod_ol_sim)
autocorr.plot(mod_ol_sim)
effectiveSize(mod_ol_sim)

summary(mod_ol_csim)

## compute DIC
dic = dic.samples(mod_ol, n.iter=1e3)

