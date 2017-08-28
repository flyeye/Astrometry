om_mod_string = " model {
    for (i in 1:length(y)) {
      y[i] ~ dnorm(mu[i], prec)
      mu[i] = x[1]*U[i] + x[2]*V[i] + x[3]*W[i]  + x[4]*Wx[i] + x[5]*Wy[i]+ x[6]*B[i] + x[7]*M13[i] + x[8]*M23[i] + x[9]*A[i] + x[10]*C[i]+ x[11]*K[i]
    }

    for (i in 1:3) {
      x[i] ~ dnorm(10.0, 1.0/1.0e3)
    }

    x[4] ~ dnorm(0.0, 1.0/1.0e2)
    x[5] ~ dnorm(0.0, 1.0/1.0e2)
    x[6] ~ dnorm(-13.0, 1.0/1.0e1)
    x[7] ~ dnorm(0.0, 1.0/1.0e2)
    x[8] ~ dnorm(0.0, 1.0/1.0e2)
    x[9] ~ dnorm(14.0, 1.0/1.0e1)
    x[10] ~ dnorm(-2.0, 1.0/1.0e2)
    x[11] ~ dnorm(-2.0, 1.0/1.0e2)

    prec ~ dgamma(1.0/2.0, 1.0*1500.0/2.0)

   sig2 = 1.0 / prec
   sig = sqrt(sig2)
} "

data_jags = as.list(cbind(as.data.frame(a), b))
#data_jags$y = b


params = c("x", "sig")

om_inits = function() {
  inits = list("x"=c(rnorm(3,0.0,1000.0), 
                     rnorm(2,0.0,100.0), 
                     rnorm(1,-13,10.0), 
                     rnorm(2,0.0,100.0), 
                     rnorm(1,14.0,10.0), 
                     rnorm(2,-2.0,100.0)), 
               "prec"=rgamma(1,1.0,1.0))
}

om_mod = jags.model(textConnection(om_mod_string), data=data_jags, inits=om_inits, n.chains=3)

om_mod = jags.model(textConnection(om_mod_string), data=data_jags, n.chains=3)
update(om_mod, 1e3) # burn-in

om_mod_sim = coda.samples(model=om_mod,
                        variable.names=params,
                        n.iter=5e3)

om_mod_csim = as.mcmc(do.call(rbind, om_mod_sim))


## convergence diagnostics
plot(om_mod_sim)

gelman.diag(om_mod_sim)
autocorr.diag(om_mod_sim)
autocorr.plot(om_mod_sim)
effectiveSize(om_mod_sim)

summary(om_mod_csim)

## compute DIC
dic = dic.samples(om_mod, n.iter=1e3)

