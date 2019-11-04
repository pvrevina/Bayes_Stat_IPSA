### LOADING DATA 

library(haven)
library(dplyr)
library(psych)
library(stargazer)
library(ggplot2)

euro <- read_sav("ESS8FR.sav") %>% dplyr:::select(cptppola, sclmeet, sclact,inprdsc,health, hlthhmp, crmvct, 
                                                  aesfdrk, hincfel)
euro <- na.omit(euro)
summary(euro)

library(broom)
library(knitr)
library(bayesplot)
library(MCMCpack)
library(rstanarm)
euro$cptppola <- factor(euro$cptppola, levels = c(1, 2, 3, 4, 5))

### FITTING MODEL 
rst_polr <- stan_polr(cptppola ~ health + 
                        sclact + inprdsc + 
                        crmvct + aesfdrk + 
                        hincfel, data = euro,
                      prior = NULL, diagnostic_file = NULL)

kable(tidyMCMC(rst_polr, 
               conf.int = T, 
               conf.method = "quantile"), 
      digits = 3)

exp(coef(rst_polr))

mcmc_areas(as.matrix(rst_polr), pars = c("health", 
                                         "sclact",
                                         "inprdsc", 
                                         "crmvct",
                                         "aesfdrk", 
                                         "hincfel"))



mcmc_areas(as.matrix(rst_reduce_inprdsc), pars = c("health", 
                                         "sclact",
                                         "crmvct",
                                         "aesfdrk", 
                                         "hincfel"))


## CONVERGENCE DIAGNOSTICS 
library(ggmcmc)
library(ggthemes)
d_gr <- ggs(rst_polr)

ggs_traceplot(d_gr[d_gr$Parameter == "beta[1]",])
ggs_traceplot(d_gr[d_gr$Parameter == "beta[2]",])
ggs_traceplot(d_gr[d_gr$Parameter == "beta[3]",])
ggs_traceplot(d_gr[d_gr$Parameter == "beta[4]",])
ggs_traceplot(d_gr[d_gr$Parameter == "beta[5]",])
ggs_traceplot(d_gr[d_gr$Parameter == "zeta[1]",])

ggs_traceplot(d_gr)

ggs_Rhat(d_gr) + xlab("R_hat")

unique(d_gr$Parameter)

### Autocorrelation diagnostics 
ggs_autocorrelation(g_dat)

# Stationarity diagnostics
## Stationarity: local estimates from different parts of the chain should not differ much

ggs_geweke(d_gr)

## PP CHECKS 
yrep <- posterior_predict(rst_polr, draws = 500) 
yrep.num <- matrix(NA, nrow = nrow(yrep), ncol = ncol(yrep))
for (i in 1:ncol(yrep)) {
  yrep.num[,i] = as.numeric(yrep[,i])
}

ppc_stat(as.numeric(euro$cptppola), yrep.num, stat = "sd") + 
  theme(legend.title=element_text(size=22),
        legend.text=element_text(size=18),
        axis.text=element_text(size=18))

ppc_stat(as.numeric(euro$cptppola), yrep.num, stat = "mean") + 
  theme(legend.title=element_text(size=22),
        legend.text=element_text(size=18),
        axis.text=element_text(size=18))

ppc_dens_overlay(as.numeric(euro$cptppola), yrep.num) +
  theme(legend.text = element_text(face="bold", 
                                   size=36)) +
  theme(axis.text.x = element_text(face="bold", 
                                   size=18)) +
  theme(axis.text.y = element_text(face="bold", 
                                   size=18))

ppc_hist(as.numeric(euro$cptppola), yrep.num[1:5, ]) +
  theme(legend.text = element_text(face="bold", 
                                   size=36)) +
  theme(axis.text.x = element_text(face="bold", 
                                   size=18)) +
  theme(axis.text.y = element_text(face="bold", 
                                   size=18))

pp_check(rst_polr, plotfun = "stat_2d", stat = c("mean", "sd")) + 
  theme(legend.title=element_text(size=22),
        legend.text=element_text(size=18),
        axis.text=element_text(size=18))

pp_check(rst_polr, plotfun='error_scatter_avg')+ 
  theme(text=element_text(size=18))

#### COMPARING MODELS 
library(loo)

rst_reduce_sclact <- stan_polr(cptppola ~ health + 
                        inprdsc +  
                        crmvct + aesfdrk + hincfel, data = euro,
                      prior = NULL)

rst_reduce_inprdsc <- stan_polr(cptppola ~ health + 
                        sclact +  
                        crmvct + aesfdrk + hincfel, data = euro,
                      prior = NULL)

loo1 <- loo(rst_polr)
print(loo1, digits = 3)

waic1 <- waic(rst_polr)                    
print(waic1, digits = 3) 

loo2 <- loo(rst_reduce_sclact)
print(loo2, digits = 3)

waic2 <- waic(rst_reduce_sclact)                    
print(waic2, digits = 3) 

print(loo:::compare(loo1, loo2), 
      digits = 3)
print(loo:::compare(waic1, waic2),
      digits =3)

loo3 <- loo(rst_reduce_inprdsc)
print(loo3, digits = 3)

waic3 <- waic(rst_reduce_inprdsc)                    
print(waic3, digits = 3) 

print(loo:::compare(loo1, loo3), 
      digits = 3)
print(loo:::compare(waic1, waic3),
      digits =3)

