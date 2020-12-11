#'---
#'title: RSF and SSF analysis of fisher 
#'author: "S. Muff, J. Signer, J. Fieberg"
#'date: "`r format(Sys.time(), '%d %B, %Y')`"
#'output:
#'  html_document:
#'    toc: yes
#'---

#+ include = FALSE
knitr::opts_chunk$set(warning = FALSE, message = FALSE, cache = TRUE)

#' ## Purpose
#' The purpose of this document is to illustrate how a simple RSF and SSF with random effects can be fitted to tracking data. We use a data set containing fisher locations from:
#' 
#' - LaPoint, S., Gallery, P., Wikelski, M. and Kays, R. (2013). Animal behavior, cost-based corridor models, and real corridors. Landscape Ecology, 28, 1615-1630.
#' 

#' ## Load libraries and prepare data
#+ echo=TRUE, message=FALSE, warning=FALSE
library(glmmTMB)
library(INLA)
library(tidyverse)
library(raster)
library(survival)
library(TwoStepCLogit)
library(amt)

#' First load the fisher data (this is the data as it was downloaded from movebank).
#' We simplified the fisher data slightly, by running the following code prior to upload the data (to save space). This code is not needed, unless you download the original data from: https://www.datarepository.movebank.org/handle/10255/move.330.
#+ eval=FALSE
dat <- read_csv("Martes pennanti LaPoint New York.csv") %>% 
  filter(!is.na(`location-lat`)) %>% 
  select(x = `location-long`, y = `location-lat`, 
         t = `timestamp`, id = `tag-local-identifier`) %>% 
  filter(id %in% c(1465, 1466, 1072, 1078, 1016, 1469))
write_csv(dat, "Fisher_analysis/fisher_data.csv")

#' Now lets start by reading in the simplified fisher data

dat <- read_csv("fisher_data.csv")

#' Include sex of each animal and create tracks with an appropriate coordinate reference system using the amt package
dat_all <- dat %>% nest(-id) 
dat_all$sex <- c("f", "f", "f", "m", "m", "m")
dat_all <- dat_all %>% 
  mutate(trk = map(data, function(d) {
    make_track(d, x, y, t, crs = sp::CRS("+init=epsg:4326")) %>% 
      transform_coords(sp::CRS("+init=epsg:5070"))
  }))

#' Summarize sampling rates, 
dat_all %>% mutate(sr = lapply(trk, summarize_sampling_rate)) %>% 
  select(id, sr) %>% unnest

#' 10 minutes seems to appropriate for all animals.
#' Resample the track to 10 minutes with a tolerance of 2 minutes.

dat1 <- dat_all %>% mutate(dat_clean = map(trk, ~ {
  .x %>% track_resample(rate = minutes(10), tolerance = seconds(120))
  }))

#' Read in the landuse raster and reclassify to two categories (wet forests and other).
landuse <- raster("landuse_study_area.tif")
wet_forests <- landuse %in% c(90, 95)
names(wet_forests) <- "forest"

#' # Resource Selection Functions (RSF)
#' 
#' ## Data development for RSF
#' 
#' Now start with an RSF by creating random points per animal and extract the covariates for the observed and random points.

dat_rsf <- dat1 %>% mutate(rp = map(dat_clean, ~ .x %>% random_points() %>% 
      extract_covariates(wet_forests))) %>% 
  select(id, rp) %>%  unnest()
#' Change id column, to 1:6
dat_rsf$id <- as.numeric(factor(dat_rsf$id))

#' Make response numeric (required for INLA)
dat_rsf$y <- as.numeric(dat_rsf$case_)

#' We use a weighted likelihood for to fit the RSF. To this end, we need to create a variable for the weights, where used points (`case_ = TRUE`) keep weight 1, and available points (`case_ = FALSE`) obtain a large weight $W$ (here $W=1000$):
#+ echo=TRUE, message=FALSE
dat_rsf$weight <- 1000^(1 - dat_rsf$case_)

#' ## Mixed RSFs
#' 
#' ### glmmTMB()
#' 
#' 
#' As explained in the manuscript (Section 3.4), we recommend to manually fix the variance of the random intercept at a large value. This can be done in glmmTMB() by first setting up the model, but do not yet fit it:
#+ echo=TRUE, message=FALSE,cache=TRUE
fisher.tmp <- glmmTMB(case_ ~ forest + (1|id) + (0 + forest |id) , family=binomial(), data = dat_rsf,
                         doFit=FALSE, weights = weight)


#' Then fix the standard deviation of the first random term, which is the `(1|id)` component  in the above model equation. We use $\sigma=10^3$, which corresponds to a variance of $10^6$:
#+ echo=TRUE, message=FALSE,cache=TRUE
fisher.tmp$parameters$theta[1] <- log(1e3)


#' We need to tell `glmmTMB` not to change the first entry of the vector of variances, and give all other variances another indicator to make sure they can be freely estimated:
#+ echo=TRUE, message=FALSE,cache=TRUE
fisher.tmp$mapArg <- list(theta=factor(c(NA, 1)))

#' Then fit the model and look at the results:
#+ echo=TRUE, message=FALSE, cache=TRUE 
fisher.rsf <- glmmTMB:::fitTMB(fisher.tmp)
summary(fisher.rsf)


#' ###  INLA 
#'
#' Let us now carry the analysis with random intercept $\mathsf{N}(0,\sigma_{id}^2)$ and fixed variance $\sigma_{id}^2=10^6$ using INLA. A peculiarity of INLA is that the same variable cannot be used more than once. So for ID we need to generate a new (but identical) variable
#+ echo=TRUE, message=FALSE
dat_rsf$id1 <-dat_rsf$id

#' For the fixed effects we use the INLA (default) priors $\beta \sim \mathsf{N}(0,\sigma_\beta^2)$ with $\sigma_\beta^2=10^4$. The precisions of the priors are thus set to:
#+ echo=TRUE, message=FALSE
prec.beta.forest  <- 1e-4  

#' We now store the INLA formula with the fixed effects `forest`, plus two random effects, namely one for the individual-specific intercept and one for the individual-specific slope for `forest`. Note that the precision (thus $1/\sigma^2$) for `id` is fixed (`fixed=TRUE`) at the value of $10^{-6}$ (thus the variance is fixed at $10^6$). The other precision is given a PC(1,0.05) prior:
#+ echo=TRUE, message=FALSE
formula.inla <- y ~  forest + 
  f(id, model="iid", hyper=list(theta = list(initial=log(1e-6),fixed=TRUE))) +
  f(id1,forest,values=1:6,model="iid",
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(1,0.05)))) 


#' The actual INLA call is then given as follows:
#+  echo=TRUE, message=FALSE, cache=TRUE
inla.setOption(enable.inla.argument.weights=TRUE)
fisher.inla  <- inla(formula.inla, family ="binomial", data=dat_rsf, weights=dat_rsf$weight,
                        control.fixed = list(
                          mean = 0,
                          prec = list(forest = prec.beta.forest)
                       )
)



#' The summary for the posterior distribution of the fixed effects is given as follows:
#+ echo=TRUE
fisher.inla$summary.fixed


#' Since variances are parameterized and treated as precisions, the summary of the respective posterior distributions is given for the precisions:
#+ echo=TRUE
fisher.inla$summary.hyperpar


#' Source R functions for calculating posterior means 
#' and medians of the precisions.
source("inla_emarginal.R")
source("inla_mmarginal.R")
inla_emarginal(fisher.inla)
inla_mmarginal(fisher.inla)



#' # Step-Selection Function (SSF)
#' 
#' ## Data development for step-selection function
#' 
#' First we need to prepare the data. We need to pair each observed point with 10 random points and extract the covariate value at the end point of each step.

#+ warning = FALSE
dat_ssf <- dat1 %>% 
  mutate(stps = map(dat_clean, ~ .x %>% steps_by_burst() %>% 
                      random_steps() %>% extract_covariates(wet_forests))) %>% 
  select(id, stps) %>% unnest() %>% 
  mutate(
    y = as.numeric(case_),
    id = as.numeric(factor(id)), 
    step_id = paste0(id, step_id_, sep = "-"))
dat_ssf


#' ## Mixed SSFs 

#' ### 2StepCLogit 

#' The two-step procedure with independent random effect (D="UN(1)"):
r.Twostep <-  Ts.estim(formula = y ~ forest + strata(step_id) + 
                     cluster(id), data = dat_ssf, random = ~ forest,
                   all.m.1=F, D="UN(1)") 

#' Slope estimates and standard errors
r.Twostep$beta
r.Twostep$se

#' Variance estimates
r.Twostep$D

#' ### glmmTMB

#' Now the same model using `glmmTMB()`. Note that we do not need an overall intercept in this model, because the stratum-specific intercepts are (almost) freely estimated due to the large, fixed variance. Again start to set up the model without fitting the model:
TMBStruc <- glmmTMB(y ~ -1 + forest + (1|step_id) + 
                     (0 + forest | id),
                   family=poisson, data = dat_ssf, doFit=FALSE) 

#' Set the value of the standard deviation of the first random effect (here (1|step_id)):
TMBStruc$parameters$theta[1] <- log(1e3) 

#' Tell glmmTMB not to change the first standard deviation, all other values are freely estimated (and are different from each other)
TMBStruc$mapArg <- list(theta=factor(c(NA,1)))

#' Fit the model and look at the summary:
glmm.TMB.random <- glmmTMB:::fitTMB(TMBStruc)
summary(glmm.TMB.random)

#' 95\% CIs for fixed and random effects (standard deviations) are obtained via the confint() function:
confint(glmm.TMB.random)


#' Note: It is currently a problem of glmmTMB that the confint() function only shows a table of the length equal to the number of parameters estimated. As the variance for str_ID was not estimated but fixed to 10^6, but is still listed, the last variance component (here the one for (0 + forest | id)) is not shown. This can be solved by moving the component (1|step_id) to the last position in the formula, and then replace the respective lines by 
TMBStruc$parameters$theta[1] = log(1e3) 
TMBStruc$mapArg = list(theta=factor(c(1, NA)))
#' in the above code.
#' 
#' ### INLA  
#'
#' Start by setting the precision for the priors of slope coefficients (here for `forest`). 
prec.beta.forest <- 1e-4 

#' In the model formula for INLA, we set the stratum-specific intercept variance to $10^6$ (or rather: the precision to $10^{−6}$) by fixing it (`fixed=T`) to an intitial value. The other precision is given a PC(1,0.05) prior:
formula.random <- y ~  -1 + forest +
  f(step_id,model="iid",hyper=list(theta = list(initial=log(1e-6),fixed=T))) +
  f(id, forest, values=1:6,model="iid",
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(1,0.05)))) 

#' Fit the model
fisher.inla.ssf <- inla(formula.random, family ="Poisson", data=dat_ssf, 
                      control.fixed = list(
                        mean = 0,
                        prec = list(default = prec.beta.forest)
                      )
)

#' The summary for the posterior distribution of the fixed effects:
fisher.inla.ssf$summary.fixed 

#' Since variances are parameterized and treated as precisions, the summary of the respective posterior distributions is given for the precisions:
fisher.inla.ssf$summary.hyperpar

#' Posterior mean and mode are obtained as
inla_emarginal(fisher.inla.ssf)
inla_mmarginal(fisher.inla.ssf)



#' ## Session Info
#'
devtools::session_info()