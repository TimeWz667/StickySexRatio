


source("Source/model.R")
load(file="Input/WHO_TB.rdata")







m <- get_ssr_model(function(ti, t0) sqrt(ti-t0), 2000)
pars <- m$new_parameters(sig_v=0, sig_h=0)
ss = simulate_ssr(m, pars, 2020, F)


ssr <- m
dat <- ss





calculate_likelihood <- function(x, ssr, dat) {
  x[4] <- max(min(x[4], 0.99), 0.01)
  pars0 <- ssr$new_parameters(lr_f=x[1], lr_m=x[2], betaT=x[3], kappa=x[4], sig_h=0)
  pars <- ssr$new_parameters(lr_f=x[1], lr_m=x[2], betaT=x[3], kappa=x[4])
  
  li <- 0
  d <- dat[1, ]
  h <- ssr$rh(d$Time, pars0)
  li <- li + ssr$dv_h(d, h, d$Time, pars)
  
  for (i in 2:nrow(dat)) {
    d <- dat[i, ]
    h <- ssr$rh_t(d$Time, pars0)
    li <- li + ssr$dv_h(d, h, d$Time, pars)
  }
  li
}


est <- nlm(function(x) - calculate_likelihood(x, ssr, dat), c(1.5, 0.25, 0.1, 0.5), hessian=T)
est$estimate





calculate_likelihood <- function(x, y) {
  x <- c(x, y, 0.1, 0.3)
  pars0 <- ssr$new_parameters(lr_f=x[1], lr_m=x[2], betaT=x[3], kappa=x[4], sig_h=0)
  pars <- ssr$new_parameters(lr_f=x[1], lr_m=x[2], betaT=x[3], kappa=x[4])
  
  li <- 0
  d <- dat[1, ]
  h <- ssr$rh(d$Time, pars0)
  li <- li + ssr$dv_h(d, h, d$Time, pars)
  
  for (i in 2:nrow(dat)) {
    d <- dat[i, ]
    h <- ssr$rh_t(d$Time, pars0)
    li <- li + ssr$dv_h(d, h, d$Time, pars)
  }
  li
}

fs <- seq(0.5, 20, by=0.5)
ms <- seq(0.5, 20, by=0.5)

lis <- matrix(0, length(fs), length(ms))
for (i in 1:length(fs)) {
  for (j in 1:length(ms)) {
    lis[i, j] <- calculate_likelihood(fs[i], ms[j])
  }
}

filled.contour(lis)
