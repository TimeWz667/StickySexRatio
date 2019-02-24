get_ssr_model <- function(fn_dt, t0) {
  model <- list()
  class(model) <- "SSR"
  
  model$new_parameters <- function(lr_f=log(2), lr_m=log(1.5), betaT=0.1, kappa=0.3, 
                                   sig_v=1e-2, sig_h=1e-4, phi=1/3) {
    c(lr_f=lr_f, lr_m=lr_m, 
      betaT=betaT, rrT=exp(betaT),
      kappa=kappa,
      sig_h=sig_h,
      sig_v=sig_v,
      phi=phi)
  }
  
  model$dt <- fn_dt
  model$t0 <- t0
  
  model$dh <- function(h, ti, pars) {
    with(as.list(c(pars)), {
      dnorm(h$LogMuT, 0, sig_h, log=T)
    })
  }
  
  model$rh <- function(ti, pars) {
    with(as.list(c(pars)), {
      c(LogMuT=rnorm(1, 0, sig_h))
    })
  }
  
  model$dh_t <- function(h, ti, pars) {
    with(as.list(c(pars)), {
      log_mu <- betaT*fn_dt(ti, t0)
      dnorm(h$LogMuT, log_mu, sig_h, log=T)
    })
  }
  
  model$rh_t <- function(ti, pars) {
    with(as.list(c(pars)), {
      log_mu <- betaT*fn_dt(ti, t0)
      c(LogMuT=rnorm(1, log_mu, sig_h))
    })
  }
  
  model$dh_h <- function(h1, h0, ti, pars) {
    with(as.list(c(pars)), {
      log_mu <- h0$LogMuT + betaT*(fn_dt(ti, t0)-fn_dt(ti-1, t0))
      dnorm(h1$LogMuT, log_mu, sig_h, log=T)
    })
  }
  
  model$rh_h <- function(h0, ti, pars) {
    with(as.list(c(pars)), {
      log_mu <- h0$LogMuT + betaT*(fn_dt(ti, t0)-fn_dt(ti-1, t0))
      c(LogMuT=rnorm(1, log_mu, sig_h))
    })
  }
  
  find_delays <- function(h, pars) {
    with(as.list(c(pars, h)), {
      log_rate <- LogMuT
      log_rate_f <- log_rate + lr_f
      log_rate_m <- log_rate + lr_m
      list(DelayF=exp(-log_rate_f), DelayM=exp(-log_rate_m))
    })
  }
  
  model$find_delays <- find_delays
  
  model$dv_h <- function(v, h, ti, pars) {
    with(as.list(c(pars, h)), {
      delays <- find_delays(h, pars)
      l.f <- delays$DelayF
      l.m <- delays$DelayM
      mu <- log(1+phi*l.m) - log(1+phi*l.f) + log(kappa) - log(1-kappa)
      dnorm(v$LogFM, mu, sig_v, log=T)
    })
  }
  
  model$rv_h <- function(h, ti, pars) {
    with(as.list(c(pars, h)), {
      delays <- find_delays(h, pars)
      l.f <- delays$DelayF
      l.m <- delays$DelayM
      mu <- log(1+phi*l.m) - log(1+phi*l.f) + log(kappa) - log(1-kappa)
      
      lfm <- rnorm(1, mu, sig_v)
      prm <- 1/(1+exp(lfm))
      prf <- 1 - prm
      gm <- 1 / (1 + phi*l.m)
      gf <- 1 / (1 + phi*l.f)
          
      c(
        LogFM=lfm,
        PrF=prf,
        PrM=prm,
        GapF=gf,
        GapM=gm,
        Gap=gf*prf+gm*prm,
        DelayF=l.f*365,
        DelayM=l.m*365
      )
    })
  }
  
  model$dv_h_ty <- function(v, h, ti, pars) {
    with(as.list(c(pars, h)), {
      delays <- find_delays(h, pars)
      l.f <- delays$DelayF
      l.m <- delays$DelayM
      
      mu <- phi*(l.f-l.m) + phi^2*(l.m^2-l.f^2)/2 + log(kappa) - log(1-kappa)
      dnorm(v$LogFM, mu, sig_v, log=T)
    })
  }
  
  model$rv_h_ty <- function(h, ti, pars) {
    with(as.list(c(pars, h)), {
      delays <- find_delays(h, pars)
      l.f <- delays$DelayF
      l.m <- delays$DelayM
      mu <- phi*(l.f-l.m) + phi^2*(l.m^2-l.f^2)/2 + log(kappa) - log(1-kappa)
      
      lfm <- rnorm(1, mu, sig_v)
      prm <- 1/(1+exp(lfm))
      prf <- 1 - prm
      gm <- 1 / (1 + phi*l.m)
      gf <- 1 / (1 + phi*l.f)
      
      c(
        LogFM=lfm,
        PrF=prf,
        PrM=prm,
        GapF=gf,
        GapM=gm,
        Gap=gf*prf+gm*prm,
        DelayF=l.f*365,
        DelayM=l.m*365
      )
    })
  }
  
  return (model)
}


print.SSR <- function(model) {
  cat("Sticky Sex Ratio Model\n")
  cat("Start Time: ", model$t0, "\n")
  cat("Time function: \n")
  print(model$dt)
}


simulate_ssr <- function(ssr, pars, t1, taylor=F) {
  rv_h <- ifelse(taylor, ssr$rv_h_ty, ssr$rv_h)
  t0 <- ssr$t0
  h <- ssr$rh(t0, pars)
  v <- rv_h(h, t0, pars)
  dat <- data.frame(t(unlist(c(Time=t0, h, v))))
  
  for (ti in (t0+1):t1) {
    h <- ssr$rh_t(ti, pars)
    v <- rv_h(h, ti, pars)
    dat <- rbind(dat, unlist(c(Time=ti, h, v)))
  }
  
  dat
}





loglikelihood_trend <- function(log.fm, times, pars) {
  times <- times - min(times)
  with(as.list(pars), {
    lrf <- beta0 + betaT*times^2
    lrm <- beta0 + betaT*times^2 + betaM
    lambda.f <- exp(-lrf)
    lambda.m <- exp(-lrm)
    
    mu <- log(1+phi*lambda.m) - log(1+phi*lambda.f) + log(kappa) - log(1-kappa)
    
    li <- dnorm(log.fm, mu, sig, log=T)
    sum(li)
  })
}

mse_trend <- function(log.fm, times, pars) {
  times <- times - min(times)
  with(as.list(pars), {
    lrf <- beta0 + betaT*times
    lrm <- beta0 + betaT*times + betaM
    lambda.f <- exp(-lrf)
    lambda.m <- exp(-lrm)
    
    mu <- log(1+phi*lambda.m) - log(1+phi*lambda.f) + log(kappa) - log(1-kappa)
    
    mean((mu - log.fm)^2)
  })
}

