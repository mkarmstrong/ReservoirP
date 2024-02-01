

fsg721 <- function(x, smth = 7) {
  
  # 2nd order polynomial sg filter
  # windowing for smoothing = smth (default 7)
  sg <- signal::sgolay(p = 2, n = smth, m = 1)
  sig <- signal::filter(sg, x)
  #plot(sig)
  return(sig)
  
}

RootSpline1 <- function (x, y, y0 = 0, verbose = FALSE) {
  
  if (is.unsorted(x)) {
    ind <- order(x)
    x <- x[ind]; y <- y[ind]
  }
  z <- y - y0
  ## which piecewise linear segment crosses zero?
  k <- which(z[-1] * z[-length(z)] <= 0)
  ## analytical root finding
  xr <- x[k] - z[k] * (x[k + 1] - x[k]) / (z[k + 1] - z[k])
  ## make a plot?
  if (verbose) {
    plot(x, y, "l"); abline(h = y0, lty = 2)
    points(xr, rep.int(y0, length(xr)))
  }
  ## return roots
  return(xr)
}

dicrotic <- function(pw, plot = FALSE) {
  
  # Get derivatives
  dp1 <- fsg721(pw)
  dp2 <- fsg721(fsg721(pw))
  dp3 <- fsg721(fsg721(fsg721(pw)))
  
  
  # FIND DICROTIC DEPRESSION ------------------------------------------------
  
  # End index
  end <- length(pw)
  
  # End index without potential perturbation at end diastole  
  end2 <- end * .9
  
  # Isolate notch area with 1st derivatives
  nni <- which.min(dp1)
  
  # Dicrotic notch from local dp2 max
  dic <- which.max(dp2[nni:end2]) + nni - 1
  
  # plot(pw, type="l", lwd=2)
  # par(new=T)
  # plot(dp2, type='o',col="grey")
  # abline(v = dic, h = 0)
  
  
  # FIND DICROTIC PEAK ------------------------------------------------------
  
  end3 <- ((end - dic) * .6) + dic # 60% of diastolic duration
  
  # Dicrotic peak from min of 2nd derivative
  # works better for subtle peaks
  if(sum(dp2[dic:end3] < 0) < 1) {
    dia <- 9999
  } else {
    dia <- which.min(dp2[dic:end3]) + dic - 1
  }
  
  # plot(pw, type="l", lwd=2)
  # par(new=T)
  # plot(dp2, type='o',col="grey")
  # abline(v = c(dic, dia), h = 0)
  
  
  # Dicrotic peak from 0 crossing of 1st derivative
  # works better for very definable peaks
  if (pw[dia] > pw[dic] & !is.na(pw[dia])) {
    hold <- RootSpline1(1:(end - nni + 1), dp1[nni:end], verbose = F)
    dia <- hold[2] + nni - 1
  }
  
  # plot(pw, type="l", lwd=2)
  # par(new=T)
  # plot(dp1, type='o',col="grey")
  # abline(v = c(dic, dia), h = 0)
  
  
  # PLOTS -------------------------------------------------------------------
  
  if(isTRUE(plot)) {
    plot(pw, type = "l", lwd=2, ylab="BP (mmHg)")
    abline(v=c(dic, dia), col="grey", lty=3, lwd=2)
    mtext(c("Ed", "P3"), side = 3, at = c(dic,dia))
  }
  
  return(data.frame(dicrotic_notch = dic, 
                    dicrotic_peak = dia))
  
}

low.pass <- function(y, fq, do.plot = FALSE) {
  # Second order low pass filter
  # Removes high frequency components below fq
  # y = a numeric vector, typically a tree-ring series.
  # fq = a numeric vector giving frequency or period of the filter.
  # Rp = a numeric value giving the dB for the passband ripple.
  
  if (any(is.na(y))) stop("y contains NA")
  
  ## n = a numeric value giving the order of the filter. 
  ## Larger numbers create steeper fall off.
  n = 4
  
  if (any(fq>1)) {
    f <- 1/fq
    p <- fq
  } else {
    p <- 1/fq
    f <- fq
  }
  
  ## sort f in case it's passed in backwards
  f <- sort(f)
  
  filt <- signal::butter(
    n = n,
    W = f * 2,
    type = "low",
    plane = "z"
  )
  
  ## remove mean
  yAvg <- mean(y)
  y <- y - yAvg
  
  ## pad the data to twice the max period
  pad <- max(p) * 2
  ny <- length(y)
  ## pad the data
  yPad <- c(y[pad:1], y, y[ny:(ny - pad)])
  ## run the filter
  yFilt <- signal::filtfilt(filt, yPad)
  ## unpad the filtered data
  yFilt <- yFilt[(pad + 1):(ny + pad)]
  ## return with mean added back in
  filt.sig <- yFilt + yAvg
  
  if(isTRUE(do.plot)){
    ## plot results
    plot(filt.sig,
         type = "l",
         lwd = 2)
  }
  
  ## return filtered signal
  return(filt.sig)
  
}

Tintersect <- function(wf, plot = FALSE) {
  
  k1 <- wf - min(wf[1:which.max(wf)])
  xvar <- (1:length(k1) - 1)
  spl <- smooth.spline(k1 ~ xvar)
  newx <- which.max(diff(k1))
  pred0 <- predict(spl, x = newx, deriv = 0)
  pred1 <- predict(spl, x = newx, deriv = 1)
  yint <- pred0$y - (pred1$y * newx)
  xint <- (-yint / pred1$y)
  
  if(isTRUE(plot)) {
    plot(xvar, k1, ylim=c(min(k1)-20, max(k1)))
    abline(h=min(k1), col="red", lty=3)
    lines(spl, col="red") 
    lines(xvar, yint + pred1$y*xvar, col="green", lwd=2)
    points(pred0,col="red", pch=8, lwd=2) 
    points(xint, 0, col="red", pch=8, lwd=2) 
    abline(v=xint, lty=3)
  }
  
  return(xint)
  
}

rp_exp_style <- function(x) {
  
  # Set up varibles
  P <- low.pass(x, .1)
  p <- P[(Tintersect(P)-1):(tail(which(diff(P) < 0), n = 1))]
  t <- 0:(length(p)-1)/200
  dt <- t[2] - t[1]
  nn <- dicrotic(p)$dicrotic_notc
  pd <- p[nn:length(p)]
  td <- 0:(length(pd)-1)/200
  
  # Fit diastolic pressure --------------------------------------------------
  
  # Initial parameter estimates (pa is per Kottenburg-Assenmacher 2009 & Schipke 2003)
  P0 <- pd[1]                               # p @ dicrotic notch
  TC <- unname(coef(lm(log(1/pd) ~ td))[2]) # slope
  pa <- 25                                  # asymptote
  
  # Diastolic fitting
  fit <- minpack.lm::nlsLM(
    pd ~ pa + (p0-pa)*exp(-td/tc), # exponential with asymptote
    control = minpack.lm::nls.lm.control(maxiter = 200, ptol = 1e-6),
    start = list(p0 = P0,
                 tc = TC))
  
  # Get diastolic fits
  prd <- predict(fit, data.frame(td = td))
  prd_extend <- predict(fit, data.frame(td = (0:1000)/200)) # extend for 5 sec
  
  # Extract optimized values
  tc <- unname(coef(fit)["tc"]) # time constant
  p0 <- unname(coef(fit)["p0"]) # starting pressure
  
  # plot(prd_extend, type="l", lty=2, ylim = c(pa, max(prd)))
  # lines(pd, col=2, lwd=3); abline(h=pa, lty=2)
  
  
  # Derive reservoir pressure -----------------------------------------------
  
  # Set up
  RP <- integer(length(p))
  RP[1] <- p[1]
  k_data <- data.frame(kratio = numeric(40), 
                       kerror = numeric(40))
  
  # Find optimal kratio by minimizing the diff between p0 and RP[nn]
  for(kratio in seq(1, 40, 1)) { # Max kratio = 22 per Mitchell, could flag kratio > 22
    for (i in 2:nn) {
      RP[i] <- RP[i-1]+((p[i]-RP[i-1])*kratio/tc-(RP[i-1]-pa)*1/tc)*dt
    }
    k_data[kratio, 1] <- kratio
    k_data[kratio, 2] <- abs(p0 - RP[nn])
  }
  
  # Find first minimum error
  kratio <- k_data$kratio[which.min(c(TRUE, diff(k_data$kerror) < 0)) - 1]
  
  #plot(k_data$kerror ~ k_data$kratio, col=3, pch=19); abline(v=kratio)
  
  # Derive reservoir pressure with optimized values
  for (i in 2:length(p)) { # could be: i in 2:nn
    RP[i] <- RP[i-1]+((p[i]-RP[i-1])* kratio /tc-(RP[i-1]-pa)*1/tc)*dt
  }
  
  ks <- tc/kratio # systolic rate constant
  kd <- 1/tc      # diastolic time constant to rate constant
  
  # Replace reservoir fit in diastole with diastolic fit
  rp <- RP
  rp[nn:length(rp)] <- prd
  
  # Calculate xsp from p minus rp
  rp  <- low.pass(rp, 0.07)
  rpld  <- rp - p[1]
  xsp <- p - rp
  pp  <- p - p[1]
  
  #plot(pp,col="lightgrey",type="l", lwd=3, ylim=c(min(xsp), max(pp)))
  #lines(rpld, lwd=3, col=2); lines(rpld[1:nn], lwd=3, col=4)
  #lines(xsp, lwd=3, col=3)
  
  # Plot results ------------------------------------------------------------
  
  # Plot the original data, fitted curve, and derived reservoir pressure
  par(mfrow = c(2, 2),
      mar = c(3.5, 3.5, .5, .5),
      mgp = c(2, 1, 0))
  plot(pp,col="grey",type="l", lwd=3, ylim=c(min(xsp), max(pp)))
  lines(rpld, lwd=3, col=2); lines(rpld[1:nn], lwd=3, col=4)
  lines(xsp, lwd=3, col=7)
  plot(prd_extend, type="l", lty=2, ylim = c(pa-5, max(prd)))
  lines(pd, col=2, lwd=3); abline(h=pa)
  plot(k_data$kerror ~ k_data$kratio, type = "l", lwd=2, col=6); abline(v=kratio, lty=2)
  
  # Save results ------------------------------------------------------------
  
  out <- data.frame(
    rp_amp = max(rpld),
    ep_amp = max(xsp),
    rp_int = sum(rpld)/200,
    ep_int = sum(xsp)/200,
    sys_k  = ks,
    dia_k  = kd,
    crit_p = pa,
    kratio = kratio
  )
  
  # Round values in "out"
  num_cols <- unlist(lapply(out, is.numeric)) # identify numeric cols
  out[num_cols] <-  round(out[num_cols], 3)   # round numeric cols
  
  # Print values to console
  for(i in 1:length(out)){
    print(paste0(names(out[i]),": ", out[1,i]), quote = F)
  }
  
  return(out)
  
}
