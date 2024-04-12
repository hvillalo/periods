#' Cyclic Descent
#' 
#' This function implements the cyclic descent with periodoc regression to find periodicities in a regular time series.
#' 
#' @param x Vector of values to be analyzed.
#'
#' @param t Vector of time (optional). If specified, it has to have the same length as \code{x} and must be monotonically increasing. 
#'  
#' @param trend Logical. If \code{TRUE}, the linear tendency is calculated and subtracted from the data.
#' 
#' @param ip Initial period to test.
#'
#' @param lp Final period to test.
#'
#' @param step Is the increment between periods. Defaults to 1, assuming that two successive data values are separated by one unit of time.
#'
#' @param hn Number of harmonics to estimate. If not specified, the function will stop when the increase in R-squared is not significative or when a maximum of 30 harmonics is reached. This behaviour is overrid when a value is provided.
#'  
#' @param neig By default, once an optimal period (op) is identified, it is excluded from subsequent searches, unless neig = −1. Positive integers will also exclude neighboring periods (op − neig and op + neig).
#'
#' @param exclude The periods given here will be excluded in the periods search.
#'
#' @param alpha The confidence level for the statistical test in the cyclic descent. Defaults to 0.05.
#'
#' @param rrss Logical. When \code{TRUE}, the reciprocal of the residual sum of squares in every step of the cyclic descent will be returned.
#'
#' @param plots Logical, when TRUE, the RRSS versus the tested periods and the cumulated harmonics fit in every step of the cyclic descent will be plotted.
#'
#' @references González-Rodríguez E, Villalobos H, Gómez-Muñoz VM and Ramos-Rodríguez A. 2015. Computational Method for Extracting and Modeling Periodicities in Time Series. Open Journal of Statistics 05(06): 604–17. \url{https://doi.org/10.4236/ojs.2015.56062}.
#'
#' @details This function can ...  
#'
#' @return A list with two elements, (1) 'harmonics' a six-column data frame giving the periods, amplitudes, phases, lags, RSS and R-squared. Note that after the first model, the following will include the new period referred to and those already found. (2) 'Stats' this data frame provides the F statistics, with the corresponding degrees of freedom and associated p-values for every model comparisons.
#'
#' @author Héctor Villalobos.   
#'
#' @seealso \code{\link{lm.harmonics}}.
#'
#' @examples
#' # load simulated data
#' data(sim)
#' # find periodicites
#' sim.cd <- cyclicDescent(x = sim)
#' sim.cd
#' 
cyclicDescent <- 
  function(x, t=NULL, trend=FALSE, ip=NULL, lp=NULL, step=NULL, hn=NULL,
           neig=0, exclude=NULL, alpha=0.05, rrss=FALSE, 
           plots=FALSE)
  {
    n <- length(x)                          # length of data series
    if (missing(t)){ 
      t <- 1:n                              # time vector
	}  
    t.unit <- unique(t[2:n]-t[1:(n-1)])
    if (length(t.unit) > 1){
      stop("only regular time series are allowed")
    }
    if (missing(ip) & t.unit != 1){
      ip <- t[3]                            # initial period to test
    } else { 
	  ip <- 3 
	} 
    if (missing(lp) & t.unit != 1){
      lp <- t[ceiling(n/2)]                 # final period to test
    } else { 
	  lp <- ceiling(n/2) 
	}
    if (missing(step))  
      step <- t.unit                        # step between test periods
    x.mean <- mean(x)                   
    xd <- x - x.mean                        # center data
    datdf <- data.frame(t=t, x=xd)
    # Removing Linear Trend
    if (trend == TRUE){ 
      xlt <- lm(xd ~ t)$fitted.values       # linear fit
      xd <- xd - xlt                        # detrended data
    } 
    
    # Periodic Regression function
    periodicRegression <- function(x, t, perio, cos.ot, sin.ot) {
      np <- length(perio)                   # number of periods to test
      vna <-rep(NA, np) 
      modp.pr <- data.frame(a1=vna, b1=vna, period=perio, rrss=vna)
      resi <- matrix(NA, nrow=length(x), ncol=np)
      for (i in 1:np) {
        mreg <- lm( x ~ 0 + cos.ot[, i] + sin.ot[, i] )
        modp.pr[i, 1:2] <- coef(mreg)
        resi[, i] <- residuals(mreg)
      }
      modp.pr[, 'rrss'] <- 1/colSums(resi^2)
      return(modp.pr)
    }
    
    # Cyclic Descent plot function   
    plot.cd <- function(pnb) {
      p1 <- as.numeric(formatC(harmonics[pnb, 'Period'], digits=3))
      Rsq <- as.numeric(formatC(harmonics[pnb, 'R.sq'], digits=3))
      if (pnb == 1){
        plev <- NA #; harmonics[pnb, 'p.value']
      } else {
        plev <- sta[pnb-1, 'p.value']
      }
      plev <- format.pval(plev, digits=max(3, getOption("digits") - 3))
      
      par(mfrow = c(2, 1), mar = c(4, 4, 2, 2))
      plot(resids[[pnb]], type="l", main=paste("op =", p1))
      plot(t, xd, type="o", col="grey30", main=substitute(paste(R^2, " = ", 
                    Rsq, " ; ", "p-value: ", plev), list(Rsq=Rsq, plev=plev)))
      lines(t, xcum, col="blue")
    }
    
    # Tables for results
    harmonics <- data.frame(a1=NA, b1=NA, Period=NA, MRRSS=NA, Amplitude=NA,
                            Phase=NA, Lag=NA, RSS=NA, R.sq=NA, F=NA, dfn=NA, 
                            dfd=NA, p.value=NA)
    sta <- data.frame(F=NA, dfn=NA, dfd=NA, p.value=NA)
    resids <- list()                        # rrss
    
    # First Model ***************************************************************
    SST <- sum((xd - mean(xd))^2)           # Total Sum of Squares
    xcum <- rep(0, n)                       # cumulated harmonic series
    per2eval <- seq(ip, lp, step)           # periods to evaluate

    if ( !missing(exclude) ){
      per2eval <- per2eval[!per2eval %in% exclude]
    }
    k1 <- 2                                 # initial number of paramaters
    
    # Periodic regression
    omega <- t(2*pi/per2eval)               # terms for periodic regressions
    cos.ot <- cos(t %*% omega)
    sin.ot <- sin(t %*% omega)
    
    rrss.pr <- periodicRegression(x=xd, t, perio=per2eval, cos.ot, sin.ot)  
    pars <- rrss.pr[which.max(rrss.pr[, 'rrss']), ]    # best fit model (MRRSS)
    harmonics[1, 1:4] <- pars[1:4]                     # save model coefs.
    a1 <- pars[1, 1]
    b1 <- pars[1, 2]
    p1 <- pars[1, 3]                              # first identified period 
    omega <- 2*pi/p1                               
    xcal <- a1 * cos(omega * t) + b1 * sin(omega * t) # fitted values
    xr <- xd - xcal                               # residual series
    xcum <- xcum + xcal                           # cumulated harmonics series
    RSS1 <- sum((xd - xcal) ^ 2)                  # sum of squares for model 1
    Rsq1 <- 1 - RSS1 / SST                        # R-squared for model 1
    # F-test for R-squared 
    dfn <- k1
    dfd <- n - k1 - 1
    Fm <- (Rsq1 / dfn) / ((1 - Rsq1) / dfd) 
    pFm <- pf(Fm, dfn, dfd, lower.tail=FALSE)
    
    harmonics[1, 'RSS'] <- RSS1
    harmonics[1, 'R.sq'] <- Rsq1 
    harmonics[1, 'F'] <- Fm
    harmonics[1, 'dfn'] <- dfn
    harmonics[1, 'dfd'] <- dfd
    harmonics[1, 'p.value'] <- pFm
    
    resids[[1]] <- rrss.pr[, 3:4]
    
    if (plots == TRUE) plot.cd(1)
    
    # Second and subsequent models... if any ***********************************
    # Initializing counters
    nn  <- 2                                # counter for harmonics
    rnn <- 2                                # counter for results in tables 
    if (is.null(hn)) {
      cstop <- 30                           # counter for stoping loop
      nsh <- cstop                          # maximum 30 harmonics when
      if (round(Rsq1, 2) == 1) {            #   not specified by the user
        nn <- cstop + 1                     # Stop if Rsq == 1
        nsh <- 1
      }
    } else {                                                        
      cstop <- hn                         
      nsh <- hn
    }
    
    while(nn <= cstop) {
      # removing periods found (and neighbors) 
      if (neig == -1) {          # an op is allowed to be eligible again
        per2eval <- per2eval
      } else {                     # ... or eliminated 
        idx <- which(per2eval >= (p1 - neig * step) & per2eval <= 
                     (p1 + neig * step))
        per2eval <- per2eval[-idx]
        cos.ot <- cos.ot[, -idx]
        sin.ot <- sin.ot[, -idx]
      }
      
      rrss.pr <- periodicRegression(t, x=xr, perio=per2eval, cos.ot, sin.ot) 
      pars <- rrss.pr[which.max(rrss.pr[, 'rrss']), ] # best fit model (MRRSS)
      harmonics[rnn, 1:4] <- pars[1:4]                # save model coefs.
      a1 <- pars[1, 1]
      b1 <- pars[1, 2]
      p1 <- pars[1, 3]                                # identified period 
      omega <- 2 * pi / p1                               
      xcal <- a1 * cos(omega * t) + b1 * sin(omega * t) # fitted values
      xr <- xr - xcal                                 # residual series
      xcum <- xcum + xcal                           # cumulated harmonics series
      RSS2 <- sum( (xd - xcum) ^ 2 )                # sum of squares for model 2
      Rsq2 <- 1 - RSS2 / SST                        # R-squared for model 2      
      k2 <- k1 + 2
      # F-test for R-squared 
      dfn <- k2
      dfd <- n - k2 - 1
      Fm <- (Rsq2 / dfn) / ((1 - Rsq2) / dfd)
      pFm <- pf(Fm, dfn, dfd, lower.tail=FALSE)
      
      harmonics[rnn, 'RSS'] <- RSS2
      harmonics[rnn, 'R.sq'] <- Rsq2 
      harmonics[rnn, 'F'] <- Fm
      harmonics[rnn, 'dfn'] <- dfn
      harmonics[rnn, 'dfd'] <- dfd
      harmonics[rnn, 'p.value'] <- pFm
      
      resids[[rnn]] <- rrss.pr[, 3:4]
      
      # F-test for successive models
      dfn <- k2 - k1                     # numerator d.f. for F-test
      dfd <- (n - k2 - 1)                # denominator d.f. for F-test
      RSS1 <- harmonics[rnn - 1, 'RSS']  # sum of squares for model 2
      Fc <- ((RSS1 - RSS2) / dfn) / (RSS2 / dfd)   # F-test 
      p.level <- pf(Fc, dfn, dfd, lower.tail=FALSE)  # associated p-level
      
      sta[rnn - 1, 'dfn'] <- dfn
      sta[rnn - 1, 'dfd'] <- dfd
      sta[rnn - 1, 'F'] <- Fc
      sta[rnn - 1, 'p.value'] <- p.level  
      
      if (plots == TRUE) plot.cd(rnn)
      
      # Controling the number of loop' tours
      if (missing(hn)) {   # When harmonic number (hn) is not specified
        if (round(Rsq2, 2) == 1) {  # Stop if Rsq == 1, again
          nsh <- nn
          nn <- cstop + 1
        } else {
          if (sta[rnn - 1, 'p.value'] < alpha) { # condition for continuing 
            k1 <- k1 + 2        # the number of parameters increases by 2
            nn <- nn + 1        # and the counters by one
            rnn <- nn
            Rsq1 <- Rsq2
          } else {              # when the condition above is not meet,
            nsh <- nn - 1
            nshm <- 0
            nn <- cstop + 1       # the counter exceeds 'cstop' by one    
          }
        }                       # and the loop stops                   
      } else {
        k1 <- k1 + 2            # When hn is specified, the number of 
        nn <- nn + 1            # paramaters and counters continue to 
        rnn <- nn               # increase until 'hn' is reached
        Rsq1 <- Rsq2
      }
    }

#    # Whole model plot 
#     if (plots %in%  c("a", "all", "l", "last")) {
#       R2 <- as.numeric(formatC(harmonics[nsh, 'R.sq'], digits=3))
#       pval <- format.pval(harmonics[nsh, 'p.value'], digits=max(3, 
#                           getOption("digits") - 3))
#       xe <- xcum
#       if ( exists("nshm") ) xe <- xcum - xcal
#       xe <- xe + x.mean
#       if ( trend == TRUE )  xe <- xe + xlt
#       main.t <- paste("Periods =", paste(formatC(harmonics$Period[1:nsh],
#                                                  digits=3), collapse=", "))
#       sub.t <- substitute(paste(R^2, " = ", R2, " ; ", "p-value: ", pval), 
#                           list( R2 = R2, pval = pval ))
#       ylims <- c( min( c(x, xe) ), max( c(x, xe) ) )
#       par( mfrow = c(1, 1) )
#       plot( t, x, type="o", ylab="", xlab = "time", col = "grey30", 
#             main = main.t, ylim = ylims ) 
#       lines(t, xe, lwd = 2, col = "blue" )
#       mtext(sub.t, side = 3)
#     } 
    
    # Formating tables and calculating Amplitude, Phase and Lag
    sta$p.value <- 
      format.pval(sta$p.value, digits = max(3, getOption("digits") - 3))  
    rownames(sta) <- paste("Models", 1:nrow(sta), "&", 2:(nrow(sta)+1), ":")
    
    rownames(harmonics) <- paste("Model", 1:nrow(harmonics), ":")
    harmonics$Amplitude <- sqrt(harmonics$a1^2 + harmonics$b1^2)
    harmonics$Phase <- atan2(harmonics$b1, harmonics$a1)
    harmonics$Lag <- harmonics$Period*harmonics$Phase/(2*pi)
    harmonics <- harmonics[, -c(1, 2, 4, 10:13)]
    
    if (rrss == TRUE){  
      CycD <- list(RRSS=resids, harmonics=harmonics, Stats=sta)
    } else { 
      CycD <- list(harmonics=harmonics, Stats=sta)
    }
    attr(CycD, "class") <- "periods"
    return(CycD)
  }