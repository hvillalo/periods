#' multiple regression
#' 
#' Fit a periodic regression model to a time series by multiple linear regression.
#' 
#' @param x Vector containing the time series to be fitted by multiple linear regression. 
#'
#' @param t Vector of time. Must be monotonically increasing of the same length of x. 
#'  
#' @param periods Periods of the harmonics (normally found by \code{cyclicDescent}.
#' 
#' @param trend Logical. If \code{TRUE}, the linear trend is included in the model. 
#'
#'
#' @details This function can ...  
#'
#' @return object of class 'lm'
#'
#' @author Héctor Villalobos.   
#'
#' @seealso \code{\link{cyclicDescent}}.
#'
#' @examples
#' # load simulated data
#' data(sim)
#' # find periodicities
#' sim.cd <- cyclicDescent(x = sim)
#' # significant or desired periods
#' p <- sim.cd$harmonics$Period[1:4]
#' multiple linear regression fit
#' sim.fit <- lm.harmonics(x = sim, periods = p, trend=FALSE)
#' summary(sim.fit)
#' 
lm.harmonics <- 
  function(x, t=NULL, periods, trend=FALSE)
  { 
    per <- periods
    n <- length(x)
    if (missing(t)) t <- 1:n 
    x.mean <- mean(x)                   
    x <- x - x.mean
    datdf <- data.frame(t, x)
    
    # Build Model
    compharm <- paste("cos(2*pi/", per, "*t) + ", "sin(2*pi/", per, "*t) ", 
                      sep="")
    if ( trend == TRUE ) { # estimate linear trend?
      modl <- as.formula(paste("x ~ t + ", paste(compharm, collapse=" + ")))
    } else {
      modl <- as.formula(paste("x ~ 0 + ", paste(compharm, collapse=" + ")))
    }
    # Multiple Regression
    fit.lm <- lm(modl, data=datdf)
    fit.lm$fitted <- fit.lm$fitted + x.mean # restore mean
    fit.lm
  }