#' multiple regression
#' 
#' Read a specific variable from several global CMEMS files, for choosen area 
#' limits (and possibly a certain depth).
#' 
#' @param x.
#'
#' @param t. 
#'  
#' @param periods.
#' 
#' @param trend.
#'
#'
#' @details This function can ...  
#'
#' @return objeto de clase 'lm'
#'
#' @author Héctor Villalobos.   
#'
#' @seealso \code{\link{cyclicDescent}}.
#'
#' @examples
#' data(sim)
#' # descenso cíclico
#' sim.cd <- cyclicDescent(x = sim)
#' # periodos significativos
#' p <- sim.cd$harmonics$Period[1:4]
#' sim.fit <- lm.harmonics(x = sim, periods = p, trend=FALSE)
#' summary(sim.fit)
#' 
#' }
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