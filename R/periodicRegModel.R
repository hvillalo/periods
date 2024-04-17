#' Periodic Regression Model
#' 
#' Build the periodic regression model according to given periods.
#' 
#' @param x Vector containing the time series to be fitted. 
#'
#' @param t Vector of time. Must be monotonically increasing of the same length of x. 
#'
#' @param periods Vector of periods, typically found with \code{cyclicDescent}. 
#'
#' @param center.x Logical. If \code{TRUE}, the vector x is centered to mean zero.
#'
#' @param trend Logical. If \code{TRUE}, the linear trend is included in the model. 
#'
#'
#' @details This function can ...  
#'
#' @return A list with two elements, a data frame with variables and the periodic regression model to be passed to \code{lm}.
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
#' model to fit with lm
#' perReg <- periodicRegModel(x = sim, periods = p, center.x = FALSE, trend = FALSE)
#' perReg
#' # model fit
#' summary(lm(perReg$model, data = perReg$data))
periodicRegModel <- 
  function(x, t = NULL, periods, center.x = TRUE, trend = FALSE)
  { 
    n <- length(x)
    if (missing(t)) {
      t <- 1:n 
    }
    if(center.x == TRUE){
      x <- x - mean(x)
    }
    df <- data.frame(t, x)
    
    # Build Model
    compHarm <- paste("cos(2*pi/", periods, "*t) + ", 
                      "sin(2*pi/", periods, "*t) ", sep = "")
    if (trend == TRUE) { # estimate linear trend?
      model <- as.formula(paste("x ~ t + ", paste(compHarm, collapse=" + ")))
    } else {
      model <- as.formula(paste("x ~ 0 + ", paste(compHarm, collapse=" + ")))
    }
    ans <- list(data = df, model = model)
    ans
  }
 