#' harmonicos
#' 
#' Read a specific variable from several global CMEMS files, for choosen area 
#' limits (and possibly a certain depth).
#' 
#' @param fit.lm
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
#' harmonics(sim.fit)
#' 
harmonics <- 
  function(fit.lm) 
  {
    if (class(fit.lm) != "lm")
      stop ("object must be of class lm")
    coefs <- coef(fit.lm)
    ncoef <- length(coefs)
    lt <- is.na(coefs['t'])
    if (lt == FALSE)   
      nhar <- ncoef / 2
    else 
      nhar <- (ncoef - 2) / 2
    
    if (is.na(coefs['t'])){
      cnam <- names(coefs)[seq(1, ncoef, 2)]
      Period <- as.numeric(substr(cnam, 12, nchar(cnam)-5))
      a1b1 <- as.data.frame(matrix(coefs[1:ncoef], ncol=2, byrow=TRUE))
    } else {
      cnam <- names(coefs)[seq(3, ncoef, 2)]
      coefs.lt <- coefs[1:2]; names(coefs.lt) <- c("alpha", "betha")
      Period <- as.numeric(substr(cnam, 12, nchar(cnam)-5))
      a1b1 <- as.data.frame(matrix(coefs[3:ncoef], ncol=2, byrow=TRUE))
    }
    a1b1 <- cbind(Period, a1b1)
    names(a1b1) <- c("Period", "a1", "b1")
    a1b1$Amplitude <- sqrt(a1b1$a1^2 + a1b1$b1^2)
    a1b1$Phase <- atan2(a1b1$b1, a1b1$a1)
    a1b1$Lag <- a1b1$Period * a1b1$Phase / (2 * pi)
    a1b1 <- a1b1[order(a1b1$Amplitude, decreasing=TRUE), ]
    a1b1 <- a1b1[, -c(2, 3)]
    if (is.na(coefs['t']))    
      a1b1 <- list(cyclic_components=a1b1)
    else
      a1b1 <- list(linear_trend=coefs.lt, cyclic_components=a1b1)
     return(a1b1)
  }
