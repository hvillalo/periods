#' Plot a Periodic Regression Model
#' 
#' Read a specific variable from several global CMEMS files, for choosen area 
#' limits (and possibly a certain depth).
#' 
#' @param fit.lm
#'
#'
#' @details This function can ...  
#'
#' @return objeto a plot is produced
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
#' # build Model
#' perReg <- periodicRegModel(x = sim, periods = p, center.x = FALSE)
#' plot(perReg)
#'
plot_periodicReg <- function(fit.lm) {
  if (class(fit.lm) != "lm")
    stop ("object must be of class lm")
  t <- fit.lm$model$t
  x <- fit.lm$model$x
  x.hat <- fit.lm$fitted.values
  R2 <- as.numeric(formatC(summary(fit.lm)$r.squared, digits=3))
  stat <- summary(fit.lm)$fstatistic
  pval <- pf(stat[1], stat[2], stat[3], lower.tail = FALSE)
  pval <- format.pval(pval, digits = max(3, getOption("digits") - 3))
  
  main.t <- paste("Periods =", paste(formatC(op, digits=3), collapse=", "))
  sub.t <- substitute(paste(R^2, " = ", R2, " ; ", "p-value: ", pval), 
                      list( R2 = R2, pval = pval ))
  
  plot(t, x, type = "n", main = main.t)
  grid()
  lines(t, x, type = "b", col = "grey30")
  lines(x.hat, col = "blue")
  mtext(sub.t, side = 3)
}
