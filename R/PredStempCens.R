#' Prediction in spatio-temporal model with censored/missing responses
#'
#' This function performs spatio-temporal prediction in a set of new \code{S} spatial locations for fixed time points.
#'
#' @param Est.StempCens an object of class \code{Est.StempCens} given as output by the \code{\link{EstStempCens}} function.
#' @param locPre a matrix of coordinates for which prediction is performed.
#' @param timePre the time point vector for which prediction is performed.
#' @param xPre a matrix of covariates for which prediction is performed.
#'
#' @return The function returns an object of class \code{Pred.StempCens} which is a list given by:
#'
#' \item{predValues}{predicted values.}
#' \item{VarPred}{predicted covariance matrix.}
#'
#' @author Katherine L. Valeriano, Victor H. Lachos and Larissa A. Matos
#'
#' @seealso \code{\link{EstStempCens}}
#'
#' @examples
#' \donttest{
#' set.seed(12345)
#' # Parameter values
#' beta <- c(-1,1.50)
#' phi  <- 5
#' rho  <- 0.60
#' tau2 <- 0.80
#' sigma2 <- 2
#' coord <- matrix(runif(34, 0, 10), ncol=2)
#' time  <- matrix(1:5, ncol=1)
#' x     <- cbind(rexp(85,2), rnorm(85,2,1))
#'
#' # Data simulation
#' data  <- rnStempCens(x, time, coord, beta, phi, rho, tau2, sigma2,
#'                      type.S="pow.exp", kappa=0.5, cens="left", pcens=0.10)
#' # Splitting the dataset
#' train <- data[11:85,]
#' test  <- data[1:10,]
#'
#' # Estimation
#' x <- cbind(train$x1, train$x2)
#' coord2 <- cbind(train$x.coord, train$y.coord)
#' est_teste <- EstStempCens(train$yObs, x, train$ci, train$time, coord2, train$lcl,
#'                           train$ucl, init.phi=3.5, init.rho=0.5, init.tau2=1, kappa=0.5,
#'                           type.S="pow.exp", IMatrix=FALSE, M=20, perc=0.25, MaxIter=300,
#'                           pc=0.20)
#' # Prediction
#' xPre <- cbind(test$x1, test$x2)
#' pre_test <- PredStempCens(est_teste, test[,1:2], test$time, xPre)}

PredStempCens = function(Est.StempCens, locPre, timePre, xPre){

  if (!inherits(Est.StempCens, "Est.StempCens")) stop("An object of the class Est.StempCens must be provided")

  if (!inherits(locPre, c("matrix", "data.frame"))) stop("locPre must be a matrix or data.frame")
  if (ncol(locPre)!=2) stop("Non-conformable dimensions in locPre")
  if (!inherits(xPre, "matrix")) stop("xPre must be a matrix")

  timePre = matrix(timePre, ncol=1)
  locPre  = as.matrix(locPre)
  if (!inherits(timePre, c("matrix", "data.frame"))) stop("timePre must be a matrix or data.frame")

  if (nrow(locPre)!=nrow(timePre)) stop("The number of rows in locPre must be equal to the length of timePre")
  if (nrow(locPre)!=nrow(xPre)) stop("The number of rows in locPre must be equal to the number of rows in xPre")

  #---------------------------------------------------------------------#
  #                                Outputs                              #
  #---------------------------------------------------------------------#

  model  = Est.StempCens
  ypred  = PredictNewValues(model, locPre, timePre, xPre)
  out.ST = ypred

 class(out.ST) = "Pred.StempCens"

 return(invisible(out.ST))

}

