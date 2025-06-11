#' Cross-Validation in spatio-temporal model with censored/missing responses
#'
#' This function performs cross-validation in spatio-temporal model with censored/missing responses, which measure
#' the performance of the predictive model on new test dataset. The cross-validation method for assessing the
#' model performance is validation set approach (or data split).
#'
#' @param Pred.StempCens an object of class \code{Pred.StempCens} given as output by the \code{\link{PredStempCens}} function.
#' @param yObs.pre a vector of the observed responses, the test data.
#'
#' @return
#' \item{Bias}{bias prediction error.}
#' \item{Mspe}{mean squared prediction error.}
#' \item{Rmspe}{root mean squared prediction error.}
#' \item{Mae}{mean absolute error.}
#'
#' @author Katherine L. Valeriano, Victor H. Lachos and Larissa A. Matos
#'
#' @seealso \code{\link{EstStempCens}}, \code{\link{PredStempCens}}
#'
#' @examples
#' \donttest{
#' set.seed(400)
#' # Parameter values
#' beta <- c(-1,1.50)
#' phi  <- 5
#' rho  <- 0.6
#' tau2 <- 0.80
#' sigma2 <- 2
#'
#' # Coordinates and covariates
#' coords <- matrix(round(runif(14, 0, 10), 9), ncol=2) # Coordinates without repetitions
#' time <- as.matrix(seq(1, 5)) # Time index without repetitions
#' x    <- cbind(rexp(35, 2), rnorm(35, 2, 1))
#'
#' # Data
#' data <- rnStempCens(x, time, coords, beta, phi, rho, tau2, sigma2,
#'                     type.S="gaussian", kappa=0, cens="left", pcens=0.2)
#'
#' # Splitting the dataset
#' train <- data[1:32,]
#' test  <- data[33:35,]
#'
#' # Estimation
#' x  <- cbind(train$x1, train$x2)
#' cc <- train$ci
#'
#' est_teste <- EstStempCens(y=train$yObs, x, cc=train$ci, time=train$time, coord=train[, 1:2],
#'                           LI=train$lcl, LS=train$ucl, init.phi=3.5, init.rho=0.5,
#'                           init.tau2=1, type.Data="unbalanced", method="nlminb", kappa=0,
#'                           type.S="gaussian", IMatrix=TRUE, M=20, perc=0.25,
#'                           MaxIter=300, pc=0.20)
#' # Prediction
#' xPre <- cbind(test$x1, test$x2)
#' pre_teste <- PredStempCens(est_teste, test[,1:2], test$time, xPre)
#' class(pre_teste)
#'
#' # Cross-validation
#' cross_teste <- CrossStempCens(pre_teste, test$yObs)
#' cross_teste$Mspe # MSPE}

CrossStempCens = function(Pred.StempCens, yObs.pre){

  if(!inherits(Pred.StempCens, "Pred.StempCens")) stop("An object of the class Pred.StempCens must be provided")

  if (length(Pred.StempCens$predValues) != length(yObs.pre)) stop("non-confortable length of parameter yObs.pre")

  #---------------------------------------------------------------------#
  #                                Outputs                              #
  #---------------------------------------------------------------------#

  ypred = Pred.StempCens
  crossvalida <- CrossValidation(yObs.pre, ypred$predValues)

  crossvalida

}
