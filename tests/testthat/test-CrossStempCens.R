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
#' # Initial parameter values
#' beta <- c(-1,1.50)
#' phi <- 5;     rho <- 0.6
#' tau2 <- 0.80; sigma2 <- 2
#' # Simulating data
#' n1 <- 7    # Number of spatial locations
#' n2 <- 5    # Number of temporal index
#' set.seed(400)
#' x.co <- round(runif(n1,0,10),9)    # X coordinate
#' y.co <- round(runif(n1,0,10),9)    # Y coordinate
#' coord <- cbind(x.co,y.co)          # Cartesian coordinates without repetitions
#' coord2 <- cbind(rep(x.co,each=n2),rep(y.co,each=n2)) # Cartesian coordinates with repetitions
#' time <- as.matrix(seq(1,n2))       # Time index without repetitions
#' time2 <- as.matrix(rep(time,n1))   # Time index with repetitions
#' x1 <- rexp(n1*n2,2)
#' x2 <- rnorm(n1*n2,2,1)
#' x  <- cbind(x1,x2)
#' media <- x%*%beta
#' # Covariance matrix
#' Ms <- as.matrix(dist(coord))   # Spatial distances
#' Mt <- as.matrix(dist(time))    # Temporal distances
#' Cov <- CovarianceM(phi,rho,tau2,sigma2,Ms,Mt,0,"gaussian")
#' # Data
#' require(mvtnorm)
#' y <- as.vector(rmvnorm(1,mean=as.vector(media),sigma=Cov))
#' data <- data.frame(coord2,time2,y,x)
#' names(data) <- c("x.coord","y.coord","time","yObs","x1","x2")
#' # Splitting the dataset
#' local.est  <- coord[c(1,2,4,5,6),]
#' data.est   <- data[data$x.coord%in%local.est[,1] & data$y.coord%in%local.est[,2],]
#' data.valid <- data[data$x.coord%in%coord[c(3,7),1] & data$y.coord%in%coord[c(3,7),2],]
#' # Censored
#' perc <- 0.20
#' y  <- data.est$yObs
#' aa <- sort(y);  bb <- aa[1:(perc*nrow(data.est))]
#' cutof <- bb[perc*nrow(data.est)]
#' cc <- matrix(1,nrow(data.est),1)*(y<=cutof)
#' y[cc==1] <- cutof
#' data.est <- data.frame(data.est[,-c(4,5,6)],y,cc,data.est[,c(5,6)])
#' names(data.est) <- c("x.coord","y.coord","time","yObs","censored","x1","x2")
#'
#' # Estimation
#' y  <- data.est$yObs
#' x  <- cbind(data.est$x1,data.est$x2)
#' cc <- data.est$censored
#' time2  <- as.data.frame(data.est$time)
#' coord2 <- data.est[,1:2]
#' LI <- y; LI[cc==1] <- -Inf     # Left-censored
#' LS <- y
#' est_teste <- EstStempCens(y, x, cc, time2, coord2, LI, LS, init.phi=3.5, init.rho=0.5,
#'                  init.tau2=1,type.Data="balanced", method="nlminb", kappa=0,
#'                  type.S="gaussian", IMatrix=FALSE, M=20, perc=0.25, MaxIter=300,
#'                  pc=0.20)
#' # Prediction
#' locPre <- data.valid[,1:2]
#' timePre <- as.data.frame(data.valid$time)
#' xPre <- cbind(data.valid$x1,data.valid$x2)
#' pre_teste <- PredStempCens(est_teste, locPre, timePre, xPre)
#' class(pre_teste)
#'
#' # Cross-validation
#' cross_teste <- CrossStempCens(pre_teste,data.valid$yObs)
#' cross_teste$Mspe # MSPE

CrossStempCens = function(Pred.StempCens, yObs.pre){

  if(class(Pred.StempCens)!="Pred.StempCens") stop("An object of the class Pred.StempCens must be provided")

  if (length(Pred.StempCens$predValues) != length(yObs.pre)) stop("non-confortable length of parameter yObs.pre")

  #---------------------------------------------------------------------#
  #                                Outputs                              #
  #---------------------------------------------------------------------#

  ypred = Pred.StempCens
  crossvalida <- CrossValidation(yObs.pre,ypred$predValues)

  crossvalida

}
