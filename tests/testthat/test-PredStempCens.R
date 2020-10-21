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
#' # Initial parameter values
#' beta <- c(-1,1.50)
#' phi  <- 5;    rho <- 0.60
#' tau2 <- 0.80; sigma2 <- 2
#' # Simulating data
#' n1 <- 17   # Number of spatial locations
#' n2 <- 5    # Number of temporal index
#' set.seed(12345)
#' x.co <- round(runif(n1,0,10),9)   # X coordinate
#' y.co <- round(runif(n1,0,10),9)   # Y coordinate
#' coord <- cbind(x.co,y.co)         # Cartesian coordinates without repetitions
#' coord2 <- cbind(rep(x.co,each=n2),rep(y.co,each=n2)) # Cartesian coordinates with repetitions
#' time <- as.matrix(seq(1,n2))      # Time index without repetitions
#' time2 <- as.matrix(rep(time,n1))  # Time index with repetitions
#' x1 <- rexp(n1*n2,2)
#' x2 <- rnorm(n1*n2,2,1)
#' x  <- cbind(x1,x2)
#' media <- x%*%beta
#' # Covariance matrix
#' Ms  <- as.matrix(dist(coord))   # Spatial distances
#' Mt  <- as.matrix(dist(time))    # Temporal distances
#' Cov <- CovarianceM(phi,rho,tau2,sigma2,Ms,Mt,0.50,"pow.exp")
#' # Data
#' require(mvtnorm)
#' y <- as.vector(rmvnorm(1,mean=as.vector(media),sigma=Cov))
#' data <- data.frame(coord2,time2,y,x)
#' names(data) <- c("x.coord","y.coord","time","yObs","x1","x2")
#' # Splitting the dataset
#' local.est  <- coord[-c(4,13),]
#' data.est   <- data[data$x.coord%in%local.est[,1]&data$y.coord%in%local.est[,2],]
#' data.valid <- data[data$x.coord%in%coord[c(4,13),1]&data$y.coord%in%coord[c(4,13),2],]
#' # Censored
#' perc <- 0.10
#' y <- data.est$yObs
#' aa <- sort(y);  bb <- aa[1:(perc*nrow(data.est))]
#' cutof <- bb[perc*nrow(data.est)]
#' cc <- matrix(1,nrow(data.est),1)*(y<=cutof)
#' y[cc==1] <- cutof
#' data.est <- cbind(data.est[,-c(4,5,6)],y,cc,data.est[,c(5,6)])
#' names(data.est) <- c("x.coord","y.coord","time","yObs","censored","x1","x2")
#'
#' # Estimation
#' y  <- data.est$yObs
#' x  <- cbind(data.est$x1,data.est$x2)
#' cc <- data.est$censored
#' time2  <- matrix(data.est$time)
#' coord2 <- data.est[,1:2]
#' LI <- y; LI[cc==1] <- -Inf    # Left-censored
#' LS <- y
#' est_teste <- EstStempCens(y, x, cc, time2, coord2, LI, LS, init.phi=3.5,
#'                  init.rho=0.5, init.tau2=1, kappa=0.5, type.S="pow.exp",
#'                  IMatrix=FALSE, M=20, perc=0.25, MaxIter=300, pc=0.20)
#' class(est_teste)
#'
#' # Prediction
#' locPre <- data.valid[,1:2]
#' timePre <- matrix(data.valid$time)
#' xPre <- cbind(data.valid$x1,data.valid$x2)
#' pre_teste <- PredStempCens(est_teste, locPre, timePre, xPre)
#' library(ggplot2)
#' Model <- rep(c("y Observed","y Predicted"),each=10)
#' station <- rep(rep(c("Station 1", "Station 2"),each=5),times=2)
#' xcoord1 <- rep(seq(1:5),4)
#' ycoord1 <- c(data.valid$yObs,pre_teste$predValues)
#' data2 <- data.frame(Model,station,xcoord1,ycoord1)
#' ggplot(data=data2,aes(x=xcoord1,y=ycoord1)) + geom_line(aes(color=Model)) +
#' facet_wrap(station~.,nrow=2) + labs(x="",y="") + theme(legend.position="bottom")

PredStempCens = function(Est.StempCens, locPre, timePre, xPre){

  if(class(Est.StempCens)!="Est.StempCens") stop("An object of the class Est.StempCens must be provided")

  if(class(locPre)!="matrix" & class(locPre)!="data.frame") stop("locPre must be a matrix or data.frame")

  if(class(xPre)!="matrix") stop("xPre must be a matrix")

  if(class(timePre)!="matrix" & class(timePre)!="data.frame") stop("timePre must be a matrix or data.frame")

  if(nrow(locPre)!=nrow(timePre)) stop("The number of rows in locPre must be equal to the length of timePre")

  if(nrow(locPre)!=nrow(xPre)) stop("The number of rows in locPre must be equal to the number of rows in xPre")

  #---------------------------------------------------------------------#
  #                                Outputs                              #
  #---------------------------------------------------------------------#

  model = Est.StempCens
  loc.Pre = locPre
  time.Pre = timePre
  x.pre = xPre
  ypred <- PredictNewValues(model,loc.Pre,time.Pre,x.pre)

  out.ST <- ypred

 class(out.ST) <- "Pred.StempCens"

 return(invisible(out.ST))

}

