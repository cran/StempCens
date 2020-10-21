#' Covariance matrix for spatio-temporal model
#'
#' It computes the spatio-temporal covariance matrix for balanced data, i.e., when we have the same temporal indexes
#' per location. To compute the spatial correlation it provides 5 functions: exponential, gaussian, matern,
#' spherical and power exponential. To compute the temporal correlation is used an autocorrelation function of an AR(1) process.
#'
#' @param phi value of the spatial scaling parameter.
#' @param rho value of the time scaling parameter.
#' @param tau2 value of the the nugget effect parameter.
#' @param sigma2 value of the partial sill.
#' @param distSpa \eqn{n x n} spatial distance matrix without considering repetitions.
#' @param disTemp \eqn{T x T} temporal distance matrix without considering repetitions.
#' @param type.S type of spatial correlation function: '\code{exponential}' for exponential, '\code{gaussian}' for gaussian,
#' '\code{matern}' for matern, '\code{pow.exp}' for power exponential and '\code{spherical}' for spherical function, respectively.
#' See the analytical form of these functions in \code{\link{EffectiveRange}}.
#' @param kappa parameter for all spatial covariance functions. In the case of exponential, gaussian and spherical function \eqn{\kappa} is equal to zero.
#' For the power exponential function \eqn{\kappa} is a number between 0 and 2. For the matern correlation function is upper than 0.
#'
#' @return The function returns the \eqn{nT x nT} spatio-temporal covariance matrix for balanced data.
#'
#' @author Katherine L. Valeriano, Victor H. Lachos and Larissa A. Matos
#'
#' @examples
#' # Initial parameter values
#' phi <- 5;     rho <- 0.45
#' tau2 <- 0.80; sigma2 <- 2
#' # Simulating data
#' n1 <- 10   # Number of spatial locations
#' n2 <- 5    # Number of temporal index
#' set.seed(1000)
#' x.co <- round(runif(n1,0,10),5)  # X coordinate
#' y.co <- round(runif(n1,0,10),5)  # Y coordinate
#' coord <- cbind(x.co,y.co)        # Cartesian coordinates without repetitions
#' time <- as.matrix(seq(1,n2))     # Time index without repetitions
#' # Covariance matrix
#' Ms <- as.matrix(dist(coord))     # Spatial distances
#' Mt <- as.matrix(dist(time))      # Temporal distances
#' Cov <- CovarianceM(phi,rho,tau2,sigma2,distSpa=Ms,disTemp=Mt,kappa=0,type.S="exponential")

CovarianceM = function(phi, rho, tau2, sigma2, distSpa, disTemp, kappa, type.S){

  if (phi <= 0){stop("The spatial parameter can not be negative or equal to zero")}

  if (rho>1 | rho<(-1)){stop("The time scaling parameter can not be >1 or < -1")}

  if (tau2 < 0){stop("The nugget effect can not be negative")}

  if (sigma2 < 0){stop("The partial sill can not be negative")}

  if (ncol(distSpa) != nrow(distSpa)) stop("Spatial distance matrix must be specified")

  if (ncol(disTemp) != nrow(disTemp)) stop("Temporal distance matrix must be specified")

  if (type.S!="matern" & type.S !="gaussian" & type.S != "spherical" & type.S != "pow.exp" & type.S != "exponential"){
    stop('type.S should be one of matern, gaussian, spherical, pow.exp, exponential')}

  if (type.S=="pow.exp" & (kappa > 2| kappa<=0)) stop("kappa must be a real in (0,2]")

  if (type.S=="matern" & kappa <= 0) stop("kappa must be a real number in (0,Inf)")

  #---------------------------------------------------------------------#
  #                                Outputs                              #
  #---------------------------------------------------------------------#

  Covariance <- CovM(phi, rho, tau2, sigma2, distSpa, disTemp, kappa, type.S)

  out.ST <- Covariance

 return(invisible(out.ST))

}
