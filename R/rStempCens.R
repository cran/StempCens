#' Censored spatio-temporal data simulation
#'
#' It simulates balanced censored spatio-temporal data with a linear structure for an established censoring rate o limit of detection.
#'
#' @param x design matrix of dimensions \eqn{nt x p}.
#' @param time vector containing the unique time points at which the observations are made, of length \eqn{t}.
#' @param coords 2D unique spatial coordinates of dimension \eqn{n x 2}.
#' @param beta linear regression parameters.
#' @param phi value of the spatial scaling parameter.
#' @param rho value of the time scaling parameter.
#' @param tau2 value of the the nugget effect parameter.
#' @param sigma2 value of the partial sill.
#' @param type.S type of spatial correlation function: '\code{exponential}' for exponential, '\code{gaussian}' for gaussian,
#' '\code{matern}' for matern, '\code{pow.exp}' for power exponential and '\code{spherical}' for spherical function, respectively.
#' See the analytical expressions of these functions in \code{\link{EffectiveRange}}.
#' @param kappa parameter for all spatial covariance functions. In the case of exponential, gaussian and spherical function \eqn{\kappa} is equal to zero.
#' For the power exponential function \eqn{\kappa} is a number between 0 and 2. For the matern correlation function is upper than 0.
#' @param cens '\code{left}' or '\code{right}' censoring. By default='\code{left}'.
#' @param pcens desired censoring rate. By default=\code{0.10}.
#' @param lod desired detection limit for censored observations. By default=\code{NULL}.
#'
#' @return The function returns a data.frame containing the simulated data.
#'
#' @author Katherine L. Valeriano, Victor H. Lachos and Larissa A. Matos
#'
#' @examples
#' set.seed(1000)
#' # Initial parameter values
#' phi  <- 5
#' rho  <- 0.45
#' tau2 <- 0.80
#' sigma2 <- 2
#' beta   <- c(1, 2.5)
#' x <- cbind(1, rnorm(50))
#'
#' # Coordinates
#' coords <- matrix(runif(20, 0, 10), ncol=2) # Cartesian coordinates without repetitions
#' time   <- 1:5  # Time index without repetitions
#'
#' # Data simulation
#' data <- rnStempCens(x, time, coords, beta, phi, rho, tau2, sigma2,
#'                     type.S="exponential", cens="left", pcens=0.10)

rnStempCens = function(x, time, coords, beta, phi, rho, tau2, sigma2, type.S="exponential", kappa=0, cens="left", pcens=0.10, lod=NULL){

  x = as.matrix(x)
  if (!inherits(coords, c("matrix", "data.frame"))) stop("coords must be a matrix or data.frame")
  if (ncol(coords)!=2) stop("coords must be a matrix with 2 columns")

  if (!inherits(x, "matrix")) stop("x must be a matrix")
  time = matrix(time, ncol=1)
  coords = as.matrix(coords)

  nt = nrow(time) * nrow(coords)
  if (nrow(x) != nt) stop("x must have as many rows as the total number of spatio-temporal observations, i.e., rows in coords x number of time points.")

  if (length(beta) != ncol(x)) stop("Non-conformable dimensions between beta and x")
  if (phi <= 0){stop("The spatial parameter can not be negative or equal to zero")}
  if (rho>1 | rho<(-1)){stop("The time scaling parameter can not be >1 or < -1")}
  if (tau2 < 0){stop("The nugget effect can not be negative")}
  if (sigma2 < 0){stop("The partial sill can not be negative")}

  if (!type.S%in%c("matern","gaussian","spherical","pow.exp","exponential")){
    stop('type.S should be one of matern, gaussian, spherical, pow.exp, exponential')}

  if (type.S=="pow.exp" & (kappa > 2| kappa<=0)) stop("kappa must be a real in (0,2]")
  if (type.S=="matern" & kappa <= 0) stop("kappa must be a real number in (0,Inf)")

  if (is.null(cens)) stop("cens should be one of left or right")
  if (!cens%in%c("left", "right")) stop("cens should be one of left or right")

  if (is.null(pcens)) stop("pcens should be a real number in (0, 1)")
  if (pcens<0 | pcens>=1) stop("pcens should be a real number in [0, 1)")


  #---------------------------------------------------------------------#
  #                                Outputs                              #
  #---------------------------------------------------------------------#

  data <- randomStempCens(x, time, coords, beta, phi, rho, tau2, sigma2, kappa,
                          type.S, cens, pcens, lod)
  return(data)

}
