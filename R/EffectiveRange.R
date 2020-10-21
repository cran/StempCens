#' Effective range for some spatial correlation functions
#'
#' It computes the effective range for an isotropic spatial correlation function, which is
#' commonly defined to be the distance from which the correlation becomes small, typically
#' below 0.05.
#'
#' @param cor effective correlation to check for. By default = 0.05.
#' @param phi spatial scaling parameter.
#' @param kappa smoothness parameter, required by the matern and power exponential functions. By default = 0.
#' @param Sp.model type of spatial correlation function: '\code{exponential}' for exponential,
#' '\code{gaussian}' for gaussian, '\code{matern}' for matern, '\code{pow.exp}' for power
#' exponential and '\code{spherical}' for spherical function, respectively. By default = \code{exponential}.
#'
#' @details The available isotropic spatial correlation functions are:
#'
#' \describe{
#' \item{\strong{Exponential}:}{\eqn{Corr(d) = exp{-d/\phi}},}
#'
#' \item{\strong{Gaussian}:}{\eqn{Corr(d) = exp{-(d/\phi)^2}},}
#'
#' \item{\strong{Matern}:}{\eqn{Corr(d) = 1/(2^(\kappa-1)\Gamma(\kappa))(d/\phi)^\kappa K_\kappa(d/\phi)},}
#'
#' \item{\strong{Power exponential}:}{\eqn{Corr(d) = exp{-(d/\phi)^\kappa}},}
#'
#' \item{\strong{Spherical}:}{\eqn{Corr(d) = 1 - 1.5 d/\phi + 0.5(d/\phi)^3},}
#'
#' where \eqn{d} is the Euclidean distance between two observations, \eqn{\phi} is the spatial scaling
#' parameter, \eqn{\Gamma(.)} is the gamma function, \eqn{\kappa} is the smoothness parameter and
#' \eqn{K_\kappa(.)} is the modified Bessel function of the second kind of order \eqn{\kappa}.
#'}
#'
#' @return The function returns the effective range, i.e., the approximate distance from which the
#' spatial correlation is lower than \code{cor}.
#'
#' @author Katherine L. Valeriano, Victor H. Lachos and Larissa A. Matos
#'
#' @examples
#' phi <- 164.60
#' range1 <- EffectiveRange(0.05, phi, kappa=0, Sp.model="exponential")
#' range2 <- EffectiveRange(0.05, phi, kappa=1, Sp.model="pow.exp")
#' # Note that these functions are equivalent.

EffectiveRange <- function(cor=0.05, phi, kappa=0, Sp.model="exponential"){

  if (cor<=0 | cor>=1) stop("The correlation must be a number in (0,1)")

  if (phi <= 0) stop("The spatial parameter can not be negative or equal to zero")

  if (Sp.model!="matern" & Sp.model !="gaussian" & Sp.model != "spherical" & Sp.model != "pow.exp" & Sp.model != "exponential"){
    stop('Sp.model should be one of matern, gaussian, spherical, pow.exp, exponential')}

  if (Sp.model=="pow.exp" & (kappa > 2| kappa<=0)) stop("kappa must be a real in (0,2]")

  if (Sp.model=="matern" & kappa <= 0) stop("kappa must be a real number in (0,Inf)")


  #---------------------------------------------------------------------#
  #                                Outputs                              #
  #---------------------------------------------------------------------#

  Ef.range <- Effective.range(cor, phi, kappa, Sp.model)

  out.ST <- Ef.range

  return(invisible(out.ST))

}
