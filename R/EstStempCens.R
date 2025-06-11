#' ML estimation in spatio-temporal model with censored/missing responses
#'
#' Return the maximum likelihood estimates of the unknown parameters of spatio-temporal model with censored/missing responses.
#' The estimates are obtained using SAEM algorithm.
#' The function also computes the observed information matrix using the method developed by \insertCite{louis1982finding;textual}{StempCens}.
#' The types of censoring considered are left, right, interval or missing values.
#'
#' @param y a vector of responses.
#' @param x a matrix or vector of covariates.
#' @param cc a vector of censoring indicators. For each observation: \code{1} if censored/missing and \code{0} if non-censored/non-missing.
#' @param time a vector of time.
#' @param coord a matrix of coordinates of the spatial locations.
#' @param LI lower limit of detection. For each observation: if non-censored/non-missing \code{=y}, if left-censored/missing \code{=-Inf} or \code{=LOD} if right/interval-censored.
#' @param LS upper limit of detection. For each observation: if non-censored/non-missing \code{=y}, if right-censored/missing \code{=Inf} or \code{=LOD} if left/interval-censored.
#' @param init.phi initial value of the spatial scaling parameter.
#' @param init.rho initial value of the time scaling parameter.
#' @param init.tau2 initial value of the the nugget effect parameter.
#' @param tau2.fixo \code{TRUE} or \code{FALSE}. Indicate if the nugget effect (\eqn{\tau^2}) parameter must be fixed.
#' By default = \code{FALSE}.
#' @param type.Data type of the data: '\code{balanced}' for balanced data and '\code{unbalanced}' for unbalanced data. By default = \code{balanced}.
#' @param method optimization method used to estimate (\eqn{\phi}, \eqn{\rho} and \eqn{\tau^2}):
#' '\code{optim}' for the function \code{\link[stats]{optim}} and '\code{nlminb}' for the function \code{\link[stats]{nlminb}}.
#' By default = \code{nlminb}.
#' @param lower.lim,upper.lim vectors of lower and upper bounds for the optimization method.
#' If unspecified, the default is \code{c(0.01,-0.99,0.01)} for the lower bound and \code{c(30,0.99,20)} for the upper bound if tau2.fixo=\code{FALSE}.
#' @param type.S type of spatial correlation function: '\code{exponential}' for exponential, '\code{gaussian}' for gaussian,
#' '\code{matern}' for matern, '\code{pow.exp}' for power exponential and '\code{spherical}' for spherical function, respectively.
#' Default is \code{exponential} function.
#' @param kappa parameter for all spatial covariance functions. In the case of exponential, gaussian and spherical function \eqn{\kappa} is equal to zero.
#' For the power exponential function \eqn{\kappa} is a number between 0 and 2. For the matern correlation function is upper than 0.
#' @param IMatrix \code{TRUE} or \code{FALSE}. Indicate if the observed information matrix will be computed. By default = \code{TRUE}.
#' @param M number of Monte Carlo samples for stochastic approximation. By default = \code{20}.
#' @param perc percentage of burn-in on the Monte Carlo sample. By default = \code{0.25}.
#' @param MaxIter the maximum number of iterations of the SAEM algorithm. By default = \code{300}.
#' @param pc percentage of iterations of the SAEM algorithm with no memory. By default = \code{0.20}.
#' @param error the convergence maximum error. By default = \code{1e-6}.
#'
#' @details The spatio-temporal Gaussian model is giving by:
#'
#' \eqn{ Y(s_i,t_j)= \mu(s_i,t_j)+ Z(s_i,t_j) +  \epsilon(s_i,t_j),}
#'
#' where the deterministic term \eqn{\mu(s_i,t_j)} and the stochastic terms \eqn{Z(s_i,t_j)},
#' \eqn{\epsilon(s_i,t_j)} can depend on the observed spatio-temporal indexes for \eqn{Y(s_i,t_j)}.
#' We assume \eqn{Z} is normally distributed with zero-mean and covariance matrix \eqn{\Sigma_z = \sigma^2 \Omega_{\phi\rho}},
#' where \eqn{\sigma^2} is the partial sill, \eqn{\Omega_{\phi\rho}} is the spatio-temporal correlation matrix,\eqn{\phi}
#' and \eqn{\rho} are the spatial and time scaling parameters; \eqn{\epsilon(s_i,t_j)} is an independent and
#' identically distributed measurement error with \eqn{E[\epsilon(s_i,t_j)]=0}, variance
#' \eqn{Var[\epsilon(s_i,t_j)]=\tau^2} (the nugget effect) and \eqn{Cov[\epsilon(s_i,t_j), \epsilon(s_k,t_l)]=0}
#' for all \eqn{s_i =! s_k} or \eqn{t_j =! t_l}.
#'
#' In particular, we define \eqn{\mu(s_i,t_j)}, the mean of the stochastic process as
#'
#' \eqn{\mu(s_i,t_j)=\sum_{k=1}^{p} x_k(s_i,t_j)\beta_k,}
#'
#' where \eqn{x_1(s_i,t_j),..., x_p(s_i,t_j)} are known functions of \eqn{(s_i,t_j)}, and \eqn{\beta_1,...,\beta_p}
#' are unknown parameters to be estimated. Equivalently, in matrix notation, we have the spatio-temporal linear model as follows:
#'
#' \eqn{Y = X \beta + Z + \epsilon,}
#'
#' \eqn{Z ~ N(0,\sigma^2 \Omega_{\phi\rho}),}
#'
#' \eqn{\epsilon ~ N(0,\tau^2 I_m).}
#'
#' Therefore the spatio-temporal process, \eqn{Y}, has normal distribution with mean \eqn{E[Y]=X\beta} and
#' variance \eqn{\Sigma=\sigma^2\Omega_{\phi\rho}+\tau^2 I_m}. We assume that \eqn{\Sigma} is non-singular
#' and \eqn{X} has full rank.
#'
#' The estimation process was computed via SAEM algorithm initially proposed by \insertCite{delyon1999convergence;textual}{StempCens}.
#'
#' @return The function returns an object of class \code{Est.StempCens} which is a list given by:
#'
#'\describe{
#'   \item{\code{m.data}}{Returns a list with all data components given in input.}
#'   \item{\code{m.results}}{A list given by:}
#' }
#'   \item{theta}{final estimation of \eqn{\theta = (\beta, \sigma^2, \tau^2, \phi, \rho)}.}
#'   \item{Theta}{estimated parameters in all iterations, \eqn{\theta = (\beta, \sigma^2, \tau^2, \phi, \rho)}.}
#'   \item{beta}{estimated \eqn{\beta}.}
#'   \item{sigma2}{estimated \eqn{\sigma^2}.}
#'   \item{tau2}{estimated \eqn{\tau^2}.}
#'   \item{phi}{estimated \eqn{\phi}.}
#'   \item{rho}{estimated \eqn{\rho}.}
#'   \item{Eff.range}{estimated effective range.}
#'   \item{PsiInv}{estimated \eqn{\Psi^-1}, where \eqn{\Psi=\Sigma/\sigma^2}.}
#'   \item{Cov}{estimated \eqn{\Sigma}.}
#'   \item{SAEMy}{stochastic approximation of the first moment for the truncated normal distribution.}
#'   \item{SAEMyy}{stochastic approximation of the second moment for the truncated normal distribution.}
#'   \item{Hessian}{Hessian matrix, the negative of the conditional expected second derivative matrix given the observed values.}
#'   \item{Louis}{the observed information matrix using the Louis' method.}
#'   \item{loglik}{log likelihood for SAEM method.}
#'   \item{AIC}{Akaike information criteria.}
#'   \item{BIC}{Bayesian information criteria.}
#'   \item{AICcorr}{corrected AIC by the number of parameters.}
#'   \item{iteration}{number of iterations needed to convergence.}
#'
#' @author Katherine L. Valeriano, Victor H. Lachos and Larissa A. Matos
#'
#' @examples
#' \donttest{
#' set.seed(12345)
#' # Initial parameter values
#' beta <- c(-1,1.50)
#' phi  <- 5
#' rho  <- 0.45
#' tau2 <- 0.80
#' sigma2 <- 1.5
#' coord <- matrix(runif(10, 0, 10), ncol=2)
#' time  <- matrix(1:5, ncol=1)
#' x     <- cbind(rexp(25,2), rnorm(25,2,1))
#' # Data simulation
#' data  <- rnStempCens(x, time, coord, beta, phi, rho, tau2, sigma2,
#'                      type.S="matern", kappa=1.5, cens="left", pcens=0.20)
#' # Estimation
#' est_teste <- EstStempCens(data$yObs, x, data$ci, data$time, cbind(data$x.coord, data$y.coord),
#'                           data$lcl, data$ucl, init.phi=3.5, init.rho=0.5, init.tau2=0.7,
#'                           tau2.fixo=FALSE, kappa=1.5, type.S="matern", IMatrix=TRUE, M=20,
#'                           perc=0.25, MaxIter=300, pc=0.2)}
#'
#' @import mvtnorm
#' @importFrom tmvtnorm rtmvnorm
#' @importFrom MCMCglmm rtnorm
#' @importFrom ggplot2 ggplot aes geom_point labs geom_line geom_text
#' @importFrom grid pushViewport viewport grid.layout
#' @importFrom stats dnorm nlminb optim optimize sd dist uniroot
#' @importFrom utils setTxtProgressBar txtProgressBar
#' @importFrom Rdpack reprompt
#' @importFrom stats quantile
#'
#' @export CovarianceM
#' @export EstStempCens
#' @export DiagStempCens
#' @export PredStempCens
#' @export CrossStempCens
#' @export EffectiveRange
#' @export rnStempCens

EstStempCens = function(y,x,cc,time,coord,LI,LS,init.phi,init.rho,init.tau2,tau2.fixo=FALSE,
                        type.Data="balanced",method="nlminb",kappa=0,type.S="exponential",
                        IMatrix=TRUE,lower.lim=c(0.01,-0.99,0.01),upper.lim=c(30,0.99,20),M=20,
                        perc=0.25,MaxIter=300,pc=0.20,error=1e-6){

  #---------------------------------------------------------------------#
  #                              Validations                            #
  #---------------------------------------------------------------------#
  tempo <- as.data.frame(time)
  ly <- length(y)

  if (!inherits(coord, c("matrix", "data.frame"))) stop("coord must be a matrix or data.frame")

  if (!inherits(x, "matrix")) stop("x must be a matrix")

  if (!inherits(tempo, "matrix") & !inherits(tempo, "data.frame")) stop("time must be a matrix or data.frame")

  if (nrow(x)!=ly) stop("Number of rows in x is different to the number of observations in y")

  if (length(cc)!=ly) stop("Length of censored vector (cc) different to the length of y")

  if (nrow(tempo)!=ly) stop("Number of elements in time vector is different to the number of observations in y")
  if (nrow(coord)!=ly|ncol(coord)!=2) stop("Dimension of coordinate matrix is non-confortable")

  if (length(LI)!=ly) stop("Number of elements in LI is different to the number of observations in y")
  if (length(LS)!=ly) stop("Number of elements in LS is different to the number of observations in y")
  if (sum(LI > LS)>0) stop("LI must be lower or equal than LS")

  if (init.phi <= 0) stop("The spatial parameter can not be negative or equal to zero")
  if (init.rho > 1 | init.rho < (-1)) stop("The time scaling parameter can not be >1 or < -1")
  if (init.tau2 < 0) stop("The nugget effect can not be negative")

  if(!is.logical(tau2.fixo) | !is.logical(IMatrix)) stop("Parameters tau2.fixo and IMatrix must be logical (TRUE/FALSE) variables")

  if (tau2.fixo==FALSE){
    if(length(lower.lim)!=3 | length(upper.lim)!=3) stop("Number of elements to the limits different of 3")
    if(upper.lim[1]<=lower.lim[1] | lower.lim[1]<0) stop("non-confortable values to the limits of the spatial parameter")
    if(upper.lim[2]<=lower.lim[2] | lower.lim[2]<(-1) | upper.lim[2]>1) stop("non-confortable values to the limits of the temporal parameter")
    if(upper.lim[3]<=lower.lim[3] | lower.lim[3]<0) stop("non-confortable values to the limits of the nugget effect")
  }

  if (tau2.fixo==TRUE){
    if(length(lower.lim)!=2 | length(upper.lim)!=2) stop("Number of elements to the limits different of 2")
    if(upper.lim[1]<=lower.lim[1] | lower.lim[1]<0) stop("non-confortable values to the limits of the spatial parameter")
    if(upper.lim[2]<=lower.lim[2] | lower.lim[2]<(-1) | upper.lim[2]>1) stop("non-confortable values to the limits of the temporal parameter")
  }

  if (!method%in%c("optim","nlminb")) stop('method should be one of optim, nlminb')

  if (!type.Data%in%c("balanced","unbalanced")) stop('type.Data should be one of balanced, unbalanced')
  if (type.Data=="balanced"){
    n1 = length(unique(tempo))
    n2 = nrow(unique(cbind(coord)))
    if (nrow(x) == n1*n2) stop("Non-conformable dimensions in the input for balanced data")
  }

  if (!type.S%in%c("matern","gaussian","spherical","pow.exp","exponential")){
    stop('type.S should be one of matern, gaussian, spherical, pow.exp, exponential')}

  if (type.S=="pow.exp" & (kappa > 2| kappa<=0)) stop("kappa must be a real in (0,2]")
  if (type.S=="matern" & kappa <= 0) stop("kappa must be a real number in (0,Inf)")

  if (perc > 1 | perc<=0) stop("Invalid value for perc")
  if (pc > 1 | pc<=0) stop("Invalid value for pc")

  # Validation, arguments support of the function
  if(MaxIter <= 0 | MaxIter%%1 != 0) stop("MaxIter must be a positive integer")
  if(M <= 0 | M%%1 != 0) stop("M must be a positive integer")
  if(error <=0 | error > 1) stop("error must belong to the interval (0,1]")


  #---------------------------------------------------------------------#
  #                                Outputs                              #
  #---------------------------------------------------------------------#


  out.ST = saemST(y,x,cc,tempo,coord,LI,LS,init.phi,init.rho,init.tau2,tau2.fixo,type.Data,method,kappa,
                  type.S,IMatrix,lower.lim,upper.lim,M,perc,MaxIter,pc,error)

  betas  = round(out.ST$m.results$beta, 4)
  sigma2 = round(out.ST$m.results$sigma2, 4)
  tau2   = round(out.ST$m.results$tau2, 4)
  phi    = round(out.ST$m.results$phi, 4)
  rho    = round(out.ST$m.results$rho, 4)
  theta  = round(out.ST$m.results$theta, 4)

  range  = round(out.ST$m.results$Eff.range, 3)

  # Criteria
  Loglik = out.ST$m.results$loglik
  AIC    = out.ST$m.results$AIC
  AICc   = out.ST$m.results$AICcorr
  BIC    = out.ST$m.results$BIC

  if(IMatrix==TRUE & tau2.fixo==FALSE){

    MI_obs = sqrt(diag(solve(out.ST$m.results$Louis)))
    SEbeta   = MI_obs[1:length(betas)]
    SEsigma2 = MI_obs[length(betas)+1]
    SEtau2   = MI_obs[length(betas)+2]
    SEphi   = MI_obs[length(betas)+3]
    SErho   = MI_obs[length(betas)+4]
    SE       = round(c(SEbeta,SEsigma2,SEtau2,SEphi,SErho), 4)

    Estimates      = cbind(theta,SE)
    colx           = ncol(as.matrix(x))

    namesx = paste0('\u03b2',1)
    if(ncol(as.matrix(x))>1)
    {
      for(i in 2:ncol(as.matrix(x))){namesx = cbind(namesx, paste0('\u03b2',i))}
    }

    greeks = c(sigma='\u03c3\u00B2', tau='\u03C4\u00B2', phi='\u03D5', rho='\u03C1')
    dimnames(Estimates) = list(c(namesx[1:colx],paste0(greeks)),c("Estimates", "SE"))

    cat('\n')
    cat('---------------------------------------------------------------\n')
    cat('     Spatio-temporal models for censored/missing responses     \n')
    cat('---------------------------------------------------------------\n')
    print(Estimates)
    cat(paste('The effective range is',range,'spatial units.\n'))
     cat('--------------------------------------------------------------\n')
     cat('\r \n')
     criteriaPCR = c(Loglik, AIC, AICc, BIC)
     criteriaFin = round(as.matrix(criteriaPCR),digits=3)
     dimnames(criteriaFin) = list(c("Loglik.", "AIC", "AICcorr.", "BIC"),c("Value"))
     cat('\n')
     cat('Model selection criteria\n')
     cat('------------------------------------\n')
     print(criteriaFin)
     cat('------------------------------------\n')
     cat('\r \n')
  }

  if(IMatrix==TRUE & tau2.fixo==TRUE){

    MI_obs = sqrt(diag(solve(out.ST$m.results$Louis)))
    SEbeta   = MI_obs[1:length(betas)]
    SEsigma2 = MI_obs[length(betas)+1]
    SEphi   = MI_obs[length(betas)+2]
    SErho   = MI_obs[length(betas)+3]
    SE       = round(c(SEbeta,SEsigma2,SEphi,SErho), 4)

    Estimates      = cbind(theta[-(length(betas)+2)],SE)
    colx           = ncol(as.matrix(x))

    namesx = paste0('\u03b2',1)
    if(ncol(as.matrix(x))>1)
    {
      for(i in 2:ncol(as.matrix(x))){namesx = cbind(namesx, paste0('\u03b2',i))}
    }

    greeks = c(sigma='\u03c3\u00B2', phi='\u03D5', rho='\u03C1')
    dimnames(Estimates) = list(c(namesx[1:colx],paste0(greeks)),c("Estimates", "SE"))

    cat('\n')
    cat('---------------------------------------------------------------\n')
    cat('     Spatio-temporal models for censored/missing responses     \n')
    cat('---------------------------------------------------------------\n')
    print(Estimates)
    cat(paste('The effective range is',range,'spatial units.\n'))
    cat('--------------------------------------------------------------\n')
    cat('\r \n')
    criteriaPCR = c(Loglik, AIC, AICc, BIC)
    criteriaFin = round(as.matrix(criteriaPCR),digits=3)
    dimnames(criteriaFin) = list(c("Loglik.", "AIC", "AICcorr.", "BIC"),c("Value"))
    cat('\n')
    cat('Model selection criteria\n')
    cat('------------------------------------\n')
    print(criteriaFin)
    cat('------------------------------------\n')
    cat('\r \n')
  }

  if(IMatrix==FALSE){

    Estimates      = matrix(theta,ncol=1)
    colx           = ncol(as.matrix(x))

    namesx = paste0('\u03b2',1)
    if(ncol(as.matrix(x))>1)
    {
      for(i in 2:ncol(as.matrix(x))){namesx = cbind(namesx, paste0('\u03b2',i))}
    }

    greeks = c(sigma='\u03c3\u00B2', tau='\u03C4\u00B2', phi='\u03D5', rho='\u03C1')
    dimnames(Estimates) = list(c(namesx[1:colx],paste0(greeks)),c("Estimates"))

    cat('\n')
    cat('---------------------------------------------------------------\n')
    cat('     Spatio-temporal models for censored/missing responses     \n')
    cat('---------------------------------------------------------------\n')
    print(Estimates)
    cat(paste('The effective range is',range,'spatial units.\n'))
    cat('--------------------------------------------------------------\n')
    cat('\r \n')
    criteriaPCR = c(Loglik, AIC, AICc, BIC)
    criteriaFin = round(as.matrix(criteriaPCR),digits=3)
    dimnames(criteriaFin) = list(c("Loglik.", "AIC", "AICcorr.", "BIC"),c("Value"))
    cat('\n')
    cat('Model selection criteria\n')
    cat('------------------------------------\n')
    print(criteriaFin)
    cat('------------------------------------\n')
    cat('\r \n')
  }

   class(out.ST) <- "Est.StempCens"

   return(invisible(out.ST))

}

