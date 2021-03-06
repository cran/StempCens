% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/DiagStempCens.R
\name{DiagStempCens}
\alias{DiagStempCens}
\title{Diagnostic in spatio-temporal model with censored/missing responses}
\usage{
DiagStempCens(Est.StempCens, type.diag = "individual", diag.plot = TRUE, ck)
}
\arguments{
\item{Est.StempCens}{an object of class \code{Est.StempCens} given as output by the \code{\link{EstStempCens}} function. In the \code{EstStempCens}function, \code{IMatrix} must be \code{TRUE}.}

\item{type.diag}{type of diagnostic: '\code{individual}' is related when one observation is deleted,
'\code{time}' is related when an entire time is deleted, '\code{location}' is related when an entire location is deleted and
'\code{all}' the three cases ('\code{individual}', '\code{time}' and '\code{location}').
By default \code{type.diag} is \code{individual}.}

\item{diag.plot}{\code{TRUE} or \code{FALSE}. It indicates if diagnostic plots must be showed. By default = \code{TRUE}.}

\item{ck}{the value for \code{ck} considered in the benchmark value for the index plot:
\eqn{mean(GD)+ck*sd(GD)}, where \eqn{GD} is the vector with all values of the diagnostic measures.}
}
\value{
The function returns a list with the diagnostic measures.

\describe{
  \item{If \code{type.diag == individual | time | location}:}{
  \code{GD} is a data.frame with the index value of the observation and the GD measure.}
  \item{If \code{type.diag == all}:}{
  \code{GDind} is a data.frame with the index value of the observation and the GD measure for individual.

  \code{GDtime} is a data.frame with the time index value and the GD measure for time.

  \code{GDloc} is a data.frame with the side index value and the GD measure for location.

  }
}
}
\description{
Return measures and graphics for diagnostic analysis in spatio-temporal model with censored/missing responses.
}
\details{
This function uses the case deletion approach to study the impact of deleting one or
more observations from the dataset on the parameters estimates, using the ideas of
\insertCite{cook1977detection;textual}{StempCens} and \insertCite{zhu2001case;textual}{StempCens}.
The measure is defined by

\eqn{GD_i(\theta*)=(\theta* - \theta*[i])'[-Q**(\theta|\theta*)](\theta* - \theta*[i]), i=1,....m,}

where \eqn{\theta*} is the estimate of \eqn{\theta} using the complete data, \eqn{\theta*[i]}
are the estimates obtained after deletion of the i-th observation (or group of observations) and
\eqn{Q**(\theta|\theta*)} is the Hessian matrix.

We can eliminate an observation, an entire location or an entire time index.
}
\examples{
\dontrun{
# Initial parameter values
beta <- c(-1,1.5)
phi <- 3;   rho <- 0.40
tau2 <- 1;  sigma2 <- 2
# Simulating data
n1 <- 5    # Number of spatial locations
n2 <- 5    # Number of temporal index
set.seed(98765)
x.co <- round(runif(n1,0,10),9)   # X coordinate
y.co <- round(runif(n1,0,10),9)   # Y coordinate
coord <- cbind(x.co,y.co)         # Cartesian coordinates without repetitions
coord2 <- cbind(rep(x.co,each=n2),rep(y.co,each=n2)) # Cartesian coordinates with repetitions
time <- as.matrix(seq(1,n2))      # Time index without repetitions
time2 <- as.matrix(rep(time,n1))  # Time index with repetitions
x1 <- rexp(n1*n2,2)
x2 <- rnorm(n1*n2,2,1)
x  <- cbind(x1,x2)
media <- x\%*\%beta
# Covariance matrix
Ms  <- as.matrix(dist(coord))  # Spatial distances
Mt  <- as.matrix(dist(time))   # Temporal distances
Cov <- CovarianceM(phi,rho,tau2,sigma2,Ms,Mt,0,"exponential")
# Data
require(mvtnorm)
y <- as.vector(rmvnorm(1,mean=as.vector(media),sigma=Cov))
perc <- 0.20
aa <- sort(y); bb <- aa[((1-perc)*n1*n2+1):(n1*n2)]; cutof <- bb[1]
cc <- matrix(1,(n1*n2),1)*(y>=cutof)
y[cc==1] <- cutof
y[17] <- abs(y[17])+2*sd(y)
LI <- y
LS <- y; LS[cc==1] <- Inf    # Right-censored

# Estimation
set.seed(74689)
est <- EstStempCens(y, x, cc, time2, coord2, LI, LS, init.phi=2.5, init.rho=0.5, init.tau2=0.8,
          type.Data="balanced", method="nlminb", kappa=0, type.S="exponential",
          IMatrix=TRUE, lower.lim=c(0.01,-0.99,0.01), upper.lim=c(30,0.99,20), M=20,
          perc=0.25, MaxIter=300, pc=0.20)

# Diagnostic
set.seed(12345)
diag <- DiagStempCens(est, type.diag="time", diag.plot = TRUE, ck=1)}
}
\seealso{
\code{\link{EstStempCens}}
}
\author{
Katherine L. Valeriano, Victor H. Lachos and Larissa A. Matos
}
