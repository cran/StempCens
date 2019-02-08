
<!-- README.md is generated from README.Rmd. Please edit that file -->
StempCens
=========

The goal of StempCens is to estimates the parameters of a censored or missing data in spatio-temporal models using the SAEM algorithm. This algorithm is a stochastic approximation of the widely used EM algorithm and an important tool for models in which the E-step does not have an analytic form. Besides the expressions obtained to estimate the parameters to the proposed model, we include the calculations for the observed information matrix using the method developed by Thomas (1982). To examine the performance of the fitted model, case-deletion measure are provided. Moreover, it computes the spatio-temporal covariance matrix.

Installation
------------

You can install the released version of StempCens from [CRAN](https://CRAN.R-project.org) with:

``` r
install.packages("StempCens")
```

Example
-------

This is a basic example which shows you how to solve a common problem:

``` r
 # Initial parameter values
 beta <- c(-1,1.50); phi <- 5; rho <- 0.45; tau2 <- 0.80; sigma2 <- 1.5
 # Simulating data
 n1 <- 5    # Number of spatial locations
 n2 <- 5    # Number of temporal index
 set.seed(1000)
 x.coord <- round(runif(n1,0,10),9)   # X coordinate
 y.coord <- round(runif(n1,0,10),9)   # Y coordinate
 coordenadas <- cbind(x.coord,y.coord) # Cartesian coordinates without repetitions
 coord2 <- cbind(rep(x.coord,each=n2),rep(y.coord,each=n2)) # Cartesian coordinates with repetitions
 time <- as.matrix(seq(1,n2,1))      # Time index without repetitions
 time2 <- as.matrix(rep(time,n1))    # Time index with repetitions
 x1 <- rexp(n1*n2,2)
 x2 <- rnorm(n1*n2,2,1)
 x <- cbind(x1,x2)
 media <- x%*%beta
 # Covariance matrix
 H <- as.matrix(dist(coordenadas)) # Spatial distances
 Mt <- as.matrix(dist(time))       # Temporal distances
 Cov <- CovarianceM(phi,rho,tau2,sigma2,distSpa=H,disTemp=Mt,kappa=0,type.S="exponential")
 # Data
 require(mvtnorm)
 y <- as.vector(rmvnorm(1,mean=as.vector(media),sigma=Cov))
 perc <- 0.2
 aa=sort(y);  bb=aa[1:(perc*n1*n2)];  cutof<-bb[perc*n1*n2]
 cc=matrix(1,(n1*n2),1)*(y<=cutof)
 y[cc==1] <- cutof
 # Estimation
 est <- EstStempCens(y, x, cc, time2, coord2, inits.phi=3.5, inits.rho=0.5, inits.tau2=0.7,
                           type.Data="balanced", cens.type="left", method="nlminb", kappa=0,
                           type.S="exponential",
                           IMatrix=TRUE, lower.lim=c(0.01,-0.99,0.01), upper.lim=c(30,0.99,20), M=20,
                           perc=0.25, MaxIter=300, pc=0.2, error = 10^-6)
 
 
```

For the diagnostic analysis in the `EstStempCens` function the parameter input `IMatrix` needs to be `TRUE`.

``` r
diag <- DiagStempCens(est, type.diag="location", diag.plot = TRUE, ck=1)
```
