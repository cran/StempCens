################################################################################
## Utilities Spatio Temporal ##
################################################################################

################################################################################
## Detecting rows and columns to eliminate in the covarince matrix ##
################################################################################

FilasEliminar <- function(tiempos,vecTiemp,nCoor){

  elim <- NULL # Vector with the index of rows and columns to delete in the covariance matrix
  total <- nrow(tiempos)    # Sample size
  tamanho <- nrow(vecTiemp) # Number of temporal observations without repetition
  nk <- tamanho
  nj <- 1

  for (k in 1:nCoor){
    verificar <- tiempos[nj:nk,]
    indices <- seq(((k-1)*tamanho + 1),(k*tamanho),1)
    cont <- 1

    if (k==nCoor){
      for (j in 1:tamanho){
        if (is.na(verificar[j])){
          verificar[j] <- 0
        } # End if
      } # End for
    } # End if

    for (j in 1:tamanho){
      if (verificar[cont]!=vecTiemp[j]){
        elim <- c(elim,indices[j])
      } else {
        cont <- cont + 1
      } # End else
    } # End for

    nj <- k*tamanho - length(elim) + 1
    nk <- min(nj + tamanho - 1,total)
  }# End for

  return (elim)
}

################################################################################
## Spatial correlation matrix ##
################################################################################
# 5 models: 1. exponential: Exponential model
#           2. gaussian: Gaussianmodel
#           3. matern: Matern model
#           4. pow.exp: Powered Exponential model
#           5. spherical: Spherical model

CorrSpatial <- function(H, phi, kappa, type){
  W <- H       # Distance matrix (spatial component)

  if (type=="exponential") {R <- exp(-(abs(H)/phi))}

  if (type=="gaussian") {R <- exp(-(abs(H)/phi)^2)}

  if (type=="matern"){
    H[H==0] <- 1
    R <- (1/(2^(kappa-1)*gamma(kappa)))*(abs(H)/phi)^(kappa)*besselK(abs(H)/phi,kappa)
  }

  if (type=="pow.exp") {R <- exp(-(abs(H)/phi)^(kappa))}

  if (type=="spherical") {
    R <- (matrix(1,nrow(H),ncol(H)) - 1.5*(abs(H)/phi) + 0.5*(abs(H)/phi)^3)
    Haux <- (abs(H)>phi) + 0
    R[Haux==1] <- 0
  }

  R[W==0] <- 1
  return(R)
}

################################################################################
## Temporal correlation matrix ##
################################################################################

CorrTemporal <- function(disTime, rho){
  W <- disTime # Distance matrix (Temporal component)

  Rt <- rho^abs(disTime) # Following the correlation matrix of an autorregressive model AR(1)
  Rt[W==0] <- 1
  return(Rt)
}

################################################################################
## First and second derivative of temporal correlation matrix ##
################################################################################

DevTempCorMatrix <- function(disTime, rho){

  dRT1 <- abs(disTime)*rho^(abs(disTime) - 1)
  dRT2 <- abs(disTime)*(abs(disTime) - 1)*rho^(abs(disTime) - 2)
  diag(dRT1) <- 0 # First derivative correlation matrix
  diag(dRT2) <- 0 # Second derivative correlation matrix

  return(list(devt1=dRT1, devt2=dRT2))
}

################################################################################
## Covariance Matrix for balanced or unbalanced model ##
################################################################################

CovarianceMatrix <- function(phi, rho, tau2, sigma2, distSpa, disTemp, elim, kappa, type.S){

  cs1 <- CorrSpatial(distSpa,phi,kappa,type.S) # Purely spatial correlation matrix of dimension n
  ct1 <- CorrTemporal(disTemp,rho)             # Purely temporal correlation matrix of dimension T
  In <- diag(nrow(cs1)*nrow(ct1))
  Psi1 <- (cs1%x%ct1) + (tau2/sigma2)*In

  if (length(elim)==0){
    Psi <- Psi1 # Balanced data
  } else {
    Psi <- Psi1[-elim,-elim] # Unbalanced data
  }

  Psi <- (Psi + t(Psi))/2 # Correct rounding problems
  Psi <- round(Psi,9)
  PsiInv <- inversa(Psi)  # Inverse of matrix Psi
  PsiInv <- (PsiInv + t(PsiInv))/2 # Correct rounding problems

  return(list(psi=Psi, psiInv=PsiInv))
}

# Spatio-Temporal covariance matrix for balanced data
CovM <- function(phi, rho, tau2, sigma2, distSpa, disTemp, kappa, type.S){

  cs1 <- CorrSpatial(distSpa,phi,kappa,type.S) # Purely spatial correlation matrix of dimension n
  ct1 <- CorrTemporal(disTemp,rho)             # Purely temporal correlation matrix of dimension T
  Cov <- sigma2*(cs1%x%ct1) + tau2*diag(nrow(cs1)*nrow(ct1)) # Covariance matrix for balanced data

  return(Cov)
}

################################################################################
## Effective range for spatial correlation ##
################################################################################

spherical.eq <- function(x, phi, cor){
  sph.eq <- (1-cor) - 1.5*x/phi + 0.5*(x/phi)^3
  return (sph.eq)
}

matern.eq <- function(x, cor, phi, kappa){
  mat.eq <- 1/(2^(kappa-1)*gamma(kappa))*(x/phi)^(kappa)*besselK(x/phi,kappa) - cor
  return (mat.eq)
}

Effective.range <- function(cor, phi, kappa, Sp.model) {

  if (Sp.model=="exponential"){
    Eff.range <- -phi*log(cor)
  }

  if (Sp.model=="gaussian"){
    Eff.range <- sqrt(-phi^2*log(cor))
  }

  if (Sp.model=="matern"){
    Eff.range <- uniroot(matern.eq,lower=0.0001,upper=(10000*phi),cor=cor,phi=phi,kappa=kappa)$root
  }

  if (Sp.model=="pow.exp"){
    Eff.range <- (-phi^kappa*log(cor))^(1/kappa)
  }

  if (Sp.model=="spherical"){
    Eff.range <- uniroot(spherical.eq,lower=0,upper=phi,phi=phi,cor=cor)$root
  }

  Eff.range = as.numeric(Eff.range)
  return (Eff.range)
}


################################################################################
## Function to estimate phi, rho and tau2 ##
################################################################################

FCiNlminb <- function(rhoG, media, sigma2, yb, yyb, dSpatial, dTemp, elimin, kappa, type.S){
  phi.est <- rhoG[1]
  rho.est <- rhoG[2]
  tau2.est <- rhoG[3]
  vero1 <- numeric()

  mat <- CovarianceMatrix(phi.est,rho.est,tau2.est,sigma2,dSpatial,dTemp,elimin,kappa,type.S)
  Psi <- mat$psi       # Matrix Psi
  PsiInv <- mat$psiInv # Inverse of matrix Psi
  detPsi <- det(Psi)   # Determinant of Psi
  if (detPsi > .Machine$double.xmax){detPsi <- .Machine$double.xmax
  }else{if (detPsi < .Machine$double.xmin){detPsi <- .Machine$double.xmin}} # Avoid Inf and zero

  AA <- (sum(diag(yyb%*%PsiInv)) - t(media)%*%PsiInv%*%yb - t(yb)%*%PsiInv%*%media + t(media)%*%PsiInv%*%media)
  vero1 <- (-0.5)*(log(detPsi) + AA/sigma2) # Log-likelihood to be maximized

  return(-vero1)
}

################################################################################
## Function to estimate phi and rho ##
################################################################################

FCi <- function(rhoG, media, sigma2, tau2.est, yb, yyb, dSpat, dTemp, elimin, kappa, type.S){
  phi.est <- rhoG[1]
  rho.est <- rhoG[2]
  vero1 <- numeric()

  mat <- CovarianceMatrix(phi.est,rho.est,tau2.est,sigma2,dSpat,dTemp,elimin,kappa,type.S)
  Psi <- mat$psi       # Matrix Psi
  PsiInv <- mat$psiInv # Inverse matrix Psi
  detPsi <- det(Psi)   # Determinant of Psi
  if (detPsi > .Machine$double.xmax){detPsi <- .Machine$double.xmax
  }else{if (detPsi < .Machine$double.xmin){detPsi <- .Machine$double.xmin}} # Avoid Inf and zero

  AA <- (sum(diag(yyb%*%PsiInv)) - t(media)%*%PsiInv%*%yb - t(yb)%*%PsiInv%*%media + t(media)%*%PsiInv%*%media)
  vero1 <- (-0.5)*(log(detPsi) + AA/sigma2) # Log-likelihood to be maximized

  return(-vero1)
}

################################################################################
## Function to estimate tau2 ##
################################################################################

FCik <- function(rhoG, media, sigma2, phi.est, rho.est, yb, yyb, dSpa, dTemp, elimin, kappa, type.S){
  tau2.est <- rhoG
  vero1 <- numeric()

  mat <- CovarianceMatrix(phi.est,rho.est,tau2.est,sigma2,dSpa,dTemp,elimin,kappa,type.S)
  Psi <- mat$psi       # Matrix Psi
  PsiInv <- mat$psiInv # Inverse matrix Psi
  detPsi <- det(Psi)   # Determinant of Psi
  if (detPsi > .Machine$double.xmax){detPsi <- .Machine$double.xmax
  }else{if (detPsi < .Machine$double.xmin){detPsi <- .Machine$double.xmin}} # Avoid Inf and zero

  AA <- (sum(diag(yyb%*%PsiInv)) - t(media)%*%PsiInv%*%yb - t(yb)%*%PsiInv%*%media + t(media)%*%PsiInv%*%media)
  vero1 <- (-0.5)*(log(detPsi) + AA/sigma2) # Log-likelihood to be maximized

  return(-vero1)
}

################################################################################
## Likelihood function ##
################################################################################

LogVerosCens<-function(cc, y, media, Psi, LI, LS){

  GB <- GenzBretz(maxpts = 5e4, abseps = 1e-9, releps = 0)
  m1 <- length(y)  # Number of spatio-temporal observations
  ver = numeric()

  if(sum(cc)==0){
    detP <- det(Psi)
    if (detP > .Machine$double.xmax){detP <- .Machine$double.xmax
    }else{if (detP < .Machine$double.xmin){detP <- .Machine$double.xmin}} # Avoid Inf and zero
    ver <- (-0.5)*(m1*log(2*pi) + log(detP) + t(y - media)%*%inversa(Psi)%*%(y - media))
  } # Non-censored observations

  if(sum(cc)>0){
    if(sum(cc)==m1){
      ver1 <- pmvnorm(lower=LI,upper=LS,mean=as.vector(media),sigma=Psi,algorithm=GB)
      if (ver1 > .Machine$double.xmax){ver1 <- .Machine$double.xmax
      }else{if (ver1 < .Machine$double.xmin){ver1 <- .Machine$double.xmin}} # Avoid Inf values or zero
      ver <- log(ver1)
    } # All observations are censored

    if (sum(cc)==(m1-1)){
      InvPsioo <- (1/Psi[cc==0,cc==0])
      muc <- matrix(media[cc==1],ncol=1) + (y[cc==0] - media[cc==0])*InvPsioo*matrix(Psi[cc==1,cc==0],ncol=1)
      Sc <- Psi[cc==1,cc==1] - InvPsioo*matrix(Psi[cc==1,cc==0],ncol=1)%*%Psi[cc==0,cc==1]
      Sc <- (Sc + t(Sc))/2 # Correct rounding problems

      ver1 <- pmvnorm(lower=LI[cc==1],upper=LS[cc==1],mean=as.vector(muc),sigma=Sc,algorithm=GB)
      if (ver1 > .Machine$double.xmax){ver1 <- .Machine$double.xmax
      }else{if (ver1 < .Machine$double.xmin){ver1 <- .Machine$double.xmin}} # Avoid Inf values or zero
      ver <- log(dnorm(x=y[cc==0],mean=as.numeric(media[cc==0]),sd=sqrt(Psi[cc==0,cc==0]))) + log(ver1)
    } # One observation is not censored

    if(sum(cc)<(m1-1)){
      InvPsioo <- inversa(Psi[cc==0,cc==0])
      muc <- matrix(media[cc==1],ncol=1) + Psi[cc==1,cc==0]%*%InvPsioo%*%(y[cc==0] - media[cc==0])
      Sc <- Psi[cc==1,cc==1] - Psi[cc==1,cc==0]%*%InvPsioo%*%Psi[cc==0,cc==1]
      Sc <- (Sc + t(Sc))/2 # Correct rounding problems

      ver1 <- pmvnorm(lower=LI[cc==1],upper=LS[cc==1],mean=as.vector(muc),sigma=Sc,algorithm=GB)
      if (ver1 > .Machine$double.xmax){ver1 <- .Machine$double.xmax
      }else{if (ver1 < .Machine$double.xmin){ver1 <- .Machine$double.xmin}} # Avoid Inf values or zero
      detP <- det(Psi[cc==0,cc==0])
      if (detP > .Machine$double.xmax){detP <- .Machine$double.xmax
      }else{if (detP < .Machine$double.xmin){detP <- .Machine$double.xmin}} # Avoid Inf and zero

      ver <- (-0.5)*(nrow(InvPsioo)*log(2*pi) + log(detP) +
                       t(y[cc==0] - media[cc==0])%*%InvPsioo%*%(y[cc==0] - media[cc==0])) + log(ver1)
    } # Otherwise
  } # At least one censored observation
  obj.out <- list(ver = ver)

  return(obj.out)
}

################################################################################
## Gibbs Sampler ##
################################################################################

amostradordegibbs <- function(M1, M0, cc1, y1, media, Gama, LI, LS){

  nj <- length(y1)  # Number of spatio-temporal observations
  gammai <- as.vector(media)      	  # Mean
  draws <- matrix(NA,nrow=M1,ncol=nj) # Matrix with the random samples
  draws[1,] <- y1
  t1 <- rep(0,nj)

  if(sum(cc1)==0){
    for(i in 2:M1){
      t1 <- y1
      draws[i,] <- t1
    } # End for
  } # Non-censored observations

  if(sum(cc1)==nj){
    for(i in 2:M1){
      t1 <- as.vector(rtmvnorm(1,mean=gammai,sigma=Gama,lower=LI,upper=LS,algorithm="gibbs",thinning=2))
      draws[i,] <- t1
    } # End for
  } # All observations are censored

  if(sum(cc1)==1){
    InvPsioo <- inversa(Gama[cc1==0,cc1==0])
    t1[cc1==0] <- y1[cc1==0]

    muc <- gammai[cc1==1] + Gama[cc1==1,cc1==0]%*%InvPsioo%*%(y1[cc1==0] - gammai[cc1==0])
    muc <- as.numeric(muc)
    Sc <- Gama[cc1==1,cc1==1] - Gama[cc1==1,cc1==0]%*%InvPsioo%*%Gama[cc1==0,cc1==1]
    Sc <- as.numeric(Sc)

    for(i in 2:M1){
      y_r <- rtnorm(1,mean=muc,sd=(sqrt(Sc)),lower=LI[cc1==1],upper=LS[cc1==1])
      t1[cc1==1] <- y_r
      draws[i,] <- t1
    }
  } # End if. One observation is censored

  if(sum(cc1)==(nj-1)){
    InvPsioo <- 1/(Gama[cc1==0,cc1==0])
    t1[cc1==0] <- y1[cc1==0]
    muc <- matrix(gammai[cc1==1],ncol=1) + (InvPsioo*(y1[cc1==0] - gammai[cc1==0]))*matrix(Gama[cc1==1,cc1==0],ncol=1)
    muc <- as.vector(muc)
    Sc <- Gama[cc1==1,cc1==1] - InvPsioo*(matrix(Gama[cc1==1,cc1==0],ncol=1)%*%t(matrix(Gama[cc1==1,cc1==0],ncol=1)))
    Sc <- (Sc + t(Sc))/2 # Correct rounding problems

    for(i in 2:M1){
      y_r <- rtmvnorm(1,mean=muc,sigma=Sc,lower=LI[cc1==1],upper=LS[cc1==1],algorithm="gibbs",thinning=2)
      t1[cc1==1] <- y_r
      draws[i,] <- t1
    } # End for
  } # End if. One observation is non-censored

  if(sum(cc1)>1 & sum(cc1)<(nj-1)){
      InvPsioo <- inversa(Gama[cc1==0, cc1==0])
      t1[cc1==0] <- y1[cc1==0]

      muc <- matrix(gammai[cc1==1],ncol=1) + Gama[cc1==1,cc1==0]%*%InvPsioo%*%matrix(y1[cc1==0] - gammai[cc1==0], ncol=1)
      muc <- as.vector(muc)
      Sc <- Gama[cc1==1,cc1==1] - Gama[cc1==1,cc1==0]%*%InvPsioo%*%Gama[cc1==0,cc1==1]
      Sc <- (Sc + t(Sc))/2 # Correct rounding problems

      for(i in 2:M1){
        y_r <- rtmvnorm(1,mean=muc,sigma=Sc,lower=LI[cc1==1],upper=LS[cc1==1],algorithm="gibbs",thinning=2)
        t1[cc1==1] <- y_r
        draws[i,] <- t1
      } # End for
  } # End if. Otherwise

  # Amostra com burnin (M0)
  amostragibbs <- draws[(M0+1):M1,] # Delete the M0 initials samples
  obj.out <- list(amostragibbs = amostragibbs)

  return(obj.out)
}

################################################################################
## Estimate model parameter values ##
################################################################################

saemST <- function(y,x,cc,tempo,coord,LI,LS,init.phi,init.rho,init.tau2,tau2.fixo=FALSE,type.Data="balanced",
                   method="nlminb",kappa=0,type.S="exponential",IMatrix=TRUE,lower.lim=c(0.01,-0.99,0.01),
                   upper.lim=c(30,0.99,20),M=20,perc=0.25,MaxIter=300,pc=0.2,error=1e-6){

  m <- length(y) # Number of spatio-temporal observations
  marcador <- 30 # Allow you to use an alternative function to maximize when using nlminb

  # Starting values
  na.y <- ifelse(is.na(y)==TRUE,1,0) # Detecting missing values
  beta1 <- solve(t(x[na.y==0,])%*%x[na.y==0,])%*%t(x[na.y==0,])%*%y[na.y==0] # Initial value for beta
  p <- length(beta1)   # Number of beta parameter
  media <- x%*%beta1   # Mean
  sigma2 <- as.numeric(t(y[na.y==0] - media[na.y==0])%*%(y[na.y==0] - media[na.y==0])/(m-sum(na.y))) # Initial value for sigma2
  tau2 <- init.tau2
  phi <- init.phi
  rho <- init.rho
  teta <- c(beta1,sigma2,tau2,phi,rho) # Initial values

  Theta <- NULL # Matrix with the values for the parameters in each iteration

  # Covariance matrix
  vec <- as.matrix(unique(cbind(coord)))              # Matrix of unique coordinates (without repetitions)
  tiempo <- as.matrix(sort(as.matrix(unique(tempo)))) # Vector of unique times index (without repetitions)
  MDist  <- crossdist(vec)     # Distance matrix of spatial coordinates (without coordinates repetitions)
  disTem <- crossdist2(tiempo) # Distance matrix of temporal index (without repetitions)

  if (type.Data=="balanced"){
    eliminar <- NULL
  }else{if (type.Data=="unbalanced"){
    eliminar <- FilasEliminar(tempo,tiempo,nrow(vec)) # Finding the number of the row where exist a missing time index
  }} # End if

  matCor <- CovarianceMatrix(phi,rho,tau2,sigma2,MDist,disTem,eliminar,kappa,type.S)
  V <- matCor$psiInv       # Inverse matrix Psi
  Cov <- sigma2*matCor$psi # Covariance matrix sigma2*Psi
  Cov <- (Cov + t(Cov))/2  # Correct rounding problems

  criterio <- 10
  count <- 0

  ################################ SAEM algorithm ################################

  MG <- round(M/(1 - perc),0) # Number of samples to gernerate
  M0 <- MG - M                # Number of burn samples

  # Sequence of decreasing positive numbers: smoothing parameter
  if(pc==1){
    seqq <- rep(1,MaxIter)
  } else {
    seqq <- c(1/((((pc*MaxIter)+1):MaxIter)-(pc*MaxIter)))
    seqq <- c(rep(1,MaxIter-length(seqq)),seqq)
  }

  SAEM_y  <- matrix(0,m,1)
  SAEM_yy <- matrix(0,m,m)
  if (IMatrix==TRUE){if(tau2.fixo==FALSE){Louis <- matrix(0,(p+4),(p+4))}else{Louis <- matrix(0,(p+3),(p+3))}
  } else {Louis <- NULL}

  while(criterio==0 | criterio > error){

    count <- count + 1

    # Simulation step (S-step): generating from the truncated conditional distribution
    t1 <- y
    gibbs <- amostradordegibbs(MG,M0,cc,t1,media,Cov,LI,LS)

    uyi <- matrix(gibbs$amostragibbs[,1:m],nrow=M,ncol=m)

    if (IMatrix==TRUE){
      ## Approximation of the second term of the observed information matrix by the Louis'method
      cs <- CorrSpatial(MDist,phi,kappa,type.S) # Spatial correlation matrix
      ct <- CorrTemporal(disTem,rho) # Temporal correlation matrix
      Omega  <- (cs%x%ct)              # Spatio-temporal correlation matrix
      d1Spat <- ((DevCorMatrix(MDist,phi,kappa,type.S)$dev1)%x%ct) # First spatial derivative of the spatio-temporal correlation matrix
      d1Temp <- (cs%x%(DevTempCorMatrix(disTem,rho)$devt1))        # First temporal derivative of the spatio-temporal correlation matrix

      if (length(eliminar)>0){
        Omega  <- Omega[-eliminar,-eliminar]
        d1Spat <- d1Spat[-eliminar,-eliminar]
        d1Temp <- d1Temp[-eliminar,-eliminar]
      } # Used if we are dealing with unbalanced data

      if (tau2.fixo==FALSE){soma <- matrix(0,ncol=(p+4),nrow=(p+4))
      }else{soma <- matrix(0,ncol=(p+3),nrow=(p+3))}
      for (i in 1:M){
        yi <- matrix(uyi[i,],nrow=m,ncol=1)
        score <- ScoreVector(yi,x,beta1,sigma2,media,V,Omega,d1Spat,d1Temp,tau2.fixo)
        soma <- soma + score%*%t(score)
      }
      Louis <- Louis + seqq[count]*(soma/M - Louis)
    } # End if. IMatrix == TRUE

    # Approximation step (A-step)
    somay <- matrix(0,m,1)
    somayy <- matrix(0,m,m)

    for(k in 1:M){
      yi <- matrix(uyi[k,],nrow=m,ncol=1)
      somay <- somay + yi
      somayy <- somayy + (yi%*%t(yi))
    }
    E_y <- (1/M)*somay
    E_yy <- (1/M)*somayy

    ## Stochastic approximation
    SAEM_y  <- SAEM_y + seqq[count]*(E_y - SAEM_y)
    SAEM_yy <- SAEM_yy + seqq[count]*(E_yy - SAEM_yy)

    ## Conditional maximization (CM-step)
    beta1 <- solve(t(x)%*%V%*%x)%*%t(x)%*%V%*%SAEM_y  # Update beta
    media <- x%*%beta1           # Update mean
    sigma2 <- (1/m)*(sum(diag(SAEM_yy%*%V)) - t(SAEM_y)%*%V%*%media - t(media)%*%V%*%SAEM_y + t(media)%*%V%*%media)
    sigma2 <- as.numeric(sigma2) # Update sigma2

    if (tau2.fixo==TRUE){

      if (method=="optim"){
        rhos <- optim(c(phi,rho),method="L-BFGS-B",fn=FCi,lower=lower.lim,upper=upper.lim,media=media,sigma2=sigma2,tau2.est=tau2,yb=SAEM_y,
                      yyb=SAEM_yy,dSpat=MDist,dTemp=disTem,elimin=eliminar,kappa=kappa,type.S=type.S)$par
        phi <- rhos[1] # Update phi
        rho <- rhos[2] # Update rho
      } # End if optim

      if (method=="nlminb"){
        rhos <- nlminb(c(phi,rho),objective=FCi,lower=lower.lim,upper=upper.lim,media=media,sigma2=sigma2,tau2.est=tau2,yb=SAEM_y,
                       yyb=SAEM_yy,dSpat=MDist,dTemp=disTem,elimin=eliminar,kappa=kappa,type.S=type.S)$par
        phi <- rhos[1] # Update phi
        rho <- rhos[2] # Update rho
      } # End if nlminb

    } # End if tau2.fixo==TRUE

    if (tau2.fixo==FALSE) {

      if (method=="optim"){
        tau2 <- optimize(f=FCik,lower=lower.lim[3],upper=upper.lim[3],media=media,sigma2=sigma2,phi.est=phi,rho.est=rho,yb=SAEM_y,yyb=SAEM_yy,dSpa=MDist,
                         dTemp=disTem,elimin=eliminar,kappa=kappa,type.S=type.S)$minimum # Update tau2

        rhos <- optim(c(phi,rho),method="L-BFGS-B",fn=FCi,lower=lower.lim[1:2],upper=upper.lim[1:2],media=media,sigma2=sigma2,tau2.est=tau2,
                      yb=SAEM_y,yyb=SAEM_yy,dSpat=MDist,dTemp=disTem,elimin=eliminar,kappa=kappa,type.S=type.S)$par
        phi <- rhos[1] # Update phi
        rho <- rhos[2] # Update rho
      } # End if optim

      if (method=="nlminb"){
        if (marcador == 30){
          rhos <- nlminb(c(phi,rho,tau2),objective=FCiNlminb,lower=lower.lim,upper=upper.lim,media=media,sigma2=sigma2,yb=SAEM_y,
                         yyb=SAEM_yy,dSpatial=MDist,dTemp=disTem,elimin=eliminar,kappa=kappa,type.S=type.S)$par
          phi <- rhos[1]  # Update phi
          rho <- rhos[2]  # Update rho
          tau2 <- rhos[3] # Update tau2

          # Allow you to use the functions nlminb+optimize
          if (abs(phi - (upper.lim[1]+lower.lim[1])/2) > (upper.lim[1]-lower.lim[1]-0.05)/2 | abs(rho - (upper.lim[2]+lower.lim[2])/2) > (upper.lim[2]-lower.lim[2]-0.02)/2 | abs(tau2 - (upper.lim[3]+lower.lim[3])/2) > (upper.lim[3]-lower.lim[3]-0.05)/2){
              beta1 <- solve(t(x[na.y==0,])%*%x[na.y==0,])%*%t(x[na.y==0,])%*%y[na.y==0] # Initial values for beta
              media <- x%*%beta1  # Mean
              sigma2 <- as.numeric(t(y[na.y==0] - media[na.y==0])%*%(y[na.y==0] - media[na.y==0])/(m-sum(na.y))) # Initial valeu for sigma2
              tau2 <- init.tau2  # Initialize tau2
              phi <- init.phi    # Initialize phi
              rho <- init.rho    # Initialize rho
              marcador <- 5       # Let use the combination nlminb+optimize
              count <- 0          # Initialize the total number of iterations
          } # Restarting initial values
        }else{
          tau2 <- optimize(f=FCik,lower=lower.lim[3],upper=upper.lim[3],media=media,sigma2=sigma2,phi.est=phi,rho.est=rho,yb=SAEM_y,yyb=SAEM_yy,dSpa=MDist,
                           dTemp=disTem,elimin=eliminar,kappa=kappa,type.S=type.S)$minimum # Update tau2

          rhos <- nlminb(c(phi,rho),objective=FCi,lower=lower.lim[1:2],upper=upper.lim[1:2],media=media,sigma2=sigma2,tau2.est=tau2,yb=SAEM_y,
                         yyb=SAEM_yy,dSpat=MDist,dTemp=disTem,elimin=eliminar,kappa=kappa,type.S=type.S)$par
          phi <- rhos[1] # Update phi
          rho <- rhos[2] # Update rho
          } # End if marcador==5
      } # End if nlminb
    } # End if tau2.fixo==FALSE

    # Update covariance matrix
    matCor <- CovarianceMatrix(phi,rho,tau2,sigma2,MDist,disTem,eliminar,kappa,type.S)
    V <- matCor$psiInv       # Inverse of matrix Psi
    Cov <- sigma2*matCor$psi # Covariace matrix
    Cov <- (Cov + t(Cov))/2 # Correct rounding problems

    teta1 <- c(beta1,sigma2,tau2,phi,rho)

    if (tau2.fixo==FALSE){
      criterio <- sqrt((teta1/teta-1)%*%(teta1/teta-1))
    }else{criterio <- sqrt((teta1[-(p+2)]/teta[-(p+2)]-1)%*%(teta1[-(p+2)]/teta[-(p+2)]-1))}

    if (count==(MaxIter-1)){criterio <- 10^(-12)}

    Theta <- rbind(Theta,teta1)
    teta <- teta1
  } # End while SAEM

  ############################## End SAEM Algorithm ##############################

 # Criterios AIC, BIC, AICc
  if (tau2.fixo==FALSE){npar <- length(c(teta1))}else{npar <- length(c(teta1))-1}
  loglik <- LogVerosCens(cc,y,media,Cov,LI,LS)$ver
  AICc <- -2*loglik + 2*npar
  AICcorr <- AICc + ((2*npar*(npar + 1))/(m - npar - 1))
  BICc <- -2*loglik + log(m)*npar

  HessianM <- NULL
  LouisM <- NULL

  if (IMatrix==TRUE){
    ## The last approximation of the observed information matrix by the Louis' method
    t1 <- y
    gibbs <- amostradordegibbs(MG,M0,cc,t1,media,Cov,LI,LS)

    uyi <- matrix(gibbs$amostragibbs[,1:m],nrow=M,ncol=m)

    cs <- CorrSpatial(MDist,phi,kappa,type.S)   # Spatial correlation matrix
    ct <- CorrTemporal(disTem,rho)              # Temporal correlation matrix
    Omega  <- (cs%x%ct)                         # Spatio-temporal correlation matrix
    dSpatial <- DevCorMatrix(MDist,phi,kappa,type.S) # Derivatives of the spatial correlation matrix
    dTemporal <- DevTempCorMatrix(disTem,rho)        # Derivatives of the temporal correlation matrix

    d1Spat <- ((dSpatial$dev1)%x%ct)   # First spatial derivative of the spatio-temporal correlation matrix
    d1Temp <- (cs%x%(dTemporal$devt1)) # First temporal derivative of the spatio-temporal correlation matrix
    dSpaTe <- ((dSpatial$dev1)%x%(dTemporal$devt1)) # Second spatio-temporal derivative of the spatio-temporal correlation matrix
    d2Spat <- ((dSpatial$dev2)%x%ct)   # Second spatial derivative of the spatio-temporal correlation matrix
    d2Temp <- (cs%x%(dTemporal$devt2)) # Second temporal derivative of the spatio-temporal correlation matrix

    if (length(eliminar)>0){
      Omega  <- Omega[-eliminar,-eliminar]
      d1Spat <- d1Spat[-eliminar,-eliminar]
      d1Temp <- d1Temp[-eliminar,-eliminar]
      dSpaTe <- dSpaTe[-eliminar,-eliminar]
      d2Spat <- d2Spat[-eliminar,-eliminar]
      d2Temp <- d2Temp[-eliminar,-eliminar]
    } # Using if we are dealing with unbalanced data

    somay <- matrix(0,m,1)
    somayy <- matrix(0,m,m)
    if (tau2.fixo==FALSE){soma <- matrix(0,ncol=(p+4),nrow=(p+4))}else{soma <- matrix(0,ncol=(p+3),nrow=(p+3))}
    for (i in 1:M){
      yi <- matrix(uyi[i,],nrow=m,ncol=1)
      somay <- somay + yi
      somayy <- somayy + (yi%*%t(yi))
      score <- ScoreVector(yi,x,beta1,sigma2,media,V,Omega,d1Spat,d1Temp,tau2.fixo)
      soma <- soma + score%*%t(score)
    }
    Louis <- Louis + seqq[(count + 1)]*(soma/M - Louis)
    SAEM_y  <- SAEM_y + seqq[(count + 1)]*(somay/M - SAEM_y) # Update the approximation of the first moment
    SAEM_yy <- SAEM_yy + seqq[(count + 1)]*(somayy/M - SAEM_yy) # Update the approximation of the second moment

    ### The negative of the conditional expected second derivative matrix given the observed values
    HessianM <- HessianMatrix(SAEM_y,SAEM_yy,x,beta1,sigma2,V,Omega,d1Spat,d1Temp,dSpaTe,d2Spat,d2Temp,tau2.fixo)$HessM

    ### The observed information matrix using the Louis' method
    LouisM <- HessianM - Louis
  } # End IMatrix == TRUE

  # Compute the effective range for spatial correlation model
  range = Effective.range(cor=0.05,phi=phi,kappa=kappa,Sp.model=type.S)

  # Results
  data.model <- list(y=y,cc=cc,x=x,time=tempo,coord=coord,LI=LI,LS=LS,kappa=kappa,type.S=type.S,method=method,lower=lower.lim,upper=upper.lim,
                     initphi=init.phi,initrho=init.rho,initau=init.tau2,tauFixo=tau2.fixo,eliminar=eliminar,M=M,perc=perc,MaxIter=MaxIter,pc=pc,error=error)

  results.model <- list(theta=teta1,Theta=Theta,beta=beta1,sigma2=sigma2,tau2=tau2,phi=phi,rho=rho,Eff.range=range,PsiInv=V,Cov=Cov,SAEMy=SAEM_y,SAEMyy=SAEM_yy,
                       Hessian=HessianM,Louis=LouisM,loglik=loglik,AIC=AICc,BIC=BICc,AICcorr=AICcorr,iteration=count)

  return(list(m.data=data.model,m.results=results.model))

} # End SpatioTemporal function

################################################################################
## Function to predict values ##
################################################################################

PredictNewValues <- function(modelo,loc.Pre,time.Pre,x.pre){

  # Data
    time.Pre <- as.matrix(time.Pre)
    loc.Obs  <- modelo$m.data$coord
    time.Obs <- as.matrix(modelo$m.data$time)
    y.obs <- modelo$m.results$SAEMy      # y (estimated in the case of censored observations)
    x.obs <- modelo$m.data$x
    phi.est <- modelo$m.results$phi     # Phi estimated
    rho.est <- modelo$m.results$rho     # Rho estimated
    tau.est <- modelo$m.results$tau2    # Tau2 estimated
    sig.est <- modelo$m.results$sigma2  # Sigma2 estimated
    bet.est <- modelo$m.results$beta    # Beta estimated
    p <- length(bet.est) # Number of beta parameters

    H1 <- crossdist(loc.Pre)   # Spatial distances for the new locations
    T1 <- crossdist2(time.Pre) # Time distances for the new observations
    CorrH1 <- CorrSpatial(H1, phi.est, modelo$m.data$kappa, modelo$m.data$type.S) # Spatial Correlation
    CorrT1 <-  CorrTemporal(T1, rho.est) # Temporal correlation
    sigma.pp <- sig.est*(CorrH1*CorrT1) + tau.est*diag(ncol(H1)) # Spatio-temporal covariance matrix

    np <- nrow(loc.Pre)
    no <- nrow(loc.Obs)
    H2 <- crossdist(as.matrix(rbind(loc.Pre,loc.Obs)))[1:np,(np+1):(np+no)]
    T2 <- crossdist2(rbind(time.Pre,time.Obs))[1:np,(np+1):(np+no)]
    CorrH2 <- CorrSpatial(H2, phi.est, modelo$m.data$kappa, modelo$m.data$type.S)
    CorrT2 <-  CorrTemporal(T2, rho.est)
    sigma.po <- sig.est*(CorrH2*CorrT2)  # Spatio-temporal covariance between the new space-time indexes and the space-time indexes used in the estimation process

    sigma.oo <- (1/sig.est)*modelo$m.results$PsiInv # Inverse of the spatio-temporal variance matrix for the indexes used in the estimation process

    y.pre <- x.pre%*%bet.est + sigma.po%*%sigma.oo%*%(y.obs-x.obs%*%bet.est)  # Predicted values
    var.pre <- sigma.pp - sigma.po%*%sigma.oo%*%t(sigma.po)                   # Variance of the predicted values

  return(list(predValues=y.pre, VarPred=var.pre))
}

################################################################################
## Mean Squared Prediction Error ##
################################################################################

CrossValidation <- function(yobs,ypred){
  nj <- length(yobs)
  bias <- (ypred - yobs)

  mspe <- sum(bias^2)/nj   # Mean squared prediction error
  rmspe <- sqrt(mspe)      # Root mean squared prediction error
  mae <- sum(abs(bias))/nj # Mean absolute error

  return(list(Bias=bias, Mspe=mspe, Rmspe=rmspe, Mae=mae))
}

################################################################################
## Global Influence Analysis ##
################################################################################

GlobalInf <- function(model, type){

  y <- model$m.data$y
  cc <- model$m.data$cc
  x <- model$m.data$x
  time <- model$m.data$time
  coord <- model$m.data$coord
  p <- length(model$m.results$beta)
  tauFixo <- model$m.data$tauFixo

  if(tauFixo==FALSE){parametro <- model$m.results$theta}else{parametro <- model$m.results$theta[-(p+2)]}

  hess <- model$m.results$Hessian

  if(type=="individual"){
    n <- length(y)
    pb <- txtProgressBar(min = 0, max =n, style = 3)
    GD <- matrix(0,nrow=n,ncol=(length(parametro)+1))
    for (i in 1:n){
      setTxtProgressBar(pb,i)
      newTheta <- saemST(y[-i],x[-i,],cc[-i],matrix(time[-i,],ncol=1),coord[-i,],model$m.data$LI[-i],model$m.data$LS[-i],model$m.data$initphi,model$m.data$initrho,model$m.data$initau,tau2.fixo=tauFixo,type.Data="unbalanced",
                         model$m.data$method,model$m.data$kappa,model$m.data$type.S,IMatrix=FALSE,model$m.data$lower,model$m.data$upper,M=model$m.data$M,perc=model$m.data$perc,MaxIter=model$m.data$MaxIter,
                         pc=model$m.data$pc,error=model$m.data$error)

      if(tauFixo==FALSE){newTheta1 <- newTheta$m.results$theta}else{newTheta1 <- newTheta$m.results$theta[-(p+2)]}

      GlobalDiag <- t(parametro - newTheta1)%*%(hess)%*%(parametro - newTheta1)
      GD[i,] <- c(newTheta1,GlobalDiag)
    } # End for
  } # End if type=="individual"

  if (type=="time"){
    tiempo <- as.matrix(sort(as.matrix(unique(time))))
    n <- length(tiempo)
    pb <- txtProgressBar(min = 0, max =n, style = 3)
    GD <- matrix(0,nrow=n,ncol=(length(parametro)+1))
    for (i in 1:n){
      setTxtProgressBar(pb,i)
      k <- time!=tiempo[i]

      newTheta <- saemST(y[k],x[k,],cc[k],matrix(time[k,],ncol=1),coord[k,],model$m.data$LI[k],model$m.data$LS[k],model$m.data$initphi,model$m.data$initrho,model$m.data$initau,tau2.fixo=tauFixo,type.Data="unbalanced",
                         model$m.data$method,model$m.data$kappa,model$m.data$type.S,IMatrix=FALSE,model$m.data$lower,model$m.data$upper,M=model$m.data$M,perc=model$m.data$perc,MaxIter=model$m.data$MaxIter,
                         pc=model$m.data$pc,error=model$m.data$error)

      if(tauFixo==FALSE){newTheta1 <- newTheta$m.results$theta}else{newTheta1 <- newTheta$m.results$theta[-(p+2)]}

      GlobalDiag <- t(parametro - newTheta1)%*%(hess)%*%(parametro - newTheta1)
      GD[i,] <- c(newTheta1,GlobalDiag)
    } # End for
  } # End if type=="time

  if (type=="location"){
    vec <- as.matrix(unique(cbind(coord)))
    n <- nrow(vec)
    pb <- txtProgressBar(min = 0, max =n, style = 3)
    GD <- matrix(0,nrow=n,ncol=(length(parametro)+1))

    for (i in 1:n){
      setTxtProgressBar(pb,i)
      k <- (coord[,1]!=vec[i,1]|coord[,2]!=vec[i,2])

      newTheta <- saemST(y[k],x[k,],cc[k],matrix(time[k,],ncol=1),coord[k,],model$m.data$LI[k],model$m.data$LS[k],model$m.data$initphi,model$m.data$initrho,model$m.data$initau,tau2.fixo=tauFixo,type.Data="unbalanced",
                         model$m.data$method,model$m.data$kappa,model$m.data$type.S,IMatrix=FALSE,model$m.data$lower,model$m.data$upper,M=model$m.data$M,perc=model$m.data$perc,MaxIter=model$m.data$MaxIter,
                         pc=model$m.data$pc,error=model$m.data$error)

      if(tauFixo==FALSE){newTheta1 <- newTheta$m.results$theta}else{newTheta1 <- newTheta$m.results$theta[-(p+2)]}

      GlobalDiag <- t(parametro - newTheta1)%*%(hess)%*%(parametro - newTheta1)
      GD[i,] <- c(newTheta1,GlobalDiag)
    }
  } # End if type=="location"

  return (GD)
}

################################################################################
## ggplot view ##
################################################################################

vplayout <- function(x, y) viewport(layout.pos.row = x, layout.pos.col = y)


################################################################################
## Simulation of spatio-temporal data (with censored responses) ##
################################################################################

randomStempCens = function(x, time, coords, beta, phi, rho, tau2, sigma, kappa, typeS, typeCens, pcens, lod){
  # Number of spatial and temporal observations
  n = nrow(coords)
  t = length(time)
  coord2 = cbind(rep(coords[,1], each=t), rep(coords[,2], each=t)) # Cartesian coordinates with repetitions
  time2  = as.matrix(rep(time, n))  # Time index with repetitions

  # Mean
  media = x%*%beta
  # Variance-covariance matrix
  Ms  = crossdist(coords) # Spatial distances
  Mt  = crossdist2(time)  # Temporal distances
  Cov = CovM(phi, rho, tau2, sigma, Ms, Mt, kappa, typeS)

  # Simulating dataset
  y    = as.vector(rmvnorm(1, mean=as.vector(media), sigma=Cov))

  if (pcens==0 & is.null(lod)){
    data = data.frame(coord2, time2, y, x)
    names(data) <- c("x.coord", "y.coord", "time", "yObs", paste0("x", 1:ncol(x)))
  }

  # Adding censoring observations
  if (pcens>0){

    if (typeCens=="left"){
      cutoff = quantile(y, pcens)
      cens   = rep(0, length(y)) + (y<cutoff)

    } else {
      if (typeCens=="right"){
        cutoff = quantile(y, 1-pcens)
        cens   = rep(0, length(y)) + (y>cutoff)
      }
    }
    ycens  = y
    ycens[cens==1] = cutoff
    LI = LS = ycens
    if (typeCens=="left"){ LI[cens==1] = -Inf } else { LS[cens==1] = Inf }

    data = data.frame(coord2, time2, cens, ycens, LI, LS, x)
    names(data) <- c("x.coord", "y.coord", "time", "ci", "yObs", "lcl", "ucl", paste0("x", 1:ncol(x)))

  } else {
    if (!is.null(lod)){
      if (typeCens=="left"){
        cens   = rep(0, nrow(data)) + (y<lod)

      } else {
        if (typeCens=="right"){
          cens   = rep(0, nrow(data)) + (y>lod)
        }
      }
      ycens  = y
      ycens[cens==1] = lod
      LI = LS = ycens
      if (typeCens=="left"){ LI[cens==1] = -Inf } else { LS[cens==1] = Inf }

      data = data.frame(coord2, time2, cens, ycens, LI, LS, x)
      names(data) <- c("x.coord", "y.coord", "time", "ci", "yObs", "lcl", "ucl", paste0("x", 1:ncol(x)))
    }
  }
  return(data)
}






