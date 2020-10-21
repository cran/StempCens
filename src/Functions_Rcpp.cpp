#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

// Inverse of a symmetric matrix

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::mat inversa(arma::mat M) {
  return (arma::inv(M));
} 


// First and second derivative of spatial correlation matrix

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
List DevCorMatrix(arma::mat H, double phi, double kappa, String type){
  Environment base("package:base");
  Function besselk("besselK");
  
  int n = H.n_rows;
  arma::mat H1(n,n); H1.zeros();
  arma::mat H2(n,n); H2.zeros();
  
  if (type=="exponential"){
    H1 = (abs(H)/pow(phi,2))%exp(-(abs(H)/phi));
    H2 = abs(H)%(abs(H)-2*phi)%exp(-(abs(H)/phi))/pow(phi,4);
  }
  
  if (type=="gaussian"){
    H1 = (2*pow(abs(H),2)/pow(phi,3))%exp(-pow(abs(H)/phi,2));
    H2 = pow(abs(H),2)%(4*pow(abs(H),2) - 6*pow(phi,2))%exp(-pow(abs(H)/phi,2))/pow(phi,6);
  }
  
  if (type=="matern"){
    arma::mat Ak(n,n); arma::mat Bk(n,n);
    H.replace(0,1);
    Ak = as<mat>(besselk(abs(H)/phi,(kappa-1))) + as<mat>(besselk(abs(H)/phi,(kappa+1)));
    Bk = as<mat>(besselk(abs(H)/phi,(kappa-2))) + 2*as<mat>(besselk(abs(H)/phi,kappa)) + as<mat>(besselk(abs(H)/phi,(kappa+2)));
    H1 = (-1/(pow(2,kappa)*pow(phi,2)*tgamma(kappa)))*pow(abs(H)/phi,kappa)%(2*kappa*phi*as<mat>(besselk(abs(H)/phi,kappa)) - abs(H)%Ak);
    H2 = pow(abs(H),kappa)/(pow(2,(kappa+1))*tgamma(kappa)*pow(phi,(kappa+4)))%(4*kappa*(kappa+1)*pow(phi,2)*as<mat>(besselk(abs(H)/phi,kappa)) - 4*(kappa+1)*phi*abs(H)%Ak + pow(abs(H),2)%Bk);
  }
  
  if (type=="pow.exp"){
    H1 = (kappa/phi)*pow(abs(H)/phi,(kappa))%exp(-pow(abs(H)/phi,(kappa)));
    H2 = H1%(kappa*pow(abs(H),kappa)/pow(phi,(kappa+1)) - (kappa+1)/phi);
  }
  
  if (type=="spherical"){
    arma::mat Haux = H;
    H = H - Haux.clean(phi);
    H1 = 1.5*(abs(H)/pow(phi,2)) - 1.5*(pow(abs(H),3)/pow(phi,4));
    H2 = 6*pow(abs(H),3)/pow(phi,5) - 3*abs(H)/pow(phi,3);
  }
  H1 = H1 - diagmat(H1);   // First derivative correlation matrix
  H2 = H2 - diagmat(H2);   // Second derivative correlation matrix
  
  List derSpatial;
  derSpatial["dev1"] = H1;
  derSpatial["dev2"] = H2;
  return(derSpatial);
}


// Score vector for the complete data

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::mat ScoreVector(arma::vec yb, arma::mat x, arma::vec beta, double sigma2, arma::vec media, arma::mat PsiInv, arma::mat Omega, arma::mat d1Spat, arma::mat d1Temp, bool tauF){
  arma::mat score;
  
  // First derivatives
  if(tauF==FALSE){
    arma::vec dbeta; arma::vec dsigma2; arma::vec dtau2; arma::vec dphi; arma::vec drho;
    arma::vec diff = yb - media;
    dbeta   = (1/sigma2)*trans(x)*PsiInv*(diff);
    dsigma2 = (-0.5/sigma2)*trace(PsiInv*Omega) + (0.5/pow(sigma2,2))*as_scalar(trans(diff)*PsiInv*Omega*PsiInv*(diff));
    dtau2   = (-0.5/sigma2)*trace(PsiInv) + (0.5/pow(sigma2,2))*as_scalar(trans(diff)*PsiInv*PsiInv*(diff));
    dphi    = (-0.5)*trace(PsiInv*d1Spat) + (0.5/sigma2)*as_scalar(trans(diff)*PsiInv*d1Spat*PsiInv*(diff));
    drho    = (-0.5)*trace(PsiInv*d1Temp) + (0.5/sigma2)*as_scalar(trans(diff)*PsiInv*d1Temp*PsiInv*(diff));
    
    score = join_vert(dbeta,join_vert(dsigma2,dtau2,dphi,drho));
  }else{
    arma::vec dbeta; arma::vec dsigma2; arma::vec dphi; arma::vec drho;
    arma::vec diff = yb - media;
    dbeta   = (1/sigma2)*trans(x)*PsiInv*(diff);
    dsigma2 = (-0.5/sigma2)*trace(PsiInv*Omega) + (0.5/pow(sigma2,2))*as_scalar(trans(diff)*PsiInv*Omega*PsiInv*(diff));
    dphi    = (-0.5)*trace(PsiInv*d1Spat) + (0.5/sigma2)*as_scalar(trans(diff)*PsiInv*d1Spat*PsiInv*(diff));
    drho    = (-0.5)*trace(PsiInv*d1Temp) + (0.5/sigma2)*as_scalar(trans(diff)*PsiInv*d1Temp*PsiInv*(diff));
    
    score = join_vert(dbeta,dsigma2,dphi,drho);
  }
  return (score);
}


// The negative of the second derivative of Q (-Hessian Matrix)

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
List HessianMatrix(arma::vec yb, arma::mat yyb, arma::mat x, arma::vec beta, double sigma2, arma::mat PsiInv, arma::mat Omega, arma::mat d1Spat, arma::mat d1Temp, arma::mat dSpaTem, arma::mat d2Spat, arma::mat d2Temp, bool tauF){
  
  arma::vec media1 = x*beta; // Mean
  arma::mat HessMatrix;
  
  arma::mat E0 = PsiInv*Omega*PsiInv;
  arma::mat E1 = PsiInv*PsiInv;
  arma::mat E2 = PsiInv*d1Spat*PsiInv;
  arma::mat E3 = PsiInv*d1Temp*PsiInv;
  
  // Second derivatives
  if (tauF==FALSE){
    arma::vec d2sigma; arma::vec DsigmaDtau2; arma::vec DsigmaDphi; arma::vec DsigmaDrho; arma::vec d2tau2;
    arma::vec Dtau2Dphi; arma::vec Dtau2Drho; arma::vec d2phi; arma::vec DphiDrho; arma::vec d2rho;
    arma::vec diff = yb - media1;
    arma::vec diff2 = 2*yb - media1;
    
    arma::mat d2beta      = (trans(x)*PsiInv*x)/sigma2;
    arma::mat DbetaDsigma = (trans(x)*E0*diff)/pow(sigma2,2);
    arma::mat DbetaDtau2  = (trans(x)*E1*diff)/pow(sigma2,2);
    arma::mat DbetaDphi   = (trans(x)*E2*diff)/sigma2;
    arma::mat DbetaDrho   = (trans(x)*E3*diff)/sigma2;
    
    d2sigma     = -0.5*trace(E0*Omega)/pow(sigma2,2) + (1/pow(sigma2,3))*(trace(yyb*E0*Omega*PsiInv) - as_scalar(trans(media1)*E0*Omega*PsiInv*diff2));
    DsigmaDtau2 = -0.5*trace(E0)/pow(sigma2,2) + (0.5/pow(sigma2,3))*(trace(2*yyb*E0*PsiInv) - as_scalar(trans(media1)*(PsiInv*E0 + E0*PsiInv)*diff2));
    DsigmaDphi  = 0.5*trace((PsiInv - E0)*d1Spat)/sigma2 + (0.5/pow(sigma2,2))*(trace(yyb*(2*PsiInv*Omega*E2 - E2)) - as_scalar(trans(media1)*(E2*Omega*PsiInv + PsiInv*Omega*E2 - E2)*diff2));
    DsigmaDrho  = 0.5*trace((PsiInv - E0)*d1Temp)/sigma2 + (0.5/pow(sigma2,2))*(trace(yyb*(2*PsiInv*Omega*E3 - E3)) - as_scalar(trans(media1)*(E3*Omega*PsiInv + PsiInv*Omega*E3 - E3)*diff2));
    
    d2tau2      = -0.5*trace(E1)/pow(sigma2,2) + (1/pow(sigma2,3))*(trace(yyb*PsiInv*E1) - as_scalar(trans(media1)*PsiInv*E1*diff2));
    Dtau2Dphi   = -0.5*trace(E2)/sigma2 + (0.5/pow(sigma2,2))*(trace(2*yyb*E2*PsiInv) - as_scalar(trans(media1)*(E2*PsiInv + PsiInv*E2)*diff2));
    Dtau2Drho   = -0.5*trace(E3)/sigma2 + (0.5/pow(sigma2,2))*(trace(2*yyb*E3*PsiInv) - as_scalar(trans(media1)*(E3*PsiInv + PsiInv*E3)*diff2));
    
    d2phi       = 0.5*trace(PsiInv*d2Spat - E2*d1Spat) + (0.5/sigma2)*(trace(yyb*(2*E2*d1Spat*PsiInv - PsiInv*d2Spat*PsiInv)) - as_scalar(trans(media1)*(2*E2*d1Spat*PsiInv - PsiInv*d2Spat*PsiInv)*diff2));
    DphiDrho    = 0.5*trace(PsiInv*dSpaTem - E3*d1Spat) + (0.5/sigma2)*(trace(yyb*(2*E3*d1Spat*PsiInv - PsiInv*dSpaTem*PsiInv)) - as_scalar(trans(media1)*(E3*d1Spat*PsiInv + PsiInv*d1Spat*E3 - PsiInv*dSpaTem*PsiInv)*diff2));
    
    d2rho       = 0.5*trace(PsiInv*d2Temp - E3*d1Temp) + (0.5/sigma2)*(trace(yyb*(2*E3*d1Temp*PsiInv - PsiInv*d2Temp*PsiInv)) - as_scalar(trans(media1)*(2*E3*d1Temp*PsiInv - PsiInv*d2Temp*PsiInv)*diff2));
    
    // Second derivative matrix
    HessMatrix = join_vert(join_horiz(join_horiz(d2beta,DbetaDsigma,DbetaDtau2,DbetaDphi),DbetaDrho), 
                           join_vert(join_horiz(join_horiz(trans(DbetaDsigma),d2sigma,DsigmaDtau2,DsigmaDphi),DsigmaDrho),
                                     join_horiz(join_horiz(trans(DbetaDtau2),DsigmaDtau2,d2tau2,Dtau2Dphi),Dtau2Drho),
                                     join_horiz(join_horiz(trans(DbetaDphi),DsigmaDphi,Dtau2Dphi,d2phi),DphiDrho),
                                     join_horiz(join_horiz(trans(DbetaDrho),DsigmaDrho,Dtau2Drho,DphiDrho),d2rho)));
  }else{
    arma::vec d2sigma; arma::vec DsigmaDphi; arma::vec DsigmaDrho; arma::vec d2phi; arma::vec DphiDrho; arma::vec d2rho;
    arma::vec diff = yb - media1;
    arma::vec diff2 = 2*yb - media1;
    
    arma::mat d2beta      = (trans(x)*PsiInv*x)/sigma2;
    arma::mat DbetaDsigma = (trans(x)*E0*diff)/pow(sigma2,2);
    arma::mat DbetaDphi   = (trans(x)*E2*diff)/sigma2;
    arma::mat DbetaDrho   = (trans(x)*E3*diff)/sigma2;
    
    d2sigma     = -0.5*trace(E0*Omega)/pow(sigma2,2) + (1/pow(sigma2,3))*(trace(yyb*E0*Omega*PsiInv) - as_scalar(trans(media1)*E0*Omega*PsiInv*diff2));
    DsigmaDphi  = 0.5*trace((PsiInv - E0)*d1Spat)/sigma2 + (0.5/pow(sigma2,2))*(trace(yyb*(2*PsiInv*Omega*E2 - E2)) - as_scalar(trans(media1)*(E2*Omega*PsiInv + PsiInv*Omega*E2 - E2)*diff2));
    DsigmaDrho  = 0.5*trace((PsiInv - E0)*d1Temp)/sigma2 + (0.5/pow(sigma2,2))*(trace(yyb*(2*PsiInv*Omega*E3 - E3)) - as_scalar(trans(media1)*(E3*Omega*PsiInv + PsiInv*Omega*E3 - E3)*diff2));
    
    d2phi       = 0.5*trace(PsiInv*d2Spat - E2*d1Spat) + (0.5/sigma2)*(trace(yyb*(2*E2*d1Spat*PsiInv - PsiInv*d2Spat*PsiInv)) - as_scalar(trans(media1)*(2*E2*d1Spat*PsiInv - PsiInv*d2Spat*PsiInv)*diff2));
    DphiDrho    = 0.5*trace(PsiInv*dSpaTem - E3*d1Spat) + (0.5/sigma2)*(trace(yyb*(2*E3*d1Spat*PsiInv - PsiInv*dSpaTem*PsiInv)) - as_scalar(trans(media1)*(E3*d1Spat*PsiInv + PsiInv*d1Spat*E3 - PsiInv*dSpaTem*PsiInv)*diff2));
    
    d2rho       = 0.5*trace(PsiInv*d2Temp - E3*d1Temp) + (0.5/sigma2)*(trace(yyb*(2*E3*d1Temp*PsiInv - PsiInv*d2Temp*PsiInv)) - as_scalar(trans(media1)*(2*E3*d1Temp*PsiInv - PsiInv*d2Temp*PsiInv)*diff2));
    
    // Second derivative matrix
    HessMatrix = join_vert(join_horiz(join_horiz(d2beta,DbetaDsigma,DbetaDphi),DbetaDrho), 
                           join_vert(join_horiz(join_horiz(trans(DbetaDsigma),d2sigma,DsigmaDphi),DsigmaDrho),
                                     join_horiz(join_horiz(trans(DbetaDphi),DsigmaDphi,d2phi),DphiDrho),
                                     join_horiz(join_horiz(trans(DbetaDrho),DsigmaDrho,DphiDrho),d2rho)));
  }
  List Hess;
  Hess["HessM"] = HessMatrix;
  return(Hess); // Return the negative of the second derivative matrix of Q function
}
