# Generated by using Rcpp::compileAttributes() -> do not edit by hand
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

crossdist <- function(m1) {
    .Call(`_StempCens_crossdist`, m1)
}

crossdist2 <- function(m1) {
    .Call(`_StempCens_crossdist2`, m1)
}

inversa <- function(M) {
    .Call(`_StempCens_inversa`, M)
}

DevCorMatrix <- function(H, phi, kappa, type) {
    .Call(`_StempCens_DevCorMatrix`, H, phi, kappa, type)
}

ScoreVector <- function(yb, x, beta, sigma2, media, PsiInv, Omega, d1Spat, d1Temp, tauF) {
    .Call(`_StempCens_ScoreVector`, yb, x, beta, sigma2, media, PsiInv, Omega, d1Spat, d1Temp, tauF)
}

HessianMatrix <- function(yb, yyb, x, beta, sigma2, PsiInv, Omega, d1Spat, d1Temp, dSpaTem, d2Spat, d2Temp, tauF) {
    .Call(`_StempCens_HessianMatrix`, yb, yyb, x, beta, sigma2, PsiInv, Omega, d1Spat, d1Temp, dSpaTem, d2Spat, d2Temp, tauF)
}

