globalVariables(c("x", "y"))
#' Diagnostic in spatio-temporal model with censored/missing responses
#'
#' Return measures and graphics for diagnostic analysis in spatio-temporal model with censored/missing responses.
#'
#' @param Est.StempCens an object of class \code{Est.StempCens} given as output by the \code{\link{EstStempCens}} function. In the \code{EstStempCens}function, \code{IMatrix} must be \code{TRUE}.
#' @param type.diag type of diagnostic: '\code{individual}' is related when one observation is deleted,
#' '\code{time}' is related when an entire time is deleted, '\code{location}' is related when an entire location is deleted and
#' '\code{all}' the three cases ('\code{individual}', '\code{time}' and '\code{location}').
#' By default \code{type.diag} is \code{individual}.
#' @param diag.plot \code{TRUE} or \code{FALSE}. It indicates if diagnostic plots must be showed. By default = \code{TRUE}.
#' @param ck the value for \code{ck} considered in the benchmark value for the index plot:
#' \eqn{mean(GD)+ck*sd(GD)}, where \eqn{GD} is the vector with all values of the diagnostic measures.
#'
#' @details This function uses the case deletion approach to study the impact of deleting one or
#' more observations from the dataset on the parameters estimates, using the ideas of
#' \insertCite{cook1977detection;textual}{StempCens} and \insertCite{zhu2001case;textual}{StempCens}.
#' The measure is defined by
#'
#' \eqn{GD_i(\theta*)=(\theta* - \theta*[i])'[-Q**(\theta|\theta*)](\theta* - \theta*[i]), i=1,....m,}
#'
#' where \eqn{\theta*} is the estimate of \eqn{\theta} using the complete data, \eqn{\theta*[i]}
#' are the estimates obtained after deletion of the i-th observation (or group of observations) and
#' \eqn{Q**(\theta|\theta*)} is the Hessian matrix.
#'
#' We can eliminate an observation, an entire location or an entire time index.
#'
#' @return The function returns a list with the diagnostic measures.
#'
#' \describe{
#'   \item{If \code{type.diag == individual | time | location}:}{
#'   \code{GD} is a data.frame with the index value of the observation and the GD measure.}
#'   \item{If \code{type.diag == all}:}{
#'   \code{GDind} is a data.frame with the index value of the observation and the GD measure for individual.
#'
#'   \code{GDtime} is a data.frame with the time index value and the GD measure for time.
#'
#'   \code{GDloc} is a data.frame with the side index value and the GD measure for location.
#'
#'   }
#' }
#'
#' @author Katherine L. Valeriano, Victor H. Lachos and Larissa A. Matos
#'
#' @seealso \code{\link{EstStempCens}}
#'
#' @examples
#' \donttest{
#' set.seed(12345)
#' # Parameter values
#' beta <- c(-1,1.5)
#' phi  <- 3
#' rho  <- 0.40
#' tau2 <- 1
#' sigma2 <- 2
#'
#' # Simulating data
#' coord <- matrix(runif(10, 0, 10), ncol=2) # Cartesian coordinates without repetitions
#' time <- as.matrix(1:5) # Time index without repetitions
#' x    <- cbind(rexp(25,2), rnorm(25,2,1))
#' data <- rnStempCens(x, time, coord, beta, phi, rho, tau2, sigma2,
#'                     type.S="exponential", cens="right", pcens=0.20)
#' data$yObs[17] <- abs(data$yObs[17]) + 2*sd(data$yObs) # perturbed observation
#'
#' # Estimation
#' est <- EstStempCens(y=data$yObs, x, cc=data$ci, data$time, cbind(data$x.coord,data$y.coord),
#'                     LI=data$lcl, LS=data$ucl, init.phi=2.5, init.rho=0.5, init.tau2=0.8,
#'                     type.Data="balanced", method="nlminb", kappa=0, type.S="exponential",
#'                     IMatrix=TRUE, lower.lim=c(0.01,-0.99,0.01), upper.lim=c(30,0.99,20), M=20,
#'                     perc=0.25, MaxIter=300, pc=0.20)
#' # Diagnostic
#' set.seed(12345)
#' diag <- DiagStempCens(est, type.diag="time", diag.plot = TRUE, ck=1)}

DiagStempCens = function(Est.StempCens, type.diag="individual", diag.plot=TRUE, ck){

  if(!inherits(Est.StempCens, "Est.StempCens")){stop("An object of the class Est.StempCens must be provided")}

  if(!is.logical(diag.plot)) stop("diag.plot must be TRUE or FALSE")

  if(!type.diag%in%c("all","time","individual","location")){
    stop("type.diag must be all, time, individual or location")}

  if(!is.numeric(ck)) stop("The constant ck must be a real number in [0,Inf)")

  if(ck<0) stop("The constant ck must be a real number in [0,Inf)")

  #---------------------------------------------------------------------#
  #                                Outputs                              #
  #---------------------------------------------------------------------#

  if(type.diag=="individual"|type.diag =="time"|type.diag =="location")
  {
    model = Est.StempCens
    GDm = GlobalInf(model, type=type.diag)
    GD = data.frame(x=seq(1,dim(GDm)[1]),y=GDm[,dim(GDm)[2]])
    line <- mean(GD$y)+ck*sd(GD$y)
    GD1 <- ggplot(aes(x=x,y=y),data=GD) +
      geom_point(colour="gray20",size=1.5) +
      labs(x=paste("Index",type.diag),y=expression(GD[i])) +
      geom_line(y=line,colour="gray40",linetype="dotted") +
      geom_text(aes(label=ifelse(y>line,rownames(GD),'')), size=3, vjust=2)

    #out.diag <- list(GD = GD, plotGD = GD1)
    out.diag <- list(GD = GD)

    if(diag.plot == TRUE){print(GD1)}

  }

  if(type.diag=="all")
  {
    model = Est.StempCens
    cat('\n')
    print("type = individual")
    GDm_1  = GlobalInf(model, type="individual")
    cat('\n')
    print("type = time")
    GDm_2  = GlobalInf(model, type="time")
    cat('\n')
    print("type = location")
    GDm_3  = GlobalInf(model, type="location")

    GD_1 = data.frame(x=seq(1,dim(GDm_1)[1]),y=GDm_1[,dim(GDm_1)[2]])
    line1 <- mean(GD_1$y)+ck*sd(GD_1$y)
    GD1 <- ggplot(aes(x=x,y=y),data=GD_1) +
      geom_point(colour="gray20",size=1.5) +
      labs(x="Index individual",y=expression(GD[i])) +
      geom_line(y=line1,colour="gray40",linetype="dotted") +
      geom_text(aes(label=ifelse(y>line1,rownames(GD_1),'')), size=3, vjust=2)
    GD_2 = data.frame(x=seq(1,dim(GDm_2)[1]),y=GDm_2[,dim(GDm_2)[2]])
    line2 <- mean(GD_2$y)+ck*sd(GD_2$y)
    GD2 <- ggplot(aes(x=x,y=y),data=GD_2) +
      geom_point(colour="gray20",size=1.5) +
      labs(x="Index time",y=expression(GD[i])) +
      geom_line(y=line2,colour="gray40",linetype="dotted") +
      geom_text(aes(label=ifelse(y>line2,rownames(GD_2),'')), size=3, vjust=2)
    GD_3 = data.frame(x=seq(1,dim(GDm_3)[1]),y=GDm_3[,dim(GDm_3)[2]])
    line3 <- mean(GD_3$y)+ck*sd(GD_3$y)
    GD3 <- ggplot(aes(x=x,y=y),data=GD_3) +
      geom_point(colour="gray20",size=1.5) +
      labs(x="Index location",y=expression(GD[i])) +
      geom_line(y=line3,colour="gray40",linetype="dotted") +
      geom_text(aes(label=ifelse(y>line3,rownames(GD_3),'')), size=3, vjust=2)

    out.diag<- list(GDind = GD_1, GDtime = GD_2, GDloc = GD_3)

    if(diag.plot == TRUE){
      pushViewport(viewport(layout = grid.layout(2, 2)))
      print(GD2, vp = vplayout(1, 1))
      print(GD3, vp = vplayout(1, 2))
      print(GD1, vp = vplayout(2, c(1,2)))
    }

  }

  class(out.diag) <- "Diag.StempCens"

  return(out.diag)

}
