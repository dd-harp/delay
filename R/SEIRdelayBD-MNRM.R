# -------------------------------------------------------------------------------- #
#
#   Modified Next Reaction Method (Anderson 2007) for SEIR w/delay and demographics
#   E->I following deterministic delay
#   Sean Wu (slwu89@berkeley.edu)
#   January 2020
#
# -------------------------------------------------------------------------------- #


#' SEIR model with fixed delay and demographics via Modified Next Reaction Method (**C** version)
#'
#' Using the Modified Next Reaction Method (MNRM) for delays, simulate a SEIR
#' model where the time spent in the E compartment is a deterministic period of
#' length \code{tau}. There are additional events for births in to S and deaths from
#' all compartments according to rate \code{mu} such that the total population is at
#' its stationary distribution.
#'
#' The MNRM is described in: Anderson, D. F. (2007). *A modified next reaction method for simulating chemical systems with time dependent propensities and delays.* Journal of Chemical Physics, 127(21). \url{https://doi.org/10.1063/1.2799998}
#'
#' @param tmax maximum time of simulation
#' @param S0 initial number of susceptibles
#' @param I0 initial number of infectives
#' @param beta effective contact rate
#' @param nu duration of infectiousness
#' @param tau duration of latent period fixed delay
#' @param mu mortality rate
#' @param verbose print diagnostic information
#' @param maxsize maximum size of output matrix; simulation will return early if this is exceeded
#'
#' @useDynLib delay SEIRdelay_BD_MNRM_Cinternal
#'
#' @examples
#'  \dontrun{
#'    out <- SEIRdelay_BD_MNRM_C(tmax = 1e3,S0 = 1000,I0 = 5,beta = 0.001,nu = 1/5,tau = 10,mu = 1/50,verbose = TRUE)
#'    outd <- discretize(out)
#'    matplot(outd[,-1],type="l",col=c("blue","orange","red","purple"),lwd=1.5,lty=1,ylab="Count",xlab="Days")
#'    legend(x = "topleft",pch=rep(16,4),col=c("blue","orange","red","purple"),legend=c("S","E","I","R"))
#' }
#' @export
SEIRdelay_BD_MNRM_C <- function(tmax,S0,I0,beta,nu,tau,mu,verbose=T,maxsize=1e6){
  out <- .Call(
    SEIRdelay_BD_MNRM_Cinternal,
    as.numeric(tmax),
    as.integer(S0),
    as.integer(I0),
    as.numeric(beta),
    as.numeric(nu),
    as.numeric(tau),
    as.numeric(mu),
    as.integer(as.logical(verbose)),
    as.integer(maxsize)
  )
  do.call(cbind,out)
}
