# -------------------------------------------------------------------------------- #
#
#   Modified Next Reaction Method (Anderson 2007) for SEIR w/delay
#   E->I as "CD" reaction
#   Sean Wu (slwu89@berkeley.edu)
#   January 2020
#
# -------------------------------------------------------------------------------- #

# stoiometric matrices
seir_pre <- matrix(
  data=c(
    1,0,1,0, # inf
    0,1,0,0, # latent
    0,0,1,0 # rec
  ),
  nrow=3,ncol=4,byrow=T,
  dimnames=list(c("infection","latent","recovery"),c("S","E","I","R"))
)

seir_post <- matrix(
  data=c(
    0,1,1,0, # inf
    0,0,1,0, # latent
    0,0,0,1 # rec
  ),
  nrow=3,ncol=4,byrow=T,
  dimnames=list(c("infection","latent","recovery"),c("S","E","I","R"))
)

seir_A <- seir_post - seir_pre
seir_S <- t(seir_A)


#' MNRM Algorithm for SEIR (CD)
#'
#' @param beta
#' @param nu
#' @param tau duration of fixed delay for latent period
#'
#' @export
sim_seir_cd <- function(tmax,S0,I0,beta,nu,tau){

  X0 <- c("S"=S0,"E"=0,"I"=I0,R=0)
  X <- X0 <- as.integer(X0)

  if(any(X0<0)){
    stop("initial conditions 'S0' and 'I0' must be non-negative")
  }

  # 1. initialize
  t <- 0
  Pk <- rep(0,3)
  Tk <- rep(0,3)
  delta_t <- rep(0,3)
  ak <- rep(0,3)

  sk <- list(c(Inf),c(Inf),c(Inf)) # completion times of delay

  # 2. draw random numbers
  r <- runif(n=3)

  # 3. set first jump times
  Pk <- log(1/r)

  while(t < tmax){

    # 4. calculate intensities
    ak[1] <- beta*X0[1]*X0[3]
    ak[2] <- beta*X0[1]*X0[3]
    ak[3] <- nu*X0[3]

    # 5.find time required for reach rxn to fire if no others fire first
    for(k in 1:3){
      if(ak[k] != 0){
        delta_t[k] = (Pk[k] - Tk[k]) / ak[k]
      } else {
        delta_t[k] = Inf
      }
    }

    # 6. find index of the minimum reaction
    mu <- which.min(pmin(do.call(rbind,sk)[,1]-t,delta_t))
    delay <- ifelse(sk[[mu]][1]-t < delta_t[mu], TRUE, FALSE)

    # 7. update time
    t <- t + delta_t[mu]

    # 8. update state
    if(delay){
      X <- X + seir_S[mu,]
      sk[[mu]] <- sk[[mu]][-1]
    } else if(mu != 2){
      X <- X + seir_S[mu,]
    } else {
      sk[[mu]] <- c(t+tau,sk[[mu]])
    }

    # 9. update integrated hazards
    Tk <- Tk + ak*delta_t

    # 10. find next jump times

  }

}
