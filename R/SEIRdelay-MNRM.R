# -------------------------------------------------------------------------------- #
#
#   Modified Next Reaction Method (Anderson 2007) for SEIR w/delay
#   E->I following deterministic delay
#   Sean Wu (slwu89@berkeley.edu)
#   January 2020
#
# -------------------------------------------------------------------------------- #

# stoiometric matrices
seir_pre <- matrix(
  data=as.integer(c(
    1,0,1,0, # inf
    0,1,0,0, # latent
    0,0,1,0 # rec
  )),
  nrow=3,ncol=4,byrow=T,
  dimnames=list(c("infection","latent","recovery"),c("S","E","I","R"))
)

seir_post <- matrix(
  data=as.integer(c(
    0,1,1,0, # inf
    0,0,1,0, # latent
    0,0,0,1 # rec
  )),
  nrow=3,ncol=4,byrow=T,
  dimnames=list(c("infection","latent","recovery"),c("S","E","I","R"))
)

seir_A <- seir_post - seir_pre
seir_S <- t(seir_A)

fequal <- function(x,y,tol=sqrt(.Machine$double.eps)){
  abs(x-y) <= tol
}

#' SEIR model with fixed delay via Modified Next Reaction Method
#'
#'
#'
#' @param tmax
#' @param S0
#' @param I0
#' @param beta
#' @param nu
#' @param tau duration of latent period fixed delay
#' @param verbose print diagnostic information
#'
#' @export
SEIRdelay_MNRM <- function(tmax,S0,I0,beta,nu,tau,verbose=T,maxsize=1e6){

  X0 <- c("S"=S0,"E"=0,"I"=I0,"R"=0)
  storage.mode(x = X0) <- "integer"
  X <- X0

  if(any(X0<0)){
    stop("initial conditions 'S0' and 'I0' must be non-negative")
  }

  t <- 0

  out <- matrix(NaN,nrow=1e4,ncol=length(X)+1,dimnames=list(NULL,c("time",names(X))))
  out[1,1] <- t
  out[1,2:ncol(out)] <- X0
  i <- 2

  # 2 reactions: infection and recovery

  # 1. initialize
  Pk <- rep(0,2)        # next internal firing time of Poisson process (Pk > Tk)
  Tk <- rep(0,2)        # internal time of Poisson process (integrated propensity)
  delta_t <- rep(0,2)   # absolute time to fire
  ak <- rep(0,2)        # propensity functions

  sk <- list(c(Inf),c(Inf)) # completion times of delay

  # 2. calculate propensities
  ak[1] <- beta*X0["S"]*X0["I"]
  ak[2] <- nu*X0["I"]

  # 3. sdraw internal jump times
  Pk[] <- log(1/runif(n=2))

  if(verbose){cat(" --- beginning simulation --- \n")}

  while(t < tmax){

    if(verbose){
      if((i %% 100) == 0){
        cat(" --- simulated ",i-1," reactions at time ",t," --- \n")
      }
    }

    # 4. set absolute times to fire
    delta_t <- (Pk - Tk)/ak

    # 5.find minimum
    mu <- which.min(pmin(do.call(rbind,sk)[,1]-t,delta_t))
    completion <- ifelse(sk[[mu]][1]-t < delta_t[mu], TRUE, FALSE)
    Delta <- ifelse(completion,sk[[mu]][1]-t,delta_t[mu])

    # 6. set t += Delta
    t <- t + Delta

    # 7 - 10. reaction updating
    if(completion){

      # ICD reaction completes (E->I)

      stopifnot(mu==1)

      # update system
      X <- X + seir_S[,2]

      # delete the first row of s_mu
      sk[[mu]] <- sk[[mu]][2:length(sk[[mu]])]

    } else if(mu==2){

      # ND reaction fires (I->R)
      X <- X + seir_S[,3]

    } else if(mu==1){

      # ICD reaction initiates (S->E)

      # update system
      X <- X + seir_S[,1]

      # update 2nd to last place of s_mu
      sk[[mu]] <- append(x = sk[[mu]],values = t + tau,after = length(sk[[mu]])-1)

    }

    # 11. update Tk
    Tk[] <- Tk + (ak*Delta)

    # 12. update P_mu
    Pk[mu] <- Pk[mu] + log(1/runif(1))

    # 13. recalculate propensities
    ak[1] <- beta*X["S"]*X["I"]
    ak[2] <- nu*X["I"]

    # STORE OUTPUT
    out[i,1] <- t
    out[i,2:ncol(out)] <- X
    i <- i + 1
    if(i > maxsize){
      cat(" --- warning: exceeded maximum output size, returning output early --- \n")
      return(out[-which(is.nan(out[,1])),])
    }
    if(all(fequal(ak,0)) & all(is.infinite(unlist(sk)))){
      cat(" --- warning: all propentities approximately zero and no delayed reactions queued up, returning output early --- \n")
      return(out[-which(is.nan(out[,1])),])
    }
    if(i > nrow(out)){
      if(verbose){
        cat(" --- extending output matrix --- \n")
      }
      out_orig <- out
      new_size <- nrow(out_orig)+1e3
      if(new_size > maxsize){
        new_size <- maxsize - nrow(out_orig)
      }
      out <- matrix(NaN,nrow=new_size,ncol=length(X)+1,dimnames=list(NULL,names(X)))
      out[1:nrow(out_orig),] <- out_orig
    }

  }

  return(out[-which(is.nan(out[,1])),])
}
