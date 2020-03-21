



#' Discretise Output
#'
#' Modified from \code{smfsb} package
#'
#' @param out a matrix
#' @param out dt the discretization lattice
#'
#' @export
discretise <- function(out, dt=1){

  events <- nrow(out)
  end <- out[events,"time"]
  start <- out[1,"time"]

  len <- (end-start)%/%dt+1
  x <- matrix(nrow=len,ncol=ncol(out))
  target <- 0
  j <- 1

  for(i in 2:events){
    while(out[i,"time"] >= target){
      x[j,-1] <- out[i,-1]
      x[j,1] <- target
      j <- j + 1
      target <- target + dt
    }
  }

  storage.mode(x) <- "integer"
  return(x)
}
