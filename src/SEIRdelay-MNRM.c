/* --------------------------------------------------------------------------------
#
#   Modified Next Reaction Method (Anderson 2007) for SEIR w/delay
#   Sean Wu (slwu89@berkeley.edu)
#   March 2020
#
-------------------------------------------------------------------------------- */

#include "SEIRdelay-MNRM.h"

// use a heap structure for the sk (queued reactions)
#include "minheap.h"


/* --------------------------------------------------------------------------------
#   SEIR w/delay via Anderson's method
-------------------------------------------------------------------------------- */
// as.integer(as.logical(verbose))
SEXP SEIRdelay_MNRM_C(
  SEXP tmax_r,
  SEXP S0_r,
  SEXP I0_r,
  SEXP beta_r,
  SEXP nu_r,
  SEXP tau_r,
  SEXP verbose_r,
  SEXP maxsize_r
){

  GetRNGstate();

  double tmax = Rf_asReal(tmax_r);
  double beta = Rf_asReal(beta_r);
  double nu = Rf_asReal(nu_r);
  double tau = Rf_asReal(tau_r);

  int verbose = Rf_asInteger(verbose_r);
  int maxsize = Rf_asInteger(maxsize_r);

  int X0[4] = {0};
  int X[4] = {0};

  X0[0] = Rf_asInteger(S0_r);
  X0[2] = Rf_asInteger(I0_r);

  memmove(X,X0,4*sizeof(int));

  for(int i=0; i<4; i++){
    if(X0[i] < 0){
      error("initial conditions 'S0' and 'I0' must be non-negative");
    }
  }

  double t = 0.;

  int out_rows = 1E4;
  SEXP out = PROTECT(Rf_allocMatrix(REALSXP,out_rows,5));
  double* out_ptr = REAL(out);
  int out_i = 0;
  out_ptr[out_i + out_rows*0] = t;
  for(int i=0; i<4; i++){
    out_ptr[out_i + out_rows*(i+1)] = X0[i];
  }
  out_i += 1;

  // 2 reactions: infection (ICD) and recovery (ND)

  // Anderson's method (2007)

  // 1. initialize
  double Pk[2] = {0.};        // next internal firing time of Poisson process (Pk > Tk)
  double Tk[2] = {0.};        // internal time of Poisson process (integrated propensity)
  double delta_t[2] = {0.};   // absolute time to fire
  double ak[2] = {0.};        // propensity functions

  // completion times of delay
  minHeap sk[2];
  sk[0] = initMinHeap(0);
  sk[1] = initMinHeap(0);
  insertNode(&sk[0], DBL_MAX);
  insertNode(&sk[1], DBL_MAX);

  // 2. calculate propensities
  ak[0] = beta*(double)X[0]*(double)X[2];
  ak[1] = nu*(double)X[2];

  // 3. draw internal jump times
  Pk[0] = log(1./runif(0.,1.));
  Pk[1] = log(1./runif(0.,1.));

  double pmin[2];
  int mu;
  int completion;
  double Delta;

  if(verbose){
    Rprintf(" --- beginning simulation --- \n");
  }

  while(t < tmax){

    if(verbose){
      if(out_i % 100 == 0){
        Rprintf(" --- simulated %d reactions at time %f --- \n",out_i,t);
      }
    }

    // 4. set absolute times to fire
    delta_t[0] = (Pk[0] - Tk[0]) / ak[0];
    delta_t[1] = (Pk[1] - Tk[1]) / ak[1];

    // 5. find minimum
    pmin[0] = fmin(sk[0].elem->data - t,delta_t[0]);
    pmin[1] = fmin(sk[1].elem->data - t,delta_t[1]);

    if(pmin[0] > pmin[1]){
      mu = 0;
    } else {
      mu = 1;
    }

    if((sk[mu].elem->data - t) < delta_t[mu]){
      completion = 1;
    } else {
      completion = 0;
    }

    if(completion){
      Delta = sk[mu].elem->data - t;
    } else {
      Delta = delta_t[mu];
    }

    // 6. set t += Delta
    t += Delta;

    // 7 - 10. reaction updating
    if(completion){

      // ICD reaction completes (E->I)

      // update system
      X[1] -= 1;
      X[2] += 1;

      // delete the first row of s_mu
      deleteNode(&sk[mu]);

    } else if(mu == 1){

      // ND reaction fires (I->R)
      X[2] -= 1;
      X[3] += 1;

    } else if(mu == 0){

      // ICD reaction initiates (S->E)

      // update system
      X[0] -= 1;
      X[1] += 1;

      // update s_mu
      insertNode(&sk[mu],t + tau);

    } else {
      error(" --- error: mu %d completion %d --- \n",mu,completion);
    }

    // 11. update Tk
    Tk[0] += (ak[0]*Delta);
    Tk[1] += (ak[1]*Delta);

    // 12. update P_mu
    Pk[mu] += log(1./runif(0.,1.));

    // 13. recalculate propensities
    ak[0] = beta*(double)X[0]*(double)X[2];
    ak[1] = nu*(double)X[2];

    // store output

  }



  PutRNGstate();

  deleteMinHeap(&sk);

  UNPROTECT(1);
  return out;
};
