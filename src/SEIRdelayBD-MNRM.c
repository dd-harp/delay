/* --------------------------------------------------------------------------------
#
#   Modified Next Reaction Method (Anderson 2007) for SEIR w/delay and demographics
#   Sean Wu (slwu89@berkeley.edu)
#   March 2020
#
-------------------------------------------------------------------------------- */

#include "SEIRdelayBD-MNRM.h"

// use a heap structure for the sk (queued reactions)
#include "minheap.h"

// comparison of floats
#include "fcmp.h"

// out <- SEIRdelay_BD_MNRM_C(tmax = 1e3,S0 = 1000,I0 = 5,beta = 0.001,nu = 1/5,tau = 10,mu = 1/50,verbose = TRUE)

/* --------------------------------------------------------------------------------
#   SEIR w/delay and demographics via Anderson's method
-------------------------------------------------------------------------------- */

SEXP SEIRdelay_BD_MNRM_Cinternal(
  SEXP tmax_r,
  SEXP S0_r,
  SEXP I0_r,
  SEXP beta_r,
  SEXP nu_r,
  SEXP tau_r,
  SEXP mu_r,
  SEXP verbose_r,
  SEXP maxsize_r
){

  GetRNGstate();

  double tmax = Rf_asReal(tmax_r);
  double beta = Rf_asReal(beta_r);
  double nu = Rf_asReal(nu_r);
  double tau = Rf_asReal(tau_r);
  double mu = Rf_asReal(mu_r);

  int verbose = Rf_asInteger(verbose_r);
  int maxsize = Rf_asInteger(maxsize_r);

  int X0[4] = {0}; // X[0]=S, X[1]=E, X[2]=I, X[4]=R
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

  int out_size = 8E4;
  int out_i = 0;

  double* t_out = (double*)Calloc(out_size,double);
  int* S_out = (int*)Calloc(out_size,int);
  int* E_out = (int*)Calloc(out_size,int);
  int* I_out = (int*)Calloc(out_size,int);
  int* R_out = (int*)Calloc(out_size,int);

  t_out[out_i] = t;
  S_out[out_i] = X0[0];
  E_out[out_i] = X0[1];
  I_out[out_i] = X0[2];
  R_out[out_i] = X0[3];
  out_i += 1;

  // 7 reactions
  // S + I -> E + I ; E -> I (infection followed by completion of infectivity)
  // S + I -> E + I ; E -> 0 (infection followed by death)
  // I -> R                  (recovery)
  // 0 -> S                  (birth)
  // S -> 0                  (death of S)
  // I -> 0                  (death of I)
  // R -> 0                  (death of R)

  // Anderson's method (2007)

  // 1. initialize
  double Pk[7] = {0.};        // next internal firing time of Poisson process (Pk > Tk)
  double Tk[7] = {0.};        // internal time of Poisson process (integrated propensity)
  double delta_t[7] = {0.};   // absolute time to fire
  double ak[7] = {0.};        // propensity functions

  // completion times of delay
  minHeap sk[2];
  sk[0] = initMinHeap(0);
  sk[1] = initMinHeap(0);
  insertNode(&sk[0], DBL_MAX);
  insertNode(&sk[1], DBL_MAX);

  // 2. calculate propensities
  double lambda = beta*(double)X[0]*(double)X[2];
  double p_surv = exp(-mu*tau);
  ak[0] = lambda * p_surv;            // (infection followed by completion of infectivity)
  ak[1] = lambda * (1. - p_surv);     // (infection followed by death)
  ak[2] = nu*(double)X[2];            // (recovery)
  ak[3] = mu * (double)(X[0]+X[1]+X[2]+X[3]); // (birth)
  ak[4] = mu * (double)X[0];                  // (death of S)
  ak[5] = mu * (double)X[2];                  // (death of I)
  ak[6] = mu * (double)X[3];                  // (death of R)

  // Rprintf(" --- initial propensities ak[0] %f, ak[1] %f, ak[2] %f, ak[3] %f, ak[4] %f,ak[5] %f, ak[6] %f --- \n",ak[0],ak[1],ak[2],ak[3],ak[4],ak[5],ak[6]);

  // 3. draw internal jump times
  for(int i=0; i<7; i++){
    Pk[i] = log(1./runif(0.,1.));
  }

  double pmin[7];
  int j; // mu in SEIRdelay-MNRM code
  int completion;
  double Delta;

  if(verbose){
    Rprintf(" --- beginning simulation --- \n");
  }

  while(t < tmax){

    // Rprintf(" --- simulation iteration %d, time %f --- \n",out_i,t);

    if(verbose){
      if(out_i % 100 == 0){
        Rprintf(" --- simulated %d reactions at time %f --- \n",out_i,t);
      }
    }

    // 4. set absolute times to fire
    for(int i=0; i<7; i++){
      delta_t[i] = (Pk[i] - Tk[i]) / ak[i];
    }

    // 5. find minimum and time of next event
    for(int i=0; i<7; i++){
      if(i < 2){
        pmin[i] = fmin(sk[i].elem->data - t, delta_t[i]);
      } else {
        pmin[i] = delta_t[i];
      }
    }

    j=0;
    for(int i=1; i<7; i++){
      if(pmin[j] > pmin[i]){
        j = i;
      }
    }

    completion = 0;
    if(j < 2){
      if((sk[j].elem->data - t) < delta_t[j]){
        completion = 1;
      }
    }

    if(completion){
      Delta = sk[j].elem->data - t;
    } else {
      Delta = delta_t[j];
    }

    // 6. set t+= Delta
    t += Delta;

    // Rprintf(" --- reaction %d fired, completion? %d --- \n",j,completion);

    // 7 - 10. reaction updating
    if(completion){
      // E->I
      if(j==0){

        // update system
        X[1] -= 1;
        X[2] += 1;

        // delete first row of s_mu
        deleteNode(&sk[j]);

      // E->0
      } else if(j==1){

        // update system
        X[1] -= 1;

        // delete first row of s_mu
        deleteNode(&sk[j]);

      } else {
        error(" --- warning, impossible combination of reaction (%d) and completion (%d), aborting simulation --- \n",j,completion);
      }
    } else {

      // S + I -> E + I (will complete)
      if(j==0){

        X[0] -= 1;
        X[1] += 1;

        insertNode(&sk[j],t + tau);

      // S + I -> E + I (will die)
      } else if(j==1){

        X[0] -= 1;
        X[1] += 1;

        // sample truncted exponential by inversion
        double u = runif(0.,pexp(tau,1./mu,1,0));
        double s = qexp(u,1./mu,1,0);

        insertNode(&sk[j],t + s);

      // (recovery)
      } else if(j==2){

        X[2] -= 1;
        X[3] += 1;

      // birth
      } else if(j==3){

        X[0] += 1;

      // S->0
      } else if(j==4){

        X[0] -= 1;

      // I->0
      } else if(j==5){

        X[2] -= 1;

      // R->0
      } else if(j==6){

        X[3] -= 1;

      } else {
        error(" --- warning, impossible combination of reaction (%d) and completion (%d), aborting simulation --- \n",j,completion);
      }

    }

    // 11. update Tk
    for(int i=0; i<7; i++){
      Tk[i] += (ak[i]*Delta);
    }

    // 12. update P_mu
    Pk[j] += log(1./runif(0.,1.));

    // 13. recalculate propensities
    double lambda = beta*(double)X[0]*(double)X[2];
    double p_surv = exp(-mu*tau);
    ak[0] = lambda * p_surv;            // (infection followed by completion of infectivity)
    ak[1] = lambda * (1. - p_surv);     // (infection followed by death)
    ak[2] = nu * (double)X[2];            // (recovery)
    ak[3] = mu * (double)(X[0]+X[1]+X[2]+X[3]); // (birth)
    ak[4] = mu * (double)X[0];                  // (death of S)
    ak[5] = mu * (double)X[2];                  // (death of I)
    ak[6] = mu * (double)X[3];                  // (death of R)

    // Rprintf(" --- recalculate propensities ak[0] %f, ak[1] %f, ak[2] %f, ak[3] %f, ak[4] %f,ak[5] %f, ak[6] %f --- \n",ak[0],ak[1],ak[2],ak[3],ak[4],ak[5],ak[6]);


    // store output
    t_out[out_i] = t;
    S_out[out_i] = X[0];
    E_out[out_i] = X[1];
    I_out[out_i] = X[2];
    R_out[out_i] = X[3];
    out_i += 1;

    // various conditions to return early
    if(out_i > maxsize){
      // Rprintf(" --- warning: exceeded maximum output size, returning output early --- \n");

      SEXP out = PROTECT(Rf_allocVector(VECSXP,5));

      SEXP t_2R = PROTECT(Rf_allocVector(REALSXP,out_i-1));
      memmove(REAL(t_2R),t_out,(out_i-1) * sizeof(double));

      SEXP S_2R = PROTECT(Rf_allocVector(INTSXP,out_i-1));
      memmove(INTEGER(S_2R),S_out,(out_i-1) * sizeof(int));

      SEXP E_2R = PROTECT(Rf_allocVector(INTSXP,out_i-1));
      memmove(INTEGER(E_2R),E_out,(out_i-1) * sizeof(int));

      SEXP I_2R = PROTECT(Rf_allocVector(INTSXP,out_i-1));
      memmove(INTEGER(I_2R),I_out,(out_i-1) * sizeof(int));

      SEXP R_2R = PROTECT(Rf_allocVector(INTSXP,out_i-1));
      memmove(INTEGER(R_2R),R_out,(out_i-1) * sizeof(int));

      SET_VECTOR_ELT(out,0,t_2R);
      SET_VECTOR_ELT(out,1,S_2R);
      SET_VECTOR_ELT(out,2,E_2R);
      SET_VECTOR_ELT(out,3,I_2R);
      SET_VECTOR_ELT(out,4,R_2R);

      SEXP names = PROTECT(Rf_allocVector(STRSXP,5));
      SET_STRING_ELT(names,0,Rf_mkChar("time"));
      SET_STRING_ELT(names,1,Rf_mkChar("S"));
      SET_STRING_ELT(names,2,Rf_mkChar("E"));
      SET_STRING_ELT(names,3,Rf_mkChar("I"));
      SET_STRING_ELT(names,4,Rf_mkChar("R"));

      Rf_namesgets(out,names);

      PutRNGstate();

      deleteMinHeap(&sk[0]);
      deleteMinHeap(&sk[1]);

      UNPROTECT(7);

      Free(t_out);
      Free(S_out);
      Free(E_out);
      Free(I_out);
      Free(R_out);

      return out;
    }
    // check to return early
    int all_zero = 1;
    for(int i=0; i<7; i++){
      if(i < 2){
        if((ak[i] > 1.E-9) || (sk[i].elem->data < DBL_MAX)){
          all_zero = 0;
        }
      } else {
        if(ak[i] > 1.E-9){
          all_zero = 0;
        }
      }
    }
    if(all_zero){
      Rprintf(" --- warning: all propensities approximately zero and no delayed reactions queued up, returning output early --- \n");

      SEXP out = PROTECT(Rf_allocVector(VECSXP,5));

      SEXP t_2R = PROTECT(Rf_allocVector(REALSXP,out_i-1));
      memmove(REAL(t_2R),t_out,(out_i-1) * sizeof(double));

      SEXP S_2R = PROTECT(Rf_allocVector(INTSXP,out_i-1));
      memmove(INTEGER(S_2R),S_out,(out_i-1) * sizeof(int));

      SEXP E_2R = PROTECT(Rf_allocVector(INTSXP,out_i-1));
      memmove(INTEGER(E_2R),E_out,(out_i-1) * sizeof(int));

      SEXP I_2R = PROTECT(Rf_allocVector(INTSXP,out_i-1));
      memmove(INTEGER(I_2R),I_out,(out_i-1) * sizeof(int));

      SEXP R_2R = PROTECT(Rf_allocVector(INTSXP,out_i-1));
      memmove(INTEGER(R_2R),R_out,(out_i-1) * sizeof(int));

      SET_VECTOR_ELT(out,0,t_2R);
      SET_VECTOR_ELT(out,1,S_2R);
      SET_VECTOR_ELT(out,2,E_2R);
      SET_VECTOR_ELT(out,3,I_2R);
      SET_VECTOR_ELT(out,4,R_2R);

      SEXP names = PROTECT(Rf_allocVector(STRSXP,5));
      SET_STRING_ELT(names,0,Rf_mkChar("time"));
      SET_STRING_ELT(names,1,Rf_mkChar("S"));
      SET_STRING_ELT(names,2,Rf_mkChar("E"));
      SET_STRING_ELT(names,3,Rf_mkChar("I"));
      SET_STRING_ELT(names,4,Rf_mkChar("R"));

      Rf_namesgets(out,names);

      PutRNGstate();

      deleteMinHeap(&sk[0]);
      deleteMinHeap(&sk[1]);

      UNPROTECT(7);

      Free(t_out);
      Free(S_out);
      Free(E_out);
      Free(I_out);
      Free(R_out);

      return out;

    }
    if(out_i >= out_size){
      if(verbose){
        Rprintf(" --- extending output memory --- \n");
      }

      int new_size = out_size+1E3;
      if(new_size > maxsize){
        new_size = maxsize - out_size;
      }

      // Rprintf(" --- NEW SIZE %d --- \n",new_size);

      double* t_out_new = (double*)Realloc(t_out,new_size,double);
      int* S_out_new = (int*)Realloc(S_out,new_size,int);
      int* E_out_new = (int*)Realloc(E_out,new_size,int);
      int* I_out_new = (int*)Realloc(I_out,new_size,int);
      int* R_out_new = (int*)Realloc(R_out,new_size,int);

      if(t_out_new == NULL){

        Rprintf(" --- ERROR: REALLOCATION FAILURE, ABORTING SIMULATION --- \n");
        Free(t_out);
        Free(S_out);
        Free(E_out);
        Free(I_out);
        Free(R_out);
        return R_NilValue;

      } else {
        t_out = t_out_new;
      }

      if(S_out_new == NULL){

        Rprintf(" --- ERROR: REALLOCATION FAILURE, ABORTING SIMULATION --- \n");
        Free(t_out);
        Free(S_out);
        Free(E_out);
        Free(I_out);
        Free(R_out);
        return R_NilValue;

      } else {
        S_out = S_out_new;
      }

      if(E_out_new == NULL){

        Rprintf(" --- ERROR: REALLOCATION FAILURE, ABORTING SIMULATION --- \n");
        Free(t_out);
        Free(S_out);
        Free(E_out);
        Free(I_out);
        Free(R_out);
        return R_NilValue;

      } else {
        E_out = E_out_new;
      }

      if(I_out_new == NULL){

        Rprintf(" --- ERROR: REALLOCATION FAILURE, ABORTING SIMULATION --- \n");
        Free(t_out);
        Free(S_out);
        Free(E_out);
        Free(I_out);
        Free(R_out);
        return R_NilValue;

      } else {
        I_out = I_out_new;
      }

      if(R_out_new == NULL){

        Rprintf(" --- ERROR: REALLOCATION FAILURE, ABORTING SIMULATION --- \n");
        Free(t_out);
        Free(S_out);
        Free(E_out);
        Free(I_out);
        Free(R_out);
        return R_NilValue;

      } else {
        R_out = R_out_new;
      }

      out_size = new_size;

    }

  }
  // end simulation loop

  if(verbose){
    Rprintf(" --- ending simulation --- \n");
  }

  SEXP out = PROTECT(Rf_allocVector(VECSXP,5));

  SEXP t_2R = PROTECT(Rf_allocVector(REALSXP,out_i-1));
  memmove(REAL(t_2R),t_out,(out_i-1) * sizeof(double));

  SEXP S_2R = PROTECT(Rf_allocVector(INTSXP,out_i-1));
  memmove(INTEGER(S_2R),S_out,(out_i-1) * sizeof(int));

  SEXP E_2R = PROTECT(Rf_allocVector(INTSXP,out_i-1));
  memmove(INTEGER(E_2R),E_out,(out_i-1) * sizeof(int));

  SEXP I_2R = PROTECT(Rf_allocVector(INTSXP,out_i-1));
  memmove(INTEGER(I_2R),I_out,(out_i-1) * sizeof(int));

  SEXP R_2R = PROTECT(Rf_allocVector(INTSXP,out_i-1));
  memmove(INTEGER(R_2R),R_out,(out_i-1) * sizeof(int));

  SET_VECTOR_ELT(out,0,t_2R);
  SET_VECTOR_ELT(out,1,S_2R);
  SET_VECTOR_ELT(out,2,E_2R);
  SET_VECTOR_ELT(out,3,I_2R);
  SET_VECTOR_ELT(out,4,R_2R);

  SEXP names = PROTECT(Rf_allocVector(STRSXP,5));
  SET_STRING_ELT(names,0,Rf_mkChar("time"));
  SET_STRING_ELT(names,1,Rf_mkChar("S"));
  SET_STRING_ELT(names,2,Rf_mkChar("E"));
  SET_STRING_ELT(names,3,Rf_mkChar("I"));
  SET_STRING_ELT(names,4,Rf_mkChar("R"));

  Rf_namesgets(out,names);

  PutRNGstate();

  deleteMinHeap(&sk[0]);
  deleteMinHeap(&sk[1]);

  UNPROTECT(7);

  Free(t_out);
  Free(S_out);
  Free(E_out);
  Free(I_out);
  Free(R_out);

  return out;

};
