/* --------------------------------------------------------------------------------
#
#   Discretize continuous time output
#   Sean Wu (slwu89@berkeley.edu)
#   March 2020
#
-------------------------------------------------------------------------------- */

#include "discretize.h"

SEXP discretize_C(SEXP out, SEXP dt_r){

  int dt = Rf_asInteger(dt_r);
  double* out_ptr = REAL(out);

  int ncol = Rf_ncols(out);
  int events = Rf_nrows(out);

  double end = out_ptr[(events-1) + events*0];
  double start = out_ptr[0 + events*0];

  int len = ((int)(end - start)) / dt + 1;
  SEXP x = PROTECT(Rf_allocMatrix(REALSXP,len,ncol));
  double* x_ptr = REAL(x);
  double target = 0.;
  int j = 0;

  for(int i=1; i<events; i++){
    while(out_ptr[i + events*0] >= target){
      x_ptr[j + len*0] = target;
      for(int k=1; k<ncol; k++){
        x_ptr[j + len*k] = out_ptr[i + events*k];
      }
      j += 1;
      target += dt;
    }
  }

  UNPROTECT(1);
  return x;
};
