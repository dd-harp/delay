/* --------------------------------------------------------------------------------
#
#   Modified Next Reaction Method (Anderson 2007) for SEIR w/delay
#   Sean Wu (slwu89@berkeley.edu)
#   March 2020
#
-------------------------------------------------------------------------------- */

#ifndef SEIR_DELAY_MNRM
#define SEIR_DELAY_MNRM

#include <stdlib.h>
#include <float.h>

#include <R.h>
#include <Rinternals.h>
#include <Rdefines.h>

#include <Rmath.h>


/* --------------------------------------------------------------------------------
#   SEIR w/delay via Anderson's method
-------------------------------------------------------------------------------- */
// as.integer(as.logical(verbose))
SEXP SEIRdelay_MNRM_C(SEXP tmax_r, SEXP S0_r, SEXP I0_r,SEXP beta_r,SEXP nu_r,SEXP tau_r, SEXP verbose_r, SEXP maxsize_r);

#endif
