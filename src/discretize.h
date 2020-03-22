/* --------------------------------------------------------------------------------
#
#   Discretize continuous time output
#   Sean Wu (slwu89@berkeley.edu)
#   March 2020
#
-------------------------------------------------------------------------------- */

#ifndef DISCRETIZE
#define DISCRETIZE

#include <stdlib.h>

#include <R.h>
#include <Rinternals.h>
#include <Rdefines.h>

SEXP discretize_C(SEXP out, SEXP dt_r);

#endif
