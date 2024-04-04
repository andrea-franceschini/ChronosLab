#pragma once

#include <math.h>
#include <algorithm>

#include "DEBUG.h"
#include "ri_sortsplit_nsy.h"
#include "heapsort.h"


void KapGrad_NSY(const int istep, const int irow, int& mrow, const int jendbloc,
                 const int lfil, const int* __restrict__ iat, const int* __restrict__ ja,
                 const double* __restrict__ coef, const double* __restrict__ coef_T,
                 const double* __restrict__ sol_L, const double* __restrict__ sol_U,
                 int* __restrict__ IWN, int* __restrict__ JWN,
                 double* __restrict__ WR_L, double* __restrict__ WR_U);
