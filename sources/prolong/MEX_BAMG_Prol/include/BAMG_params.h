

#pragma once


#define VLEV_NONE 0
#define VLEV_LOW 1
#define VLEV_MEDIUM 2
#define VLEV_HIGH 3


#define RELAX_FAC 1.5

#include "precision.h"


struct BAMG_params {int verbosity; rExt maxrownrm; rExt maxcond; iReg itmax_vol;
                    rExt tol_vol; rExt eps; iReg dist_min; iReg dist_max; iReg mmax; };
