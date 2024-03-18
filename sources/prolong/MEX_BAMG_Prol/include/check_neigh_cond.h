

#pragma once


#include <omp.h>


#include "precision.h"
#include "linsol_error.h"

#include "add_new_neighs.h"
#include "get_cond.h"


void check_neigh_cond(const rExt maxcond,const iReg inod, const iReg ntv, const iReg nn_S,
                      const iExt *const iat_S, const iReg *const ja_S,
                      const rExt *const *const TV,
                      rExt **TVcomp, iReg *neigh, iReg *WI, rExt &local_cond);

