

#pragma once

#include "precision.h"
#include "linsol_error.h"


void add_new_neighs(const iExt *const iat_S, const iReg *const ja_S,
                    const iReg istart_neigh, const iReg iend_neigh,
                    const iReg nmax, iReg &n_neigh, iReg *neigh, iReg *WI);
