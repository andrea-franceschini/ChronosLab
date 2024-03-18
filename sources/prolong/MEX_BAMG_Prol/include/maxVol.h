

#pragma once

#include <cmath>
#include <algorithm>

#include "precision.h"
#include "linsol_error.h"

#include "move_row_front.h"
#include "mk_HouHolVec.h"
#include "Apply_HouHol_Rot.h"
#include "backsolve.h"
#include "SWAP.h"
#include "maxVol_inner.h"


void maxVol(const iReg mmax, const rExt condmax, const iReg itmax, const rExt delta,
            const iReg n, const iReg m, rExt **A, iReg &rank, iReg *list);


