

#pragma once

#include <cmath>

#include "precision.h"
#include "linsol_error.h"

#include "SWAP.h"

using namespace std;


void maxVol_inner(const iReg itmax, const rExt delta, const iReg n, const iReg r, iReg *list, rExt** Z);
