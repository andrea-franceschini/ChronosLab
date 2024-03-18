

#pragma once

#include <cmath>
#include <algorithm>

#include "precision.h"

#include "move_row_front.h"
#include "mk_HouHolVec.h"
#include "Apply_HouHol_Rot.h"

using namespace std;


void get_cond(const rExt condmax, const iReg n, const iReg m, rExt **A,
              iReg &rank, rExt &cond);
