

#pragma once

#include "precision.h"
#include "inl_blas1.h"


void move_row_front(const iReg jcol, const iReg n, const iReg m, rExt **A,
                    iReg &irow, rExt &maxnorm);
