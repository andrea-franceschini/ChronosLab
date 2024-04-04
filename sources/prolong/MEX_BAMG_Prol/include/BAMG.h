
#pragma once

#include <vector>

#include "precision.h"
#include "BAMG_params.h"


int BAMG ( const BAMG_params &params, const iReg nthreads, iReg nn_L, iReg nn_C,
           iReg nn_S, const iExt *const iat_S, const iReg *const ja_S, const iReg ntv,
           const iReg *const fcnodes, const rExt *const *const TV,
           iExt &nt_I, std::vector<iExt>& vec_iat_I, std::vector<iReg>& vec_ja_I,
           std::vector<rExt>& vec_coef_I, std::vector<iReg> &vec_c_mark );
