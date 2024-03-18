

#pragma once


void ProlStripe_BAMG(const BAMG_params& params, iReg firstrow_0, iReg firstrow, iReg lastrow,
                     iReg nn_S, iReg ntvecs, const iExt *const iat_S, const iReg *const ja_S,
                     const iReg *const fcnodes,
                     const rExt *const *const TV, iExt &nt_P, iExt *iat_P, iReg *ja_P,
                     rExt *coef_P, iReg *c_mark, iReg *dist_count);
