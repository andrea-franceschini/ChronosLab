

#pragma once

#include "precision.h"


template <typename TYPE>
inline void SWAP(TYPE & __restrict__ i1, TYPE & __restrict__ i2){TYPE tmp = i1; i1 = i2; i2 = tmp;}


#if !IREG_LONG==IEXT_LONG
template void SWAP<iReg>(iReg & __restrict__ , iReg & __restrict__ );
#endif
template void SWAP<iExt>(iExt & __restrict__ , iExt & __restrict__ );
template void SWAP<rExt>(rExt & __restrict__ , rExt & __restrict__ );
