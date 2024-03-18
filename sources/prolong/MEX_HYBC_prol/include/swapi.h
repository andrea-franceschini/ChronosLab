

#pragma once

typedef int iReg;


inline void swapi(iReg & __restrict__ i1, iReg & __restrict__ i2){iReg tmp = i1; i1 = i2; i2 = tmp;}

