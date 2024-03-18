

#pragma once

typedef double rExt;


inline void swapr(rExt & __restrict__ i1, rExt & __restrict__ i2){rExt tmp = i1; i1 = i2; i2 = tmp;}

