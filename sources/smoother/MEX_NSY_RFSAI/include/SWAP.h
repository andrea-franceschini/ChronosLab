

#pragma once


template <typename TYPE>
inline void SWAP(TYPE & __restrict__ i1, TYPE & __restrict__ i2){TYPE tmp = i1; i1 = i2; i2 = tmp;}


template void SWAP<int>(int & __restrict__ , int & __restrict__ );
