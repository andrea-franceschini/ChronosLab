#pragma once


void mkiat_Tloc (const int nrows, const int nequ, const int nthreads, const int firstrow,
                 int** __restrict__ WI, int* __restrict__ iat_T, int& nnz  );
