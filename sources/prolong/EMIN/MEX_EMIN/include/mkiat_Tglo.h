#pragma once

void mkiat_Tglo (const int myid, const int nrows, const int nequ, const int nthreads,
                 const int firstrow, int** __restrict__ WI,
                 const int* __restrict__  nnz, int* __restrict__ iat_T  );
