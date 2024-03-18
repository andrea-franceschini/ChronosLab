
#include "mkiat_Tloc.h"




void mkiat_Tloc (const int nrows, const int nequ, const int nthreads, const int firstrow,
                 int** __restrict__ WI, int* __restrict__ iat_T, int& nnz  ){


   for ( int i = 0; i < nrows; i++ ) {
      iat_T[i] = WI[1][firstrow+i];
   }
   for ( int j = 2; j < nthreads+1; j++ ) {
      for ( int i = 0; i < nrows; i++ ) {
         iat_T[i] += WI[j][firstrow+i];
      }
   }


   int tmp = iat_T[0];
   iat_T[0] = 0;
   for ( int i = 1; i < nrows; i++ ) {
      int tmp2 = iat_T[i];
      iat_T[i] = iat_T[i-1] + tmp;
      tmp = tmp2;
   }


   nnz = iat_T[nrows-1] + tmp;

}




