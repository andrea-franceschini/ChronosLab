#include "transpose.h"

int transpose(const int nrows, const int ncols, const int *const iat, const int *const ja,
              const double *const coef, int *&iat_T, int *&ja_T, double *&coef_T){


   iat_T = (int*) malloc((ncols+1) * sizeof(int));
   int nterm = iat[nrows];
   ja_T = (int*) malloc((nterm) * sizeof(int));
   coef_T = (double*) malloc((nterm) * sizeof(double));
   int *ISCR = (int*) malloc((ncols+1) * sizeof(int));
   if (iat_T == nullptr || ja_T == nullptr || coef_T == nullptr || ISCR == nullptr){

      cout << "Allocation Error in transpose" << endl;
      return 1;
   }


   fill_n(iat_T,ncols+1,0);


   for ( int i = 0; i < nrows; i++ ){
      for ( int j = iat[i]; j < iat[i+1]; j++ ) iat_T[ja[j]]++;
   }


   ISCR[0] = 0;
   for ( int i = 1; i < ncols+1; i++ ) ISCR[i] = ISCR[i-1] + iat_T[i-1];
   for ( int i = 0; i < ncols+1; i++ ) iat_T[i] = ISCR[i];


   for ( int i = 0; i < nrows; i++ ){
      for ( int j = iat[i]; j < iat[i+1]; j++ ){
         int ind  = ISCR[ja[j]];
         ja_T[ind] = i;
         coef_T[ind] = coef[j];
         ISCR[ja[j]] = ind+1;
      }
   }


   free(ISCR);

   return 0;

}
