#include "SymmetrizePattern.h"


int SymmetrizePattern(const int nrows, int *& iat, int *& ja, double *& coef,
                      double *&coef_T){


   int *iat_T = nullptr;
   int *ja_T = nullptr;
   int ierr = transpose(nrows,nrows,iat,ja,coef,iat_T,ja_T,coef_T);
   if (ierr != 0) return 1;


   int nterm = iat[nrows];
   int nterm_tmp = 2*nterm - nrows;
   int *iat_tmp = (int*) malloc((nrows+1) * sizeof(int));
   int *ja_tmp = (int*) malloc((nterm_tmp) * sizeof(int));
   double *coef_tmp = (double*) malloc((nterm_tmp) * sizeof(double));
   if (iat_tmp == nullptr || ja_tmp == nullptr || coef_tmp == nullptr) return 1;


   iat_tmp[0] = 0;
   int ind = 0;
   int iend   = iat[0];
   int iend_T = iat_T[0];
   for (int i = 0; i < nrows; i++){

      int len_out;

      int istart   = iend;
      iend     = iat[i+1];
      int len      = iend-istart;
      int istart_T = iend_T;
      iend_T   = iat_T[i+1];
      int len_T    = iend_T-istart_T;
      merge_row_patt(len,&(ja[istart]),&(coef[istart]),
                     len_T,&(ja_T[istart_T]),&(coef_T[istart_T]),
                     len_out,&(ja_tmp[ind]),&(coef_tmp[ind]));


      ind += len_out;
      iat_tmp[i+1] = ind;

   }


   free(iat_T);
   free(ja_T);
   free(coef_T);


   int *iat_sav = iat;
   int *ja_sav = ja;
   double *coef_sav = coef;


   nterm_tmp = ind;
   iat = iat_tmp;
   ja    = (int*) realloc( ja_tmp, (nterm_tmp) * sizeof(int));
   coef  = (double*) realloc( coef_tmp , (nterm_tmp) * sizeof(double));
   if ( ja == nullptr || coef == nullptr ) return 1;


   free(iat_sav);
   free(ja_sav);
   free(coef_sav);


   ierr = transpose(nrows,nrows,iat,ja,coef,iat_T,ja_T,coef_T);
   if (ierr != 0) return 1;


   free(iat_T);
   free(ja_T);

   return 0;

}
