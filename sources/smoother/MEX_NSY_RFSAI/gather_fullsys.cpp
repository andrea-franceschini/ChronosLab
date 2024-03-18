#include "gather_fullsys.h"

void gather_fullsys(const int irow, const int mrow, const int *const vecinc,
                    const int nn, const int *const iat, const int *const ja,
                    const double *const coef, double *full_A, double *rhs_L,
                    double *rhs_U, bool &null_L, bool &null_U){


   null_L = true;
   null_U = true;
   int ind_row = 0;
   for (int i = 0; i < mrow; i++){

      int ii = 0;
      int krow = vecinc[i];
      int jj = iat[krow];
      int endrow = iat[krow+1];
      while (ii < mrow){

         while (ja[jj] < vecinc[ii]){
            jj++;
            if (jj == endrow){

               for (int k = ii; k < mrow; k++) full_A[ind_row+k] = 0.0;
               rhs_U[i] = 0.0;
               goto next_row;
            }
         }
         if (vecinc[ii] == ja[jj]){

            full_A[ind_row+ii] = coef[jj];
            ii++;
         } else {

            full_A[ind_row+ii] = 0.0;
            ii++;
         }
      }







      while(ja[jj] < irow){
         jj++;
         if (jj == endrow){
            rhs_U[i] = 0.0;
            goto next_row;
         }
      }
      if (irow == ja[jj]){

         rhs_U[i] = -coef[jj];
         null_U = false;
      } else {

         rhs_U[i] = 0.0;
      }
      next_row:;
      ind_row += mrow;
   }


   int ii = 0;
   int diag_row = irow;
   int jj = iat[diag_row];
   int endrow = iat[diag_row+1];
   while (ii < mrow){

      while (ja[jj] < vecinc[ii]){
         jj++;
         if (jj == endrow){

            for (int k = ii; k < mrow; k++) rhs_L[k] = 0.0;
            goto end_rhs_L;
         }
      }
      if (vecinc[ii] == ja[jj]){

         rhs_L[ii] = -coef[jj];
         null_L = false;
         ii++;
      } else {

         rhs_L[ii] = 0.0;
         ii++;
      }
   }
   end_rhs_L:;

}
