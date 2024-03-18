#include "omp.h"

#include <iostream>
using namespace std;


typedef int iReg;
typedef int iExt;
typedef double rExt;

#define ZERO 0.0


void KP_spmat(const iReg nthreads, const iReg nrows_A, const iExt* iat_A,
              const iReg *ja_A, const rExt *coef_A, const iReg ncols_B,
              const iExt *iat_B, const iReg *ja_B, const rExt *coef_B, rExt *coef_C,
              iReg* WNALL){


   #pragma omp parallel num_threads(nthreads)
   {


      int mythid = omp_get_thread_num();


      iReg* WN = WNALL + mythid*ncols_B;


      for (iReg i = 0; i < ncols_B; i++){
         WN[i] = -1;
      }


      #pragma omp for
      for (iReg irow = 0; irow < nrows_A; irow++){


         iExt ind_C = iat_B[irow];


         iReg nn_row = 0;
         iExt Cstop = iat_B[irow+1];
         for (iExt j = ind_C; j < Cstop; j++){
            iReg jcol = ja_B[j];
            coef_C[ind_C+nn_row] = ZERO;
            WN[jcol] = nn_row;
            nn_row++;
         }


         iExt ISTRT = iat_A[irow];
         iExt ISTOP = iat_A[irow+1];
         for (iExt j = ISTRT; j < ISTOP; j++){

            iReg jjA = ja_A[j];
            rExt WRA = coef_A[j];


            iExt JSTRT = iat_B[jjA];
            iExt JSTOP = iat_B[jjA+1];
            for (iExt k = JSTRT; k < JSTOP; k++){


               iReg pos = WN[ja_B[k]];

               if (pos >= 0) {

                  coef_C[ind_C+pos] += WRA * coef_B[k];
               }

            }

         }


         for ( iReg j = 0; j < nn_row; j ++ ) {
            WN[ja_B[ind_C+j]] = -1;
         }

      }

   }

}
