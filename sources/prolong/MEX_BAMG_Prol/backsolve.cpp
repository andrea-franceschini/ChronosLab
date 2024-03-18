
#include "backsolve.h"





void backsolve(const iReg n, const iReg m, const rExt *const *const R, rExt **Z){


   for (iReg i = n-1; i > -1; i--){

      for (iReg k = m-1; k > -1; k--){

         Z[k][i] /= R[i][i];

         for (iReg j = i-1; j > -1; j--){
            Z[k][j] -= R[i][j]*Z[k][i];
         }
      }
   }

   return;
}
