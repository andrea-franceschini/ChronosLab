
#include "maxVol_inner.h"



inline void find_zij(const iReg n, const iReg r, const rExt *const *const Z, iReg &irow,
                     iReg &jcol, rExt &zij){
   irow = 0;
   jcol = 0;
   zij = 0.0;
   for (iReg j = 0; j < r; j++){
      for (iReg i = 0; i < n; i++){
         if (abs(Z[i][j]) > zij){
            zij = abs(Z[i][j]);
            irow = i;
            jcol = j;
         }
      }
   }
  return;
}




void maxVol_inner(const iReg itmax, const rExt delta, const iReg n, const iReg r,
                  iReg *list, rExt** Z){

   iReg irow, jcol;
   rExt zij_max;
   rExt vcol[n];
   rExt vrow[r];


   iReg iter = 0;


   find_zij(n,r,Z,irow,jcol,zij_max);


   while ( (zij_max > 1.0+delta) & (iter < itmax) ){

      iter++;

      for (iReg i = 0; i < n; i++) vcol[i] = Z[i][jcol];
      vcol[irow] -= 1.0;
      rExt fac = 1.0 / Z[irow][jcol];
      for (iReg j = 0; j < r; j++) vrow[j] = Z[irow][j]*fac;
      vrow[jcol] += fac;
      for (iReg i = 0; i < n; i++){
         for (iReg j = 0; j < r; j++){
            Z[i][j] -=  vcol[i]*vrow[j];
         }
      }
      SWAP(list[irow+r],list[jcol]);

      find_zij(n,r,Z,irow,jcol,zij_max);
   }


   if (zij_max > 1.0+delta){
     throw linsol_error ("maxVol_inner","convergence not achieved");
   }
}
