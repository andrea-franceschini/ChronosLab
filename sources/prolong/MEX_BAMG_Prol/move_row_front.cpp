#include "move_row_front.h"




void move_row_front(const iReg jcol, const iReg n, const iReg m, rExt **A,
                    iReg &irow, rExt &maxnorm){


   irow = 0;
   maxnorm = inl_dnrm2(m-jcol,&(A[0][jcol]),1);
   for (iReg i = 1; i < n; i++){
      rExt norm = inl_dnrm2(m-jcol,&(A[i][jcol]),1);
      if (norm > maxnorm){
         maxnorm = norm;
         irow = i;
      }
   }
   if (irow != 0){

      for (iReg j = 0; j < m; j++){
         rExt tmp = A[irow][j];
         A[irow][j] = A[0][j];
         A[0][j] = tmp;
      }
   }

   return;
}
