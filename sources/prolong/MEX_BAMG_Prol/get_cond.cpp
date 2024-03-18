
#include "get_cond.h"





void get_cond(const rExt condmax, const iReg n, const iReg m, rExt **A,
              iReg &rank, rExt &cond){

   iReg maxrank = min(n,m);
   rExt wrot[m];


   rank = 1;



   iReg irow;
   rExt max_sv_est;
   cond = 1.0;
   move_row_front(0,n,m,A,irow,max_sv_est);


   mk_HouHolVec(m,A[0],wrot);


   for (iReg i=1; i<n; i++) Apply_HouHol_Rot(m,wrot,A[i]);

   while (rank < maxrank){

      rExt norm;


      rank += 1;


      move_row_front(rank-1,n-rank+1,m,&(A[rank-1]),irow,norm);


      rExt cond_new = max_sv_est / fabs(norm);
      if (cond_new > condmax){
         rank -= 1;
         break;
      }
      cond = cond_new;


      mk_HouHolVec(m-rank+1,&(A[rank-1][rank-1]),wrot);


      for (iReg i=rank; i<n; i++) Apply_HouHol_Rot(m-rank+1,wrot,&(A[i][rank-1]));

   }

}
