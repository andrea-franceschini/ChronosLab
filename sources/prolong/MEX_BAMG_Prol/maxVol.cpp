
#include "maxVol.h"





void maxVol(const iReg mmax, const rExt condmax, const iReg itmax, const rExt delta,
            const iReg n, const iReg m, rExt **A, iReg &rank, iReg *list){

   iReg maxrank = min(n,m);
   maxrank = min(maxrank,mmax);
   rExt wrot[m];


   rank = 1;








   iReg irow;
   rExt max_sv_est;
   move_row_front(0,n,m,A,irow,max_sv_est);
   SWAP(list[0],list[irow]);


   mk_HouHolVec(m,A[0],wrot);


   for (iReg i=0; i<n; i++) Apply_HouHol_Rot(m,wrot,A[i]);

   while (rank < maxrank){

      rExt norm;


      rank += 1;


      move_row_front(rank-1,n-rank+1,m,&(A[rank-1]),irow,norm);
      SWAP(list[rank-1],list[irow+rank-1]);


      rExt cond_new = max_sv_est / abs(norm);
      if (cond_new > condmax){
         rank -= 1;
         break;
      }


      mk_HouHolVec(m-rank+1,&(A[rank-1][rank-1]),wrot);


      for (iReg i=rank-1; i<n; i++) Apply_HouHol_Rot(m-rank+1,wrot,&(A[i][rank-1]));

   }







   backsolve(rank,n-rank,&(A[0]),&(A[rank]));


   try {
      maxVol_inner(itmax,delta,n-rank,rank,list,&(A[rank]));
   } catch (linsol_error) {
      throw linsol_error ("maxVol","permuting row in order to find the best basis");
   }

}
