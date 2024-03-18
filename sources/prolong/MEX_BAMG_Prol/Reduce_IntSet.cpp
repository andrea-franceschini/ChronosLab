#include "lapacke.h"

#include "precision.h"
#include "maxVol.h"
#include "heapsort.h"




void Reduce_IntSet(const rExt maxrownrm, const iReg itmax_vol, const rExt tol_vol,
                   const rExt maxcond, const iReg optimal_lwork, const iReg inod,
                   const iReg n_neigh, const iReg ntvecs, iReg &row_rank,
                   const iReg *const fcnodes, iReg *int_list, const iReg *neigh,
                   const rExt *const *const TV, rExt **TVcomp, rExt *WR, rExt *coef_P,
                   rExt &row_nrm){

   while (row_nrm > maxrownrm) {
      if (row_rank == 1){
         row_rank--;
         break;
      }


      iReg n_int = 0;
      for (iReg i = 1; i < n_neigh; i++){
         iReg i_neigh = neigh[i];
         if (fcnodes[i_neigh] >= 0){
            int_list[n_int] = i_neigh;
            for (iReg j = 0; j < ntvecs; j++)
               TVcomp[n_int][j] = TV[i_neigh][j];
            n_int++;
         }
      }


      int cmax = row_rank-1;


      maxVol(cmax,maxcond,itmax_vol,tol_vol,n_int,ntvecs,TVcomp,
             row_rank,int_list);

      heapsort(int_list,row_rank);


      for (int i = 0; i < row_rank; i++){
         int i_neigh = int_list[i];
         for (int j = 0; j < ntvecs; j++)
            TVcomp[i][j] = TV[i_neigh][j];
      }


      for (int j = 0; j < ntvecs; j++)
         coef_P[j] = TV[inod][j];


      lapack_int info = LAPACKE_dgels_work(LAPACK_COL_MAJOR,'N',ntvecs,row_rank,1,&(TVcomp[0][0]),
                                ntvecs,coef_P,ntvecs,WR,optimal_lwork);
      if(info != 0){throw linsol_error ("ProlStripe_BAMG","error in LAPACKE_dgels");}


      row_nrm = inl_dnrm2(row_rank,coef_P,1);
   }
}
