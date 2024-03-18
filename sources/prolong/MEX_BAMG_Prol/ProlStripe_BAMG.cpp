#include <vector>
#include <math.h>
#include <algorithm>
#include <lapacke.h>
using namespace std;

#include "precision.h"
#include "BAMG_params.h"
#include "linsol_error.h"

#include "check_neigh_cond.h"
#include "add_new_neighs.h"
#include "inl_blas1.h"
#include "maxVol.h"
#include "heapsort.h"
#include "Reduce_IntSet.h"



const rExt ONE = 1.0;


void ProlStripe_BAMG(const BAMG_params& params, iReg firstrow_0, iReg firstrow, iReg lastrow,
                     iReg nn_S, iReg ntvecs, const iExt *const iat_S, const iReg *const ja_S,
                     const iReg *const fcnodes,
                     const rExt *const *const TV, iExt &nt_P, iExt *iat_P, iReg *ja_P,
                     rExt *coef_P, iReg *c_mark, iReg *dist_count){


   bool VERB_FLAG = params.verbosity >= VLEV_MEDIUM;


   iReg itmax_vol = params.itmax_vol;
   iReg dist_min  = params.dist_min;
   iReg dist_max  = params.dist_max;
   iReg mmax      = params.mmax;
   rExt maxcond   = params.maxcond;
   rExt maxrownrm = params.maxrownrm;
   rExt tol_vol   = params.tol_vol;
   rExt eps       = params.eps;


   rExt *dummy_double = nullptr;
   iReg optimal_lwork;
   rExt db_lwork;
   lapack_int info = LAPACKE_dgels_work(LAPACK_COL_MAJOR,'N',ntvecs,ntvecs,1,
                                        dummy_double,ntvecs,dummy_double,ntvecs,
                                        &db_lwork,-1);
   if (info != 0)
      throw linsol_error ("ProlStripe_BAMG","Quering LAPACK DGELS workspace");
   optimal_lwork = static_cast<int>(db_lwork);


   vector<iReg> vec_int_list;
   vector<iReg> vec_neigh;
   vector<iReg> vec_WI;
   vector<rExt> vec_WR;
   try {
      vec_int_list.resize(nn_S);
      vec_neigh.resize(nn_S);
      vec_WI.resize(nn_S);
      vec_WR.resize(optimal_lwork);
   } catch (linsol_error) {
      throw linsol_error ("ProlStripe_BAMG","allocating scratches");
   }
   iReg *int_list = vec_int_list.data();
   iReg *neigh = vec_neigh.data();
   iReg *WI = vec_WI.data();
   rExt *WR = vec_WR.data();


   vector<rExt*> vec_TVcomp;
   vector<rExt>  vec_TVbuf;
   try {
      vec_TVcomp.resize(nn_S+1);
      vec_TVbuf.resize(ntvecs*(nn_S+1));
   } catch (linsol_error) {
      throw linsol_error ("ProlStripe_BAMG","allocating compressed TV");
   }
   rExt **TVcomp = vec_TVcomp.data();
   rExt  *TVbuf  = vec_TVbuf.data();

   iReg kk = 0;
   for (iReg i = 0; i < nn_S+1; i++){
      TVcomp[i] = &(TVbuf[kk]);
      kk += ntvecs;
   }


   for (iReg i = 0; i < nn_S; i++) WI[i] = 0;


   iExt ind_P = 0;
   iat_P[0] = ind_P;


   iReg shift = firstrow - 1;
   for (iReg inod = firstrow; inod < lastrow ; inod ++){

      iReg inod_coarse = fcnodes[inod];


      if (inod_coarse >= 0){


         ja_P[ind_P] = inod_coarse;
         coef_P[ind_P] = ONE;
         ind_P++;

      } else {


         iReg row_rank;


         iReg distance = 0;
         iReg n_neigh = 1;
         neigh[0] = inod;
         WI[inod] = 1;


         iReg istart_neigh = 0;
         iReg iend_neigh = 1;
         while (distance < dist_min){
            distance++;
            try {
               add_new_neighs(iat_S,ja_S,istart_neigh,iend_neigh,nn_S,n_neigh,neigh,WI);
            } catch (linsol_error) {
               throw linsol_error ("ProlStripe_BAMG",
                                   "add the list of neighbour up to distance dist_min");
            }
            if (n_neigh == iend_neigh){

               break;
            }
            istart_neigh = iend_neigh;
            iend_neigh = n_neigh;
         }


         rExt bnorm = inl_dnrm2(ntvecs,TV[inod],1);
         rExt res_ass = eps*bnorm;
         rExt res_row;
         rExt row_nrm;

         bool UPD_WEIGHTS = true;
         while ( UPD_WEIGHTS ){


            iReg n_int = 0;
            for (iReg i = 1; i < n_neigh; i++){
               iReg i_neigh = neigh[i];
               if (fcnodes[i_neigh] >= 0){
                  int_list[n_int] = i_neigh;
                  for (iReg j = 0; j < ntvecs; j++) TVcomp[n_int][j] = TV[i_neigh][j];
                  n_int++;
               }
            }

            if (n_int > 0){


               maxVol(mmax,maxcond,itmax_vol,tol_vol,n_int,ntvecs,TVcomp,row_rank,int_list);

               heapsort(int_list,row_rank);


               for (iReg i = 0; i < row_rank; i++){
                  iReg i_neigh = int_list[i];
                  for (iReg j = 0; j < ntvecs; j++) TVcomp[i][j] = TV[i_neigh][j];
               }


               for (iReg j = 0; j < ntvecs; j++) coef_P[ind_P+j] = TV[inod][j];


               info = LAPACKE_dgels_work(LAPACK_COL_MAJOR,'N',ntvecs,row_rank,1,
                                         &(TVcomp[0][0]),ntvecs,&(coef_P[ind_P]),ntvecs,
                                         WR,optimal_lwork);
               if(info != 0){throw linsol_error ("ProlStripe_BAMG","error in LAPACKE_dgels");}


               res_row = inl_dnrm2(ntvecs-row_rank,&(coef_P[ind_P+row_rank]),1);


               row_nrm = inl_dnrm2(row_rank,&(coef_P[ind_P]),1);

            } else {


               row_rank = 0;
               res_row = bnorm;
               row_nrm = 0.0;

            }


            if ( row_nrm < maxrownrm && ((res_row < res_ass) || (row_rank == ntvecs)) ){



               UPD_WEIGHTS = false;

               if (VERB_FLAG){

                  #pragma omp atomic update
                  dist_count[distance-1]++;
               }


            } else {


               if (distance < dist_max){
                  distance++;
                  try {
                     add_new_neighs(iat_S,ja_S,istart_neigh,iend_neigh,nn_S,n_neigh,neigh,WI);
                  } catch (linsol_error) {
                     throw linsol_error ("ProlStripe_BAMG","try to increase distance if possible");
                  }
                  if (n_neigh == iend_neigh){

                     UPD_WEIGHTS = false;
                     if (n_int == 0){
                        if (VERB_FLAG){

                           #pragma omp atomic update
                           dist_count[dist_max+2]++;
                        }

                     } else {


                        if (row_nrm > RELAX_FAC*maxrownrm){
                           Reduce_IntSet(RELAX_FAC*maxrownrm,itmax_vol,tol_vol,maxcond,
                                         optimal_lwork,inod,n_neigh,ntvecs,row_rank,
                                         fcnodes,int_list,neigh,TV,TVcomp,WR,
                                         &(coef_P[ind_P]),row_nrm);
                        }
                        res_row = inl_dnrm2(ntvecs-row_rank,&(coef_P[ind_P+row_rank]),1);
                        if (res_row < res_ass){
                           if (VERB_FLAG){

                              #pragma omp atomic update
                              dist_count[dist_max+1]++;
                           }
                        } else {

                           c_mark[inod] = 1;
                           if (VERB_FLAG){

                              #pragma omp atomic update
                              dist_count[dist_max]++;
                           }
                        }
                     }

                  }
                  istart_neigh = iend_neigh;
                  iend_neigh = n_neigh;

               } else {


                  UPD_WEIGHTS = false;


                  if (row_nrm > RELAX_FAC*maxrownrm)
                     Reduce_IntSet(RELAX_FAC*maxrownrm,itmax_vol,tol_vol,maxcond,
                                   optimal_lwork,inod,n_neigh,ntvecs,row_rank,fcnodes,
                                   int_list,neigh,TV,TVcomp,WR,&(coef_P[ind_P]),row_nrm);

                  res_row = inl_dnrm2(ntvecs-row_rank,&(coef_P[ind_P+row_rank]),1);

                  if (res_row > res_ass) c_mark[inod] = 1;

                  if (VERB_FLAG){
                     if (res_row > res_ass){

                        #pragma omp atomic update
                        dist_count[dist_max]++;
                     } else {

                        #pragma omp atomic update
                        dist_count[dist_max+1]++;
                     }
                  }


               }

            }

         }


         for (iReg i = 0; i < row_rank; i++){
            ja_P[ind_P] = fcnodes[int_list[i]];
            ind_P++;
         }


         for (iReg i = 0; i < n_neigh; i++) WI[neigh[i]] = 0;

      }


      iat_P[inod-shift] = ind_P;

   }


   nt_P = ind_P;

}
