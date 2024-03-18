#include <stdlib.h>
#include <omp.h>
#include <limits>
#include "cblas.h"
#include "lapacke.h"

#include "inl_blas1.h"
#include <iostream>
#include <stdio.h>

//using namespace std;

#include "parm_EMIN.h"

#define PRINT_LOC_INFO 0
#define PRINT_LOC_INFO_ALL 0


int gather_B_QR(const int np, const double condmax, const double maxwgt, const int nn,
                const int nn_C, const int ntv, const int *fcnode, const int *iat_patt,
                const int *ja_patt, const double *const *TV, double *&mat_Q, double *coef_P0)
{
   int ierr = 0;

   int nrows_Q = iat_patt[nn];
   int nterm_Q = ntv*nrows_Q;
   mat_Q = (double*) malloc( nterm_Q*sizeof(double) );
   if (mat_Q == nullptr) return ierr = 3;


   int *c2glo = (int*) malloc( nn_C*sizeof(int) );
   double *vec_g = (double*) malloc( nrows_Q*sizeof(double) );
   if (c2glo == nullptr || vec_g == nullptr)  return ierr = 1;

   #pragma omp parallel num_threads(np)
   {

      int mythid = omp_get_thread_num();
      int bsize = nn/np;
      int resto = nn%np;
      int firstcol, ncolth, lastcol;
      if (mythid <= resto) {
         ncolth = bsize+1;
         firstcol = mythid*ncolth;
         if (mythid == resto) ncolth--;
      } else {
         ncolth = bsize;
         firstcol = mythid*bsize + resto;
      }
      lastcol = firstcol + ncolth;


      int pos_g = iat_patt[firstcol];
      int pos_Q = pos_g*ntv;


      double *g_scr = &(vec_g[pos_g]);
      double *BB_scr = &(mat_Q[pos_Q]);


      int nrmax_blk = 0;
      int istart, iend;
      iend = iat_patt[firstcol];
      for (int i = firstcol+1; i <= lastcol; i++){
         istart = iend;
         iend = iat_patt[i];
         nrmax_blk = std::max(nrmax_blk,iend-istart);
      }


      lapack_int ierr_lapack;
      lapack_int l_nn = static_cast<lapack_int>(std::max(ntv,nrmax_blk));
      lapack_int l_mm = static_cast<lapack_int>(ntv);
      lapack_int l_kk = static_cast<lapack_int>(ntv);
      double query_work_1;
      double query_work_2;
      double query_work_3;
      double *SIGMA = nullptr;
      double *dummy_U = nullptr;
      double *VT = nullptr;
      double *tau = nullptr;
      double *work = nullptr;

      lapack_int lwork = -1;
      ierr_lapack = LAPACKE_dgeqrf_work(LAPACK_COL_MAJOR,l_nn,l_mm,BB_scr,l_nn,tau,
                                        &query_work_1,lwork);
      if (ierr_lapack != 0){
         #pragma omp atomic write
         ierr = 4;
      }
      ierr_lapack = LAPACKE_dorgqr_work(LAPACK_COL_MAJOR,l_nn,l_mm,l_kk,BB_scr,l_nn,tau,
                                        &query_work_2,lwork);
      if (ierr_lapack != 0){
         #pragma omp atomic write
         ierr = 4;
      }
      ierr_lapack = LAPACKE_dgesvd_work(LAPACK_COL_MAJOR,'o','s',l_nn,l_mm,BB_scr,l_nn,
                                        SIGMA,dummy_U,l_nn,VT,l_mm,&query_work_3,lwork);
      if (ierr_lapack != 0){
         #pragma omp atomic write
         ierr = 4;
      }
      if (ierr > 0) goto exit_pragma;
      query_work_1 = std::max(query_work_1,query_work_2);
      query_work_3 = std::max(query_work_3,static_cast<double>(ntv));
      lwork = static_cast<lapack_int>(std::max(query_work_1,query_work_3));


      SIGMA = (double*) malloc( std::max(ntv,nrmax_blk)*sizeof(double) );
      VT = (double*) malloc( ntv*ntv*sizeof(double) );
      tau = (double*) malloc( (ncolth*ntv)*sizeof(double) );
      work = (double*) malloc( lwork*sizeof(double) );
      if (SIGMA == nullptr || VT == nullptr || tau == nullptr || work == nullptr){
         #pragma omp atomic write
         ierr = 2;
      }
      if (ierr > 0) goto exit_pragma;


      #pragma omp for
      for (int i = 0; i < nn; i++){
         int k = fcnode[i];
         if (k >= 0) c2glo[k] = i;
      }


      int ind_g, ind_BB, ind_tau, nnz_BB;
      ind_g = 0;
      ind_BB = 0;
      ind_tau = 0;
      int istart_patt, iend_patt;
      iend_patt = iat_patt[firstcol];
      for (int icol = firstcol; icol < lastcol; icol++){

         if (fcnode[icol] < 0){

            istart_patt = iend_patt;
            iend_patt = iat_patt[icol+1];
            int nr_BB_loc = iend_patt-istart_patt;


            if (nr_BB_loc > 0){

               int k = ind_BB;
               for (int i = istart_patt; i < iend_patt; i++){
                  int i_F = c2glo[ja_patt[i]];
                  int kk = k;
                  for (int j = 0; j < ntv; j++){
                     BB_scr[kk] = TV[i_F][j];
                     kk += nr_BB_loc;
                  }
                  k++;
               }

               for (int j = 0; j < ntv; j++) g_scr[ind_g+j] = TV[icol][j];


               cblas_dgemv(CblasColMajor,CblasTrans,nr_BB_loc,ntv,-1.0,&(BB_scr[ind_BB]),
                           nr_BB_loc,&(coef_P0[istart_patt]),1,1.0,&(g_scr[ind_g]),1);

               bool FAIL_QR = false;
               int i_lpk_1, i_lpk_2, i_lpk_3;
               if (nr_BB_loc >= ntv){


                  l_nn = static_cast<lapack_int>(nr_BB_loc);
                  lapack_int l_ll = l_nn;
                  i_lpk_1 = LAPACKE_dgeqrf_work(LAPACK_COL_MAJOR,l_nn,l_mm,
                                    &(BB_scr[ind_BB]),l_ll,&(tau[ind_tau]),work,lwork);


                  double max_DR = 0.0;
                  double min_DR = std::numeric_limits<double>::max();
                  for (int kk = 0; kk < l_mm; kk++){
                     max_DR = std::max(max_DR,fabs(BB_scr[ind_BB+kk*l_ll+kk]));
                     min_DR = std::min(min_DR,fabs(BB_scr[ind_BB+kk*l_ll+kk]));
                  }

                  if (max_DR / min_DR > condmax){

                     FAIL_QR = true;

                     int k = ind_BB;
                     for (int i = istart_patt; i < iend_patt; i++){
                        int i_F = c2glo[ja_patt[i]];
                        int kk = k;
                        for (int j = 0; j < ntv; j++){
                           BB_scr[kk] = TV[i_F][j];
                           kk += nr_BB_loc;
                        }
                        k++;
                     }

                     //if ( PRINT_LOC_INFO ) cout << icol <<
                     //   " conditioning larger than threshold: "
                     //   << max_DR / min_DR << " > " << condmax << endl;

                  } else {

                     i_lpk_2 = LAPACKE_dtrtrs_work(LAPACK_COL_MAJOR,'U','T','N',l_mm,1,
                                       &(BB_scr[ind_BB]),l_ll,&(g_scr[ind_g]),l_mm);


                     i_lpk_3 = LAPACKE_dorgqr_work(LAPACK_COL_MAJOR,l_nn,l_mm,l_kk,
                                       &(BB_scr[ind_BB]),l_ll,&(tau[ind_tau]),work,lwork);
                     if (i_lpk_1 || i_lpk_2 || i_lpk_3){
                        #pragma omp atomic write
                        ierr = 4;
                        goto exit_loop_icol;
                     }

                     cblas_dgemv(CblasColMajor,CblasNoTrans,nr_BB_loc,ntv,1.0,&(BB_scr[ind_BB]),
                                 nr_BB_loc,&(g_scr[ind_g]),1,1.0,&(coef_P0[istart_patt]),1);


                     nnz_BB = nr_BB_loc*ntv;

                  }

               }

               if ( (nr_BB_loc < ntv) || FAIL_QR){

                  //if ( PRINT_LOC_INFO )
                  //   cout << icol << " Solution for RANK deficient system " << endl;


                  l_nn = static_cast<lapack_int>(nr_BB_loc);
                  lapack_int l_ll = l_nn;
                  i_lpk_1 = LAPACKE_dgesvd_work(LAPACK_COL_MAJOR,'o','s',l_nn,l_mm,
                                                &(BB_scr[ind_BB]),l_ll,SIGMA,dummy_U,
                                                l_ll,VT,l_mm,work,lwork);
                  if (i_lpk_1){
                     #pragma omp atomic write
                     ierr = 4;
                     goto exit_loop_icol;
                  }


                  int rank_BB = 1;
                  if (nr_BB_loc > 1){
                     while (SIGMA[0] < condmax *abs(SIGMA[rank_BB])){
                        rank_BB++;
                        if (rank_BB == std::min(ntv,nr_BB_loc)){
                           break;
                        }
                     }
                  }

                  cblas_dgemv(CblasColMajor,CblasNoTrans,rank_BB,ntv,1.0,VT,
                              ntv,&(g_scr[ind_g]),1,0.0,work,1);


                  for (int i = 0; i < rank_BB; i++) work[i] /= SIGMA[i];


                  cblas_dgemv(CblasColMajor,CblasNoTrans,nr_BB_loc,rank_BB,1.0,
                              &(BB_scr[ind_BB]),nr_BB_loc,work,1,1.0,
                              &(coef_P0[istart_patt]),1);

                  int npad = (ntv-rank_BB)*nr_BB_loc;
                  for (int k = 0; k < npad; k++) BB_scr[ind_BB+nr_BB_loc*rank_BB+k] = 0.0;

                  nnz_BB = nr_BB_loc*ntv;

               }

               ind_BB += nnz_BB;
               ind_g += ntv;
               ind_tau += ntv;

            }

         }
      }
      exit_loop_icol: ;

      free(work);
      free(tau);
      free(VT);
      free(SIGMA);
      exit_pragma: ;

   }

   free(c2glo);
   free(vec_g);

   return ierr;

}
