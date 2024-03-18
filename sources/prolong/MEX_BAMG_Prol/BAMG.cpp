#include <omp.h>
#include <vector>
#include <iomanip>
//using namespace std;

#include "precision.h"
#include "BAMG_params.h"
#include "linsol_error.h"

#include "ProlStripe_BAMG.h"

int BAMG ( const BAMG_params &params, const iReg nthreads, iReg nn_L, iReg nn_C,
           iReg nn_S, const iExt *const iat_S, const iReg *const ja_S, const iReg ntv,
           const iReg *const fcnodes, const rExt *const *const TV,
           iExt &nt_I, vector<iExt> &vec_iat_I, vector<iReg> &vec_ja_I,
           vector<rExt> &vec_coef_I, vector<iReg> &vec_c_mark ) {


   int ierr = 0;


   vector<iExt> vec_ridv_i;
   try {
      vec_ridv_i.resize(nthreads+1);
   } catch (linsol_error) {
      linsol_error ("BAMG","allocationg ridv_i");
      return 1;
   }
   iExt* ridv_i = vec_ridv_i.data();


   try {
      vec_c_mark.assign(nn_C,0);
   } catch (linsol_error) {
      linsol_error ("BAMG","allocating c_mark");
      return 1;
   }
   iReg* c_mark = vec_c_mark.data();


   iReg len_count = params.dist_max+3;
   vector<iReg> dist_count;

   if (params.verbosity >= VLEV_MEDIUM){
      dist_count.assign(len_count,0);
   }




   #pragma omp parallel num_threads(nthreads)

   {


      int ierr_L = 0;


      iReg shift;
      iExt istart_I, istart_scr, iend_scr;
      iExt* iat_scr;
      iReg* ja_scr;
      rExt* coef_scr;
      iExt* iat_I;
      iReg* ja_I;
      rExt* coef_I;


      iReg myid = static_cast<iReg>(omp_get_thread_num());


      iReg firstrow_0;
      iReg firstrow;
      iReg nn_loc;
      iReg bsize = nn_C/nthreads;
      iReg resto = nn_C%nthreads;
      if (myid <= resto) {
         nn_loc = bsize+1;
         firstrow = myid*nn_loc;
         if (myid == resto) nn_loc--;
      } else {
         nn_loc = bsize;
         firstrow = myid*bsize + resto;
      }
      firstrow_0    = nn_L;
      firstrow     += nn_L;
      iReg lastrow  = firstrow + nn_loc;
      iExt nt_Imax  = max(params.mmax,ntv)*nn_loc;
      iExt nt_I_loc;


      vector<iExt> vec_iat_scr;
      vector<iReg> vec_ja_scr;
      vector<rExt> vec_coef_scr;
      try {
         vec_iat_scr.resize(nn_loc+1);
         vec_ja_scr.resize(nt_Imax);
         vec_coef_scr.resize(nt_Imax);
      } catch (linsol_error) {
         linsol_error ("BAMG","allocating local scratches");
         ierr_L = 1;
         goto exit_omp;
      }
      iat_scr  = vec_iat_scr.data();
      ja_scr   = vec_ja_scr.data();
      coef_scr = vec_coef_scr.data();


      try {
         ProlStripe_BAMG(params,firstrow_0,firstrow,lastrow,nn_S,ntv,iat_S,ja_S,
                         fcnodes,TV,nt_I_loc,iat_scr,ja_scr,coef_scr,c_mark,
                         dist_count.data());
      } catch (linsol_error) {
         linsol_error ("BAMG","in computing the local stripe of BAMG");
         ierr_L = 2;
         goto exit_omp;
      }


      ridv_i[myid+1] = nt_I_loc;
      #pragma omp barrier


      #pragma omp single
      {
         ridv_i[0] = 0;
         for (iReg ip = 0; ip < nthreads; ip++){
             ridv_i[ip+1] = ridv_i[ip] + ridv_i[ip+1];
         }
         nt_I   = ridv_i[nthreads];
         try {
            vec_iat_I.resize(nn_S+1);
            vec_ja_I.resize(nt_I);
            vec_coef_I.resize(nt_I);
         } catch (linsol_error) {
            linsol_error ("BAMG","allocating prolongation matrix");
            ierr = 1;
         }
      }
      if (ierr != 0) goto exit_omp;
      iat_I  = vec_iat_I.data();
      ja_I   = vec_ja_I.data();
      coef_I = vec_coef_I.data();


      istart_I = ridv_i[myid];
      shift = firstrow;
      iend_scr = iat_scr[0];
      for (iReg irow = 0; irow < nn_loc; irow++){
         iat_I[shift+irow] = istart_I;
         istart_scr = iend_scr;
         iend_scr = iat_scr[irow+1];
         iReg len_I = iend_scr - istart_scr;
         for (iReg k = 0; k < len_I; k++){
            ja_I[istart_I+k] = ja_scr[istart_scr+k];
            coef_I[istart_I+k] = coef_scr[istart_scr+k];
         }
         istart_I += len_I;
      }


      if ( myid == nthreads-1 ) {
         for (iReg i = 0; i < nn_L; i++){
            iat_I[0] = 0;
         }
         for (iReg i = nn_L+nn_C; i <= nn_S; i++){
            iat_I[i] = ridv_i[nthreads];
         }
      }

      exit_omp: ;
      #pragma omp atomic update
      ierr += ierr_L;

   }


   //if (params.verbosity >= VLEV_MEDIUM){
   //   iReg tot_count = 0;
   //   for (iReg i = 0; i <= params.dist_max+2; i++) tot_count += dist_count[i];
   //   rExt r_tot_count = 100.0 / static_cast<rExt>(tot_count);
   //   for (iReg i = 0; i < params.dist_max; i++)
   //      cout << "# of nodes interp at dist " << setw(4) << i+1 << ":  " <<
   //      setw(12) << dist_count[i] << " | " << setprecision(2) << setw(6) <<
   //      static_cast<rExt>(dist_count[i]) * r_tot_count << "%" << endl;
   //   cout << "# of nodes with high error:      " <<
   //   setw(12) << dist_count[params.dist_max] << " | " << setprecision(2) <<
   //   setw(6) << static_cast<rExt>(dist_count[params.dist_max]) * r_tot_count <<
   //   "%" << endl;
   //   cout << "# of nodes with large weights:   " << setw(12) <<
   //   dist_count[params.dist_max+1] << " | " << setprecision(2) << setw(6) <<
   //   static_cast<rExt>(dist_count[params.dist_max+1]) * r_tot_count << "%" << endl;
   //   cout << "# of nodes without neighbours:   " << setw(12) <<
   //   dist_count[params.dist_max+2] << " | " << setprecision(2) << setw(6) <<
   //   static_cast<rExt>(dist_count[params.dist_max+2]) * r_tot_count << "%" << endl;
   //   cout << "---------------------------" << endl;
   //}

   return ierr;

}


