

#include <stdlib.h>
#include "omp.h"
using namespace std;

#include "MIS.h"
#include "ProlStripe_HybridClas.h"

int HybridClas_prolongation(const int level, const int nthreads,
                            const int *vecstart, const int nn_A, const int nt_A,
                            const int *const iat_A, const int *const ja_A,
                            const double *const coef_A, const int *const coef_S,
                            const int *const fcnodes, const int nr_I, const int nc_I,
                            int &nt_I, int *&iat_I, int *&ja_I, double *&coef_I){


   int *ridv_i;


   int ierr = 0;


   ridv_i = new int [nthreads+1](); if (ridv_i == NULL) exit(ierr += 1);



   #pragma omp parallel num_threads(nthreads)

   {

      int istart_I, shift, istart_scr, iend_scr;
      int *iat_scr, *ja_scr;
      double *coef_scr;


      int ierr_L = 0;


      int myid = omp_get_thread_num();


      int firstrow = vecstart[myid];
      int lastrow  = vecstart[myid + 1];
      int nn_loc = lastrow - firstrow;
      int nt_Imax = 100*nn_loc;
      int nt_I_loc;


      iat_scr = new int [nn_loc+1]();
      if (iat_scr == NULL) {ierr_L = 1; goto exit_pragma;}
      ja_scr = new int [nt_Imax]();
      if (ja_scr == NULL) {ierr_L = 1; goto exit_pragma;}
      coef_scr = new double [nt_Imax]();
      if (coef_scr == NULL) {ierr_L = 1; goto exit_pragma;}

      ierr_L = ProlStripe_HybridClas(firstrow,lastrow,nn_loc,nn_A,nt_A,nt_Imax,
                                     iat_A,ja_A,coef_A,coef_S,fcnodes,nt_I_loc,
                                     iat_scr,ja_scr,coef_scr);
      #pragma omp atomic update
      ierr += ierr_L;


      ridv_i[myid+1] = nt_I_loc;
      #pragma omp barrier

      if (ierr != 0) goto exit_pragma;


      #pragma omp single
      {
         ridv_i[0] = 0;
         for (int ip = 0; ip < nthreads; ip++){
             ridv_i[ip+1] = ridv_i[ip] + ridv_i[ip+1];
	     }
	     nt_I = ridv_i[nthreads];
         iat_I = new int [nn_A+1](); if (iat_I == NULL) ierr += 1;
         ja_I = new int [nt_I](); if (ja_I == NULL) ierr += 1;
         coef_I = new double [nt_I](); if (coef_I == NULL) ierr += 1;
      }
      if (ierr != 0) goto exit_pragma;


      istart_I = ridv_i[myid];
      shift = firstrow;
      iend_scr = iat_scr[0];
      for (int irow = 0; irow < nn_loc; irow++){
         iat_I[shift+irow] = istart_I;
         istart_scr = iend_scr;
         iend_scr = iat_scr[irow+1];
         int len_I = iend_scr - istart_scr;
	 for (int k = 0; k < len_I; k++){
            ja_I[istart_I+k] = ja_scr[istart_scr+k];
            coef_I[istart_I+k] = coef_scr[istart_scr+k];
	 }
         istart_I += len_I;
      }
      if (myid == nthreads-1) iat_I[nn_A] = ridv_i[nthreads];


      delete [] iat_scr;
      delete [] ja_scr;
      delete [] coef_scr;


      exit_pragma: ;


      #pragma omp atomic
      ierr += ierr_L;

   }

   delete [] ridv_i;

   return ierr;
}
