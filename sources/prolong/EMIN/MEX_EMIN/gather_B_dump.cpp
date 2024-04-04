#include <stdlib.h>
#include <omp.h>
#include <limits>

#include "parm_EMIN.h"

int gather_B_dump(const int np, const int nn, const int nn_C, const int ntv,
                  const int *fcnode, const int *iat_patt, const int *ja_patt,
                  const double *const *TV, double *&mat_B){

   int ierr = 0;

   int nrows_B = iat_patt[nn];
   int nterm_B = ntv*nrows_B;
   mat_B = (double*) malloc( nterm_B*sizeof(double) );
   if (mat_B == nullptr) return ierr = 3;

   int *c2glo = (int*) malloc( nn_C*sizeof(int) );
   if (c2glo == nullptr)  return ierr = 1;

   #pragma omp parallel num_threads(np)
   {

      int loc_nthreads = omp_get_num_threads();
      int mythid = omp_get_thread_num();
      int bsize = nn/loc_nthreads;
      int resto = nn%loc_nthreads;
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
      int pos_B = pos_g*ntv;
      double *BB_scr = &(mat_B[pos_B]);

      #pragma omp for
      for (int i = 0; i < nn; i++){
         int k = fcnode[i];
         if (k >= 0) c2glo[k] = i;
      }


      int ind_BB, nnz_BB;
      ind_BB = 0;
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
               nnz_BB = nr_BB_loc*ntv;
               ind_BB += nnz_BB;
            }
         }
      }
      exit_pragma: ;
   }

   free(c2glo);

   return ierr;
}
