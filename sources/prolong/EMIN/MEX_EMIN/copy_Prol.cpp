#include "omp.h"





void copy_Prol(const int np, const int nn, const int *fcnode, const int *iat_in,
               const int *ja_in, const double *coef_in, const int *iat_out,
               const int *ja_out, double *coef_out){

   #pragma omp parallel num_threads(np)
   {

      int mythid = omp_get_thread_num();
      int bsize = nn/np;
      int resto = nn%np;
      int firstrow, nrowth, lastrow;
      if (mythid <= resto) {
         nrowth = bsize+1;
         firstrow = mythid*nrowth;
         if (mythid == resto) nrowth--;
      } else {
         nrowth = bsize;
         firstrow = mythid*bsize + resto;
      }
      lastrow = firstrow + nrowth;


      int iend = iat_in[firstrow];
      for (int irow = firstrow; irow < lastrow; irow++){
         int ii = iend;
         iend = iat_in[irow+1];
         if (fcnode[irow] < 0){
            int jj = iat_out[irow];
            while (ii < iend){


               while (ja_out[jj] < ja_in[ii]) jj++;

               if (ja_out[jj] == ja_in[ii]){
                  coef_out[jj] = coef_in[ii];
                  jj++;
               }
               ii++;

            }
         }
      }

   }

}
