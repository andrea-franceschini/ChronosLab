#include <stdlib.h>
#include <omp.h>

#include <iostream>


#include "count_rowterms.h"
#include "mkiat_Tloc.h"
#include "mkiat_Tglo.h"
#include "mvjcol.h"





int Transp_Patt(const int nthreads, const int nrows, const int ncols, const int nterm,
                const int* __restrict__ iat, const int* __restrict__ ja,
                int* __restrict__ iat_T, int* __restrict__ ja_T, int* __restrict__ perm,
                int* __restrict__ iperm){


   int ierr = 0;


   int size_buff = (nthreads+1) * ncols;
   int *WI1_buff = (int*) malloc( size_buff*sizeof(int) );
   if (WI1_buff == nullptr) ierr = 1;
   int **WI1 = (int**) malloc( (nthreads+1)*sizeof(int*) );
   if (WI1 == nullptr) ierr = 1;
   int *WI2 = (int*) malloc( nthreads*sizeof(int) );
   if (WI2 == nullptr) ierr = 1;
   if (ierr != 0) return ierr;
   int j = 0;
   for ( int i = 0; i < nthreads+1; i++){
      WI1[i] = &(WI1_buff[j]);
      j += ncols;
   }


   #pragma omp parallel num_threads(nthreads)
   {

      int loc_nthreads = omp_get_num_threads();
      int mythid = omp_get_thread_num();
      int bsize = nrows/loc_nthreads;
      int resto = nrows%loc_nthreads;
      int firstrow, nrowsth, lastrow;
      if (mythid <= resto) {
         nrowsth = bsize+1;
         firstrow = mythid*nrowsth;
         if (mythid == resto) nrowsth--;
      } else {
         nrowsth = bsize;
         firstrow = mythid*bsize + resto;
      }
      lastrow = firstrow + nrowsth;


      int istart1 = iat[firstrow];
      int nnzth = iat[lastrow] - istart1;
      if (nnzth > 0) {
         count_rowterms(ncols,nnzth,&ja[istart1],WI1[mythid+1]);
      } else {
         for ( int i = 0; i < ncols; i++ ) {
            WI1[mythid+1][i] = 0;
         }
      }
      #pragma omp barrier


      bsize = ncols/loc_nthreads;
      resto = ncols%loc_nthreads;
      int firstrow_T, nrowsth_T;
      if (mythid <= resto) {
         nrowsth_T = bsize+1;
         firstrow_T = mythid*nrowsth_T;
         if (mythid == resto) nrowsth_T--;
      } else {
         nrowsth_T = bsize;
         firstrow_T = mythid*bsize + resto;
      }


      mkiat_Tloc(nrowsth_T,ncols,loc_nthreads,firstrow_T,WI1,&iat_T[firstrow_T],WI2[mythid]);
      #pragma omp barrier


      mkiat_Tglo(mythid,nrowsth_T,ncols,loc_nthreads,firstrow_T,WI1,WI2,&iat_T[firstrow_T]);
      #pragma omp barrier


      mvjcol(firstrow,nrowsth,ncols,nterm,nterm,&iat[firstrow],ja,ja_T,WI1[mythid],perm);


      #pragma omp barrier

      #pragma omp for
      for (int i = 0; i < nterm; i++) iperm[perm[i]] = i;

   }


   free(WI1_buff);
   free(WI1);
   free(WI2);

   return ierr;

}
