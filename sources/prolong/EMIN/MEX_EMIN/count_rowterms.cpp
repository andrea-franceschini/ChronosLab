
#include "count_rowterms.h"




void count_rowterms( const int nequ, const int nterm, const int* __restrict__ ja,
                     int* __restrict__ WI ){


   for ( int i = 0; i < nequ; i++ ) {
      WI[i] = 0;
   }


   for ( int i = 0; i < nterm; i++ ) {
      WI[ja[i]]++;
   }

}
