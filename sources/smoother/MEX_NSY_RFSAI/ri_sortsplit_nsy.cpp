#include "ri_sortsplit_nsy.h"







void ri_sortsplit_nsy(const int n, const int ncut, double* R_vec, int* I_vec){

   int first = 1;
   int last = n;


   if (ncut < first || ncut > last) return;

   while(true){

      int mid = first;
      double val = R_vec[mid-1];

      for (int j = first+1; j < last+1; j ++){
         if (R_vec[j-1] > val){
            mid += 1;

            double tmp  = R_vec[mid-1];
            int itmp = I_vec[mid-1];
            R_vec[mid-1] = R_vec[j-1];
            I_vec[mid-1] = I_vec[j-1];
            R_vec[j-1]   = tmp;
            I_vec[j-1]   = itmp;
         }
      }


      double tmp      = R_vec[mid-1];
      R_vec[mid-1]    = R_vec[first-1];
      R_vec[first-1]  = tmp;

      int itmp      = I_vec[mid-1];
      I_vec[mid-1]   = I_vec[first-1];
      I_vec[first-1] = itmp;


      if (mid == ncut) return;
      if (mid >  ncut){
         last = mid-1;
      }else{
         first = mid+1;
      }

   }

}
