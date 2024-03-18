#include "mk_HouHolVec.h"




void mk_HouHolVec(const iReg n, const rExt *const v, rExt *w){

   rExt norm,norm_2,norm_w;


   norm_2 = 0.0;
   for (iReg i = 0; i < n; i++){
      w[i] = v[i];
      norm_2 += w[i]*w[i];
   }
   norm = sqrt(norm_2);

   if (w[0] == 0.0){
      w[0] = norm;
      norm_w = sqrt(2.0*norm_2);
   } else {
      if (w[0]>0.0){
         norm_w = sqrt(2.0*(w[0]*norm+norm_2));
         w[0] += norm;
      } else {
         norm_w = sqrt(2.0*(-w[0]*norm+norm_2));
         w[0] -= norm;
      }
   }


   rExt w_inv = 1.0 / norm_w;
   for (iReg i = 0; i < n; i++) w[i] *= w_inv;

}
