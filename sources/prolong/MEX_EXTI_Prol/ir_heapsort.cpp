
#include "ir_heapsort.h"





void ir_heapsort(iReg* __restrict__ x1, rExt* __restrict__ x2, const iReg n){

   for (iReg node = 2; node < n+1; node ++){
      iReg i = node;
      iReg j = i/2;
      while( x1[j-1] < x1[i-1] ){
         swapi(x1[j-1],x1[i-1]);
         swapr(x2[j-1],x2[i-1]);
         i = j;
         j = i/2;
         if (i == 1) break;
      }
   }

   for (iReg i = n; i > 1; i --){
      swapi(x1[i-1],x1[0]);
      swapr(x2[i-1],x2[0]);
      iReg k = i - 1;
      iReg ik = 1;
      iReg jk = 2;
      if (k >= 3){
         if (x1[2] > x1[1]) jk = 3;
      }
      bool cont_cycle = false;
      if (jk <= k){
         if (x1[jk-1] > x1[ik-1]) cont_cycle = true;
      }
      while (cont_cycle){
         swapi(x1[jk-1],x1[ik-1]);
         swapr(x2[jk-1],x2[ik-1]);
         ik = jk;
         jk = ik*2;
         if (jk+1 <= k){
            if (x1[jk] > x1[jk-1]) jk = jk+1;
         }
         cont_cycle = false;
         if (jk <= k){
            if (x1[jk-1] > x1[ik-1]) cont_cycle = true;
         }
      }
   }

}


