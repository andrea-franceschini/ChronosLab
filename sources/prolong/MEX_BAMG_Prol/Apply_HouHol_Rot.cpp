
#include "Apply_HouHol_Rot.h"




void Apply_HouHol_Rot(iReg n, const rExt *const w, rExt *v){

   rExt fac = inl_ddot(n,w,1,v,1);
   inl_daxpy(n,-2.*fac,w,1,v,1);

}
