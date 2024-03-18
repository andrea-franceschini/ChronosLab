
#include "check_neigh_cond.h"



void check_neigh_cond(const rExt maxcond, const iReg inod, const iReg ntv, const iReg nn_S,
                      const iExt *const iat_S, const iReg *const ja_S,
                      const rExt *const *const TV,
                      rExt **TVcomp, iReg *neigh, iReg *WI, rExt &local_cond){


   iReg n_neigh = 1;
   neigh[0] = inod;
   WI[inod] = 1;
   try {
      add_new_neighs(iat_S,ja_S,0,1,nn_S,n_neigh,neigh,WI);
   } catch (linsol_error) {
      throw linsol_error ("check_neigh_cond","in adding new neighbours");
   }


   for (iReg i = 1; i < n_neigh; i++){
      for (iReg j = 0; j < ntv; j++) TVcomp[i][j] = TV[neigh[i]][j];
   }


   iReg rank;
   get_cond(maxcond,n_neigh-1,ntv,&(TVcomp[1]),rank,local_cond);


   for (iReg i = 0; i < n_neigh; i++) WI[neigh[i]] = 0;

}
