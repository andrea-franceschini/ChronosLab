
#include "KapGrad_NSY.h"




void KapGrad_NSY(const int istep, const int irow, int& mrow, const int jendbloc,
                 const int lfil, const int* __restrict__ iat, const int* __restrict__ ja,
                 const double* __restrict__ coef, const double* __restrict__ coef_T,
                 const double* __restrict__ sol_L, const double* __restrict__ sol_U,
                 int* __restrict__ IWN, int* __restrict__ JWN,
                 double* __restrict__ WR_L, double* __restrict__ WR_U){

int ind_WR = 0;


int j = iat[irow];
int jjcol = ja[j];
while (jjcol < jendbloc){

   if (JWN[jjcol] >= 0) {
      ind_WR ++;
      JWN[jjcol] = ind_WR;
      IWN[mrow+ind_WR-1] = ja[j];
      WR_L[ind_WR-1] = coef[j];
      WR_U[ind_WR-1] = coef_T[j];
   }
   j ++;
   jjcol = ja[j];
}


for (int i = 0; i < mrow; i ++){
   int iirow = IWN[i];
   for (int j = iat[iirow]; j < iat[iirow+1]; j ++){
      jjcol = ja[j];
      if (jjcol < jendbloc){
         int ind = JWN[jjcol];
         if (ind == 0){

            ind_WR++;
            JWN[jjcol] = ind_WR;
            IWN[mrow+ind_WR-1] = ja[j];
            WR_L[ind_WR-1] = sol_L[i]*coef[j];
            WR_U[ind_WR-1] = sol_U[i]*coef_T[j];
	 }else if (ind > 0){

            WR_L[ind-1] += sol_L[i]*coef[j];
            WR_U[ind-1] += sol_U[i]*coef_T[j];
         }
      }
   }
}


for (int i = 0; i < ind_WR; i ++){
   WR_L[i] = fabs(WR_L[i]) + fabs(WR_U[i]);
}


if (DEBUG){
   fprintf(dbfile,"IWN: ");
   for (int i = 0; i < ind_WR; i ++){
      fprintf(dbfile," %6d",IWN[mrow+i]);
   }
   fprintf(dbfile,"\n");
   fprintf(dbfile,"GRAD: ");
   for (int i = 0; i < ind_WR; i ++){
      fprintf(dbfile," %15.6e",0.5*WR_L[i]);
   }
   fprintf(dbfile,"\n");
}



int ncut = min(lfil,ind_WR);
ri_sortsplit_nsy(ind_WR,ncut,WR_L,&IWN[mrow]);


for (int i = mrow; i < mrow+ncut; i ++){
   JWN[IWN[i]] = -istep;
}


for (int i = mrow+ncut; i < mrow+ind_WR; i ++){
   JWN[IWN[i]] = 0;
}


mrow += ncut;
heapsort(IWN,mrow);

}
