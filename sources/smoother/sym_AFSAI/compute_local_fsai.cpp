#include "mex.h"
#include "matrix.h"

#include <vector>
#include <cmath>
#include <algorithm>
#include <lapacke.h>
#include <omp.h>

typedef int iReg;
typedef int64_t iExt;
typedef double rExt;
//typedef int iBlas;
//typedef int iLapack;

rExt cpt_ddot(iReg n, rExt *x1, rExt *x2){
   rExt ddot = 0.;
   for ( iReg i = 0; i < n; i++ ){
      ddot += x1[i] * x2[i];
   }
   return ddot;
}

rExt cpt_dnrm2(iReg n, rExt *x1){
   rExt dnrm2 = cpt_ddot(n,x1,x1);
   dnrm2 = sqrt(dnrm2);
   return dnrm2;
}

void copy_ja(iExt n, iReg *x1, iReg *x2){
   for ( iExt i = 0; i < n; i++ ){
      x2[i] = x1[i];
   }
}

void copy_coef(iExt n, rExt *x1, rExt *x2){
   for ( iExt i = 0; i < n; i++ ){
      x2[i] = x1[i];
   }
}

void ri_sortsplit(iReg n, iReg ncut, rExt * const R_vec, iReg * const I_vec){

iReg first = 1;
iReg last = n;

if (ncut < first || ncut > last) return;

while(true){

   iReg mid = first;
   rExt absval = std::abs( R_vec[mid-1] );

   for (iReg j = first+1; j < last+1; j ++){

      if (std::abs(R_vec[j-1]) > absval){
         mid += 1;

         rExt tmp  = R_vec[mid-1];
         iReg itmp = I_vec[mid-1];
         R_vec[mid-1] = R_vec[j-1];
         I_vec[mid-1] = I_vec[j-1];
         R_vec[j-1]   = tmp;
         I_vec[j-1]   = itmp;
      }
   }

   rExt tmp        = R_vec[mid-1];
   R_vec[mid-1]    = R_vec[first-1];
   R_vec[first-1]  = tmp;

   iReg itmp      = I_vec[mid-1];
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

void swapi(iReg &i1, iReg &i2){
iReg tmp = i1;
i1 = i2;
i2 = tmp;
}

void iheapsort(iReg * const x1, iReg n){

for (iReg node = 2; node < n+1; node ++){
   iReg i = node;
   iReg j = i/2;
   while( x1[j-1] < x1[i-1] ){
      swapi(x1[j-1],x1[i-1]);
      i = j;
      j = i/2;
      if (i == 1) break;
   }
}

for (iReg i = n; i > 1; i --){
   swapi(x1[i-1],x1[0]);
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

void gather_fullsys(bool &nulrhs, iReg irow, iReg mrow, iReg jendbloc, iReg nequ,
                    iExt nterm, iReg mmax, const iExt * const iat, const iReg * const ja,
                    const iReg * const vecinc, const rExt * const mat_A,
                    rExt * const full_A, rExt * const rhs){

for (iReg i = 1; i < mrow+1; i++){

   iReg ii     = i;
   iReg row    = vecinc[i-1];
   iExt jj     = iat[row-1];
   iExt endrow = iat[row];

   while (ii <= mrow){

      while (ja[jj-1] < vecinc[ii-1]){
         jj += 1;
         if (jj == endrow) {

            for (iReg k = ii; k < mrow+1; k++){
               full_A[(i-1)*mmax+k-1] = 0.;
            }
            goto exit_column_loop;
         }
      }

      if (vecinc[ii-1] == ja[jj-1]){

         full_A[(i-1)*mmax+ii-1] = mat_A[jj-1];
         ii += 1;
      }else{

         full_A[(i-1)*mmax+ii-1] = 0.;
         ii += 1;
      }

   }
   exit_column_loop: ;

}

iReg ii = 1;
iExt jj = iat[irow-1];
nulrhs = ja[jj-1] >= jendbloc;
if (nulrhs) return;
while (ja[jj-1] < jendbloc){
   if (vecinc[ii-1] > ja[jj-1]){
      jj += 1;
   }else if (vecinc[ii-1] == ja[jj-1]){
      rhs[ii-1] = - mat_A[jj-1];
      jj += 1;
      ii += 1;
      if (ii > mrow) return;
   }else{
      rhs[ii-1] = 0.;
      ii += 1;
      if (ii > mrow) return;
   }
}
for (iReg k = ii; k < mrow+1; k++){
   rhs[k-1] = 0.;
}

}

void kap_grad(iReg irow, iReg nequ, iExt nterm, iReg mmax, iReg &mrow, iReg jendbloc,
              iReg lfil, const iExt * const iat, const iReg * const ja,
              const rExt * const mat_A, const rExt * const rhs, iReg * const IWN,
              iReg * const JWN, rExt * const WR){

iReg ind_WR = 0;

iExt j = iat[irow-1];
iReg jjcol = ja[j-1];
while (jjcol < jendbloc){

   if (JWN[jjcol-1] >= 0) {
      ind_WR += 1;
      JWN[jjcol-1] = ind_WR;
      IWN[mrow+ind_WR-1] = ja[j-1];
      WR[ind_WR-1] = mat_A[j-1];
   }
   j += 1;
   jjcol = ja[j-1];
}

for (iReg i = 1; i < mrow+1; i ++){
   iReg iirow = IWN[i-1];
   for (iExt j = iat[iirow-1]; j < iat[iirow]; j ++){
      jjcol = ja[j-1];
      if (jjcol < jendbloc){
         iReg ind = JWN[jjcol-1];
         if (ind == 0){

            ind_WR += 1;
            JWN[jjcol-1] = ind_WR;
            IWN[mrow+ind_WR-1] = ja[j-1];
            WR[ind_WR-1] = rhs[i-1]*mat_A[j-1];
         }else if (ind > 0){

            WR[ind-1] += rhs[i-1]*mat_A[j-1];
         }
      }
   }
}

iReg ncut = std::min(lfil,ind_WR);
ri_sortsplit(ind_WR,ncut,WR,&IWN[mrow]);
for (iReg i = mrow+1; i < mrow+ncut+1; i ++){
   JWN[IWN[i-1]-1] = -1;
}

for (iReg i = mrow+ncut+1; i < mrow+ind_WR+1; i ++){
   JWN[IWN[i-1]-1] = 0;
}

mrow += ncut;
iheapsort(IWN,mrow);

}

void cpt_afsai_coef(iReg chunk_size, iReg n_step, iReg step_size, rExt tau, rExt eps,
                    iReg shift, iReg nrows, iReg nequ, iExt nterm, iExt &nterm_G,
                    const iExt * const iat, const iReg * const ja, iExt * const istart_G,
                    iExt * const istop_G, iReg * const ja_G, const rExt * const coef_A,
                    rExt * const coef_G){
iReg mrow_min = 5;
iReg mmax = n_step * step_size;
std::vector<iReg> IWN(nequ,0);
std::vector<iReg> JWN(nequ,0);
std::vector<rExt> WR(nequ,0);
std::vector<rExt> full_A(mmax*mmax,0);
std::vector<rExt> rhs(mmax+1,0);
std::vector<rExt> rhs_sav(mmax+1,0);
nterm_G = 0;
#pragma omp for nowait schedule(dynamic,chunk_size)

for( iReg irow = 1; irow < nrows+1; irow++){
   iReg mrow = 0;
   iReg irow_glo = irow + shift;
   iReg istep = 0;
   rExt DKap_old = 0.;
   bool Refine = n_step >= 1;
   while (Refine){

      istep += 1;

      if ((tau > 0.) && (mrow > mrow_min)){

         rhs[mrow] = 1.;

         rExt asstol = tau * cpt_dnrm2(mrow+1,rhs.data());

         iReg ind = 0;
         for (iReg i = 0; i < mrow; i ++){
            if (std::abs(rhs[i]) > asstol) {

               rhs[ind] = rhs[i];
               rhs_sav[ind] = rhs_sav[i];
               IWN[ind] = IWN[i];
               ind += 1;
            }else{

               JWN[IWN[i]-1] = -1;
            }
         }

         mrow = ind;
      }

      iReg mrow_old = mrow;
      kap_grad(irow_glo,nequ,nterm,mmax,mrow,irow_glo,step_size,iat,ja,coef_A,rhs.data(),
               IWN.data(),JWN.data(),WR.data());

      if (mrow > mrow_old){

         bool nulrhs;
         gather_fullsys(nulrhs,irow_glo,mrow,irow_glo,nequ,nterm,mmax,iat,ja,IWN.data(),
                        coef_A,full_A.data(),rhs.data());

         if (nulrhs == false){

            lapack_int info;
            lapack_int l_mrow = static_cast<lapack_int>( mrow );
            lapack_int l_mmax = static_cast<lapack_int>( mmax );
            lapack_int l_one = static_cast<lapack_int>( 1 );
            info = LAPACKE_dpotrf(LAPACK_COL_MAJOR,'L',l_mrow,full_A.data(),l_mmax);


            rhs_sav = rhs;

            info = LAPACKE_dpotrs(LAPACK_COL_MAJOR,'L',l_mrow,l_one,full_A.data(),l_mmax,rhs.data(),l_mrow);

         }

         rExt DKap_new = cpt_ddot(mrow,rhs.data(),rhs_sav.data());

         if (istep == n_step){
            Refine = false;
         }else{
            Refine = (std::abs(DKap_new-DKap_old) >= eps*DKap_old) && (DKap_new != 0.);
            DKap_old = DKap_new;
         }

      }else{

         Refine = false;

      }

   }

   iExt ind = iat[irow_glo-1];
   while (ja[ind-1] < irow_glo){
      ind += 1;
   }

   rExt scal_fac = coef_A[ind-1] - cpt_ddot(mrow,rhs_sav.data(),rhs.data());

   if (scal_fac < 0.){

      bool nulrhs;
      gather_fullsys(nulrhs,irow_glo,mrow,irow_glo,nequ,nterm,mmax,iat,ja,IWN.data(),
                     coef_A,full_A.data(),rhs.data());

      lapack_int info;
      lapack_int l_mrow = static_cast<lapack_int>( mrow );
      lapack_int l_mmax = static_cast<lapack_int>( mmax );
      lapack_int l_one = static_cast<lapack_int>( 1 );
      info = LAPACKE_dpotrf(LAPACK_COL_MAJOR,'L',l_mrow,full_A.data(),l_mmax);
      rhs_sav = rhs;

      info = LAPACKE_dpotrs(LAPACK_COL_MAJOR,'L',l_mrow,l_one,full_A.data(),l_mmax,rhs.data(),l_mrow);

      scal_fac = coef_A[ind-1] - cpt_ddot(mrow,rhs_sav.data(),rhs.data());
   }

   scal_fac = 1. / sqrt(scal_fac);

   iExt ind_G  = istart_G[irow-1];
   iExt ind_G0 = ind_G;

   for (iReg i = 1; i < mrow+1; i++){
      coef_G[ind_G-1] = scal_fac*rhs[i-1];
      ja_G[ind_G-1] = IWN[i-1];
      ind_G += 1;

      JWN[IWN[i-1]-1] = 0;
   }

   coef_G[ind_G-1] = scal_fac;
   ja_G[ind_G-1]   = irow_glo;

   istop_G[irow-1] = ind_G;

   nterm_G += (ind_G - ind_G0) + 1;

}
}

void compute_local_fsai(iReg nthread, iReg n_step, iReg step_size, rExt tau, rExt eps,
                        iReg nrows, iReg nrows_M, iExt nterm_M, const iExt * const iat_M,
                        const iReg * const ja_M, const rExt * const coef_M, iExt *nterm_G,
                        iExt *iat_G, iReg *ja_G, rExt *coef_G){

iReg chunk_size = nrows / (20*nthread); chunk_size = std::max(1,chunk_size);

iExt kmax    = 1 + (iExt) n_step * (iExt) step_size;
iExt nzmax_G = (iExt) nrows * kmax;

std::vector<iExt> nt_slice; nt_slice.resize(nthread);
std::vector<iExt> istart_scr; istart_scr.resize(nrows);
std::vector<iExt> istop_scr; istop_scr.resize(nrows);
std::vector<iReg> ja_scr; ja_scr.resize(nzmax_G);
std::vector<rExt> coef_scr; coef_scr.resize(nzmax_G);

iExt ind = 1;
for ( iReg irow = 0; irow < nrows; irow++ ){
   istart_scr[irow] = ind;
   ind += kmax;
}

nterm_G[0] = 0;
iReg shift = nrows_M - nrows;

#pragma omp parallel num_threads(nthread)
{

   iReg myid = omp_get_thread_num();

   iExt loc_nt_G = 0;
   double time_1 = omp_get_wtime();
   cpt_afsai_coef(chunk_size,n_step,step_size,tau,eps,shift,nrows,nrows_M,
                  nterm_M,loc_nt_G,iat_M,ja_M,istart_scr.data(),istop_scr.data(),
                  ja_scr.data(),coef_M,coef_scr.data());
   double time_2 = omp_get_wtime();

   #pragma omp atomic
   nterm_G[0] += loc_nt_G;

}

#pragma omp parallel num_threads(nthread)
{
   iReg myid = omp_get_thread_num();

   iReg my_nrows    = nrows / nthread;
   iReg my_firstrow = myid * my_nrows + 1;
   iReg my_lastrow  = my_firstrow + my_nrows - 1;
   if ( myid + 1 == nthread ){
      iReg resto = nrows%nthread;
      my_lastrow += resto;
      my_nrows   += resto;
   }
   iExt kcount = my_nrows;
   for ( iReg irow = my_firstrow-1; irow < my_lastrow; irow++ ) {
      kcount += istop_scr[irow] - istart_scr[irow];
   }
   nt_slice[myid] = kcount;

   #pragma omp barrier

   iExt ind_G  = 1;
   for ( iReg id = 0; id < myid; id++ ) {
      ind_G += nt_slice[id];
   }
   for ( iReg irow = my_firstrow-1; irow < my_lastrow; irow++ ) {

      iat_G[irow] = ind_G;
      iExt istart = istart_scr[irow];
      iExt istop  = istop_scr[irow];
      iExt nadd   = istop-istart;

      iExt j = ind_G-1;
      for ( iExt i = istart-1; i < istop; i++ ) {
         ja_G[j]   = ja_scr[i];
         coef_G[j] = coef_scr[i];
         j += 1;
      }

      ind_G += nadd+1;
   }

   if (myid + 1 == nthread){ iat_G[my_lastrow] = ind_G;}

}

}

void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[]){

   //mexPrintf("- get input scalars\n");

   iReg nthread   = mxGetScalar(prhs[0]);
   iReg n_step    = mxGetScalar(prhs[1]);
   iReg step_size = mxGetScalar(prhs[2]);
   rExt tau       = mxGetScalar(prhs[3]);
   rExt eps       = mxGetScalar(prhs[4]);
   iReg nrows     = mxGetScalar(prhs[5]);
   iReg nrows_M   = mxGetScalar(prhs[6]);
   iExt nterm_M   = mxGetScalar(prhs[7]);

   //mexPrintf("- get input arrays\n");
   iExt *iat_M  = (iExt*) mxGetData(prhs[8]);
   iReg *ja_M   = (iReg*) mxGetData(prhs[9]);
   rExt *coef_M = (rExt*) mxGetData(prhs[10]);

   //mexPrintf("- allocate output scalars\n");

   plhs[0] = mxCreateNumericMatrix((mwSize)1,(mwSize)1,mxINT64_CLASS,mxREAL);
   iExt *nterm_G = (iExt*) mxGetData(plhs[0]);

   //mexPrintf("- allocate work arrays\n");
   iExt kmax    = 1 + (iExt) n_step * (iExt) step_size;
   iExt nzmax_G = (iExt) nrows * kmax;

   plhs[1] = mxCreateNumericMatrix((mwSize)1,(mwSize)nrows+1,mxINT64_CLASS,mxREAL);
   iExt *iat_G = (iExt*) mxGetData(plhs[1]);

   mxArray *ja_W = mxCreateNumericMatrix((mwSize)1,(mwSize)nzmax_G,mxINT32_CLASS,mxREAL);
   iReg *ja_G = (iReg*) mxGetData(ja_W);

   mxArray *coef_W = mxCreateNumericMatrix((mwSize)1,(mwSize)nzmax_G,mxDOUBLE_CLASS,mxREAL);
   rExt *coef_G = (rExt*) mxGetData(coef_W);

   //mexPrintf("- compute G terms\n");
   compute_local_fsai(nthread,n_step,step_size,tau,eps,nrows,nrows_M,nterm_M,iat_M,
                      ja_M,coef_M,nterm_G,iat_G,ja_G,coef_G);

   //mexPrintf("- resize output arrays\n");

   mxArray *ja_R = mxCreateNumericMatrix((mwSize)1,(mwSize)nterm_G[0],mxINT32_CLASS,mxREAL);
   iReg *ptr_ja_R = (iReg*) mxGetData(ja_R);
   copy_ja(nterm_G[0],ja_G,ptr_ja_R);
   mxDestroyArray(ja_W);

   mxArray *coef_R = mxCreateNumericMatrix((mwSize)1,(mwSize)nterm_G[0],mxDOUBLE_CLASS,mxREAL);
   rExt *ptr_coef_R = (rExt*) mxGetData(coef_R);
   copy_coef(nterm_G[0],coef_G,ptr_coef_R);
   mxDestroyArray(coef_W);

   plhs[2] = ja_R;
   plhs[3] = coef_R;

}
