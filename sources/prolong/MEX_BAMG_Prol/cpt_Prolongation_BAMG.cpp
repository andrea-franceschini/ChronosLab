

#include <stdlib.h>

#include <iostream>

#include "mex.h"
#include "BAMG.h"
#include "BAMG_params.h"


int cpt_Prolongation_BAMG(const int level, const BAMG_params params, const int nthreads,
                          const int nn_S, const int nt_S, const int *const iat_S,
                          const int *const ja_S, const int *const coef_S, const int ntv,
                          const int *const fcnodes, const double *const *const TV,
                          const int nc_I, int &nt_I, std::vector<int> &vec_iat_I,
                          std::vector<int> &vec_ja_I, std::vector<double> &vec_coef_I,
                          std::vector<int> &vec_c_mark){


   iReg nn_L = 0;
   iReg nn_C = nn_S;
   int ierr = BAMG(params,nthreads,nn_L,nn_C,nn_S,iat_S,ja_S,ntv,fcnodes,TV,
                   nt_I,vec_iat_I,vec_ja_I,vec_coef_I,vec_c_mark);

   return ierr;
}


void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{








   if (nrhs != 20) {
       mexErrMsgIdAndTxt("ERROR in cpt_Prolongation","22 inputs required.");
   }
   if(nlhs!=5) {
       mexErrMsgIdAndTxt("ERROR in cpt_Prolongations","5 outputs required.");
   }


   int    level      =           mxGetScalar(prhs[ 0]);
   int    np         =           mxGetScalar(prhs[ 1]);
   int    itmax_vol  =           mxGetScalar(prhs[ 2]);
   int    dist_min   =           mxGetScalar(prhs[ 3]);
   int    dist_max   =           mxGetScalar(prhs[ 4]);
   int    mmax       =           mxGetScalar(prhs[ 5]);
   double maxcond    =           mxGetScalar(prhs[ 6]);
   double maxrownrm  =           mxGetScalar(prhs[ 7]);
   double tol_vol    =           mxGetScalar(prhs[ 8]);
   double eps_prol   =           mxGetScalar(prhs[ 9]);
   int    nn_S       =           mxGetScalar(prhs[10]);
   int    nt_S       =           mxGetScalar(prhs[11]);
   int    *iat_S     = (int*)    mxGetPr(    prhs[12]);
   int    *ja_S      = (int*)    mxGetPr(    prhs[13]);
   int    *coef_S    = (int*)    mxGetPr(    prhs[14]);
   int    ntv        =           mxGetScalar(prhs[15]);
   int    *fcnodes   = (int*)    mxGetPr(    prhs[16]);
   double *TV        = (double*) mxGetPr(    prhs[17]);
   int    nn_I       =           mxGetScalar(prhs[18]);
   int    nc_I       =           mxGetScalar(prhs[19]);


   BAMG_params params;
   params.verbosity = 0;
   params.itmax_vol = itmax_vol;
   params.dist_min = dist_min;
   params.dist_max = dist_max;
   params.mmax = mmax;
   params.maxcond = maxcond;
   params.maxrownrm = maxrownrm;
   params.tol_vol = tol_vol;
   params.eps = eps_prol;


   double **TV_2D = (double**) malloc(nn_S*sizeof(double*));
   if (TV_2D == NULL) mexErrMsgIdAndTxt("ERROR in cpt_Prolongation","TV_2D allocation");
   int kk = 0;
   for (int i = 0; i < nn_S; i++){
      TV_2D[i] = (double*) malloc(ntv*sizeof(double));
      if (TV_2D[i] == NULL) mexErrMsgIdAndTxt("ERROR in cpt_Prolongation","TV_2D allocation");
      for (int j = 0; j < ntv; j++) {
         TV_2D[i][j] = TV[kk];
         kk++;
      }
   }

   int nt_I;
   std::vector<int> vec_iat_I,vec_ja_I;
   std::vector<double> vec_coef_I;
   std::vector<int> vec_c_mark;
   int ierr = cpt_Prolongation_BAMG(level,params,np,nn_S,nt_S,iat_S,ja_S,coef_S,
                                    ntv,fcnodes,TV_2D,nc_I,nt_I,vec_iat_I,vec_ja_I,
                                    vec_coef_I,vec_c_mark);




   int    *iat_I  = vec_iat_I.data();
   int    *ja_I   = vec_ja_I.data();
   double *coef_I = vec_coef_I.data();
   int    *c_mark = vec_c_mark.data();


   plhs[0] = mxCreateDoubleMatrix(1,1,mxREAL);
   double *out_nt_I = mxGetPr(plhs[0]);
   *out_nt_I = (double) nt_I;


   plhs[1] = mxCreateDoubleMatrix((mwSize) nn_I+1,1,mxREAL);
   double *out_iat_I = mxGetPr(plhs[1]);
   for (int i = 0; i < nn_I+1; i++) out_iat_I[i] = (double) (iat_I[i]+1);


   plhs[2] = mxCreateDoubleMatrix((mwSize) nt_I,1,mxREAL);
   double *out_ja_I = mxGetPr(plhs[2]);
   for (int i = 0; i < nt_I; i++) out_ja_I[i] = (double) (ja_I[i]+1);


   plhs[3] = mxCreateDoubleMatrix((mwSize) nt_I,1,mxREAL);
   double *out_coef_I = mxGetPr(plhs[3]);
   for (int i = 0; i < nt_I; i++) out_coef_I[i] = coef_I[i];


   plhs[4] = mxCreateDoubleMatrix((mwSize) nn_I,1,mxREAL);
   double *out_c_mark = mxGetPr(plhs[4]);
   for (int i = 0; i < nn_I; i++) out_c_mark[i] = (double) (c_mark[i]);


   for (int i = 0; i < nn_S; i++) free(TV_2D[i]); free(TV_2D);

}
