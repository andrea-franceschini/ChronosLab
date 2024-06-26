

#include <stdlib.h>

#include "mex.h"
#include "Classical_prolongation.h"


int cpt_Prolongation_Classical(const int level, const int nthreads,
                               const int *const vecstart, const int nn_A, const int nt_A,
                               const int *const iat_A, const int *const ja_A,
                               const double *const coef_A,const int *const coef_S,
                               const int *const fcnodes, const int nr_I, const int nc_I,
                               int &nt_I, int *&iat_I, int *&ja_I, double *&coef_I){

   int ierr = 0;

   ierr = Classical_prolongation(level,nthreads,vecstart,nn_A,nt_A,iat_A,ja_A,
                            coef_A,coef_S,fcnodes,nr_I,nc_I,nt_I,iat_I,ja_I,coef_I);

   return ierr;
}


void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{








   if (nrhs != 12) {
       mexErrMsgIdAndTxt("ERROR in cpt_Prolongation","12 inputs required.");
   }
   if(nlhs!=4) {
       mexErrMsgIdAndTxt("ERROR in cpt_Prolongations","4 outputs required.");
   }


   int    level      =           mxGetScalar(prhs[ 0]);
   int    np         =           mxGetScalar(prhs[ 1]);
   int    *vecstart  = (int*)    mxGetPr(    prhs[ 2]);
   int    nn_A       =           mxGetScalar(prhs[ 3]);
   int    nt_A       =           mxGetScalar(prhs[ 4]);
   int    *iat_A     = (int*)    mxGetPr(    prhs[ 5]);
   int    *ja_A      = (int*)    mxGetPr(    prhs[ 6]);
   double *coef_A    = (double*) mxGetPr(    prhs[ 7]);
   int    *coef_S    = (int*)    mxGetPr(    prhs[ 8]);
   int    *fcnodes   = (int*)    mxGetPr(    prhs[ 9]);
   int    nn_I       =           mxGetScalar(prhs[10]);
   int    nc_I       =           mxGetScalar(prhs[11]);


   int nt_I;
   int *iat_I,*ja_I;
   double *coef_I;
   int ierr = cpt_Prolongation_Classical(level,np,vecstart,nn_A,nt_A,iat_A,ja_A,coef_A,
		                         coef_S,fcnodes,nn_I,nc_I,nt_I,iat_I,ja_I,coef_I);




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


   free(iat_I);
   free(ja_I);
   free(coef_I);
}
