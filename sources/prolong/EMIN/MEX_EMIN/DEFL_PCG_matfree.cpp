#include "omp.h"
#include <math.h>
#include <algorithm>
#include <iostream>
#include <iomanip>
using namespace std;

#include "parm_EMIN.h"
#include "apply_perm.h"
#include "Orth_Q.h"
#include "ddot_par.h"
#include "dnrm2_par.h"
#include "KP_spmat.h"


int DEFL_PCG_matfree(const int np, const int nn, const int nn_C, const int nn_K,
                     const int ntv, const int *perm, const int *iperm,
                     const double *D_inv, const int *iat_A, const int *ja_A,
                     const double *coef_A, const double Tr_A, const int *iat_patt,
                     const int *ja_patt, const double *mat_Q, const double *vec_P0,
                     const double *vec_f, const int itmax, const double energy_tol,
                     int &iter, double *vec_DP){
   int ierr = 0;
   double init_energy = 0.0;
   double DE = 0, DEk, DE0;

   double *rhs = (double*) malloc( nn_K*sizeof(double) );
   double *vscr = (double*) malloc( nn_K*sizeof(double) );
   double *wscr = (double*) malloc( nn_K*sizeof(double) );
   double *res = (double*) malloc( nn_K*sizeof(double) );
   double *zvec = (double*) malloc( nn_K*sizeof(double) );
   double *pvec = (double*) malloc( nn_K*sizeof(double) );
   double *QKpvec = (double*) malloc( nn_K*sizeof(double) );
   double *ridv = (double*) malloc( np*sizeof(double) );
   double *v_ntv = (double*) malloc( (np*ntv)*sizeof(double) );

   if (rhs == nullptr || vscr == nullptr || wscr == nullptr || res == nullptr ||
       zvec == nullptr || pvec == nullptr || QKpvec == nullptr ||
       ridv == nullptr || v_ntv == nullptr) return ierr = 1;

   int *WNALL = nullptr;
   WNALL = (int*) malloc( (np*nn_C)*sizeof(int) );
   if (WNALL == nullptr) return ierr = 1;

   #pragma omp parallel for num_threads(np)
   for (int i = 0; i < nn_K; i++)
   {
     vec_DP[i] = 0.0;
   }

   KP_spmat(np,nn,iat_A,ja_A,coef_A,nn_C,iat_patt,ja_patt,vec_P0,vscr,WNALL);

   #if COMP_ENRG
   apply_perm(np,nn_K,perm,vec_f,wscr);

   init_energy = Tr_A + ddot_par(np,nn_K,vec_P0,vscr,ridv) -
                        2.0*ddot_par(np,nn_K,vec_P0,wscr,ridv);
   //cout << setprecision(6) << scientific;
   //cout << "Initial Energy:  " << init_energy << endl;
   //cout << "Trace of A:      " << Tr_A << endl;
   #endif

   #pragma omp parallel for num_threads(np)
   for (int i = 0; i < nn_K; i++){
      vscr[i] -= vec_f[perm[i]];
   }

   Orth_Q(np,nn,nn_K,ntv,iat_patt,ja_patt,mat_Q,vscr,v_ntv,rhs);

   #pragma omp parallel for num_threads(np)
   for (int i = 0; i < nn_K; i++) res[i] = rhs[i];
   double bnorm = dnrm2_par(np,nn_K,rhs,ridv);

   iter = 0;
   bool exit_test = (itmax <= 0);
   double alpha, beta, gamma, gamma_old;
   while (!exit_test){
      iter++;
      apply_perm(np,nn_K,iperm,res,vscr);

      #pragma omp parallel for num_threads(np)
      for (int i = 0; i < nn_K; i++) wscr[i] = D_inv[i]*vscr[i];

      apply_perm(np,nn_K,perm,wscr,vscr);

      Orth_Q(np,nn,nn_K,ntv,iat_patt,ja_patt,mat_Q,vscr,v_ntv,zvec);

      gamma = ddot_par(np,nn_K,res,zvec,ridv);

      if (iter == 1){
         #pragma omp parallel for num_threads(np)
         for (int i = 0; i < nn_K; i++) pvec[i] = zvec[i];
      } else {
         beta = gamma / gamma_old;
         #pragma omp parallel for num_threads(np)
         for (int i = 0; i < nn_K; i++) pvec[i] = zvec[i] + beta*pvec[i];
      }

      gamma_old = gamma;

      KP_spmat(np,nn,iat_A,ja_A,coef_A,nn_C,iat_patt,ja_patt,pvec,vscr,WNALL);

      Orth_Q(np,nn,nn_K,ntv,iat_patt,ja_patt,mat_Q,vscr,v_ntv,QKpvec);

      alpha = ddot_par(np,nn_K,QKpvec,pvec,ridv);

      DEk = gamma*gamma / alpha;
      alpha = gamma / alpha;

      #pragma omp parallel for num_threads(np)
      for (int i = 0; i < nn_K; i++){
         vec_DP[i] += alpha*pvec[i];
         res[i] -= alpha*QKpvec[i];
      }

      if (iter==1){
         DE0 = DEk;
         //printf("%4s %15s %15s\n","iter","Energy","DE");
      }
      double const dDE = DEk/DE0;
      DE -= DEk;
      //printf("%4d %15.6e %15.6e\n",iter,init_energy+DE,dDE);
      exit_test = (iter == itmax) || (dDE < energy_tol);
   }

   free(v_ntv);
   free(rhs);
   free(res);
   free(vscr);
   free(wscr);
   free(zvec);
   free(pvec);
   free(QKpvec);
   free(ridv);
   free(WNALL);

   return ierr;

}
