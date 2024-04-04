#include <iostream>
#include <omp.h>
#include <cmath>
#include <chrono>

#if defined PRINT
#define dump true
#else
#define dump false
#endif

#include "parm_EMIN.h"
#include "wrCSRmat.h"
#include "Transp_Patt.h"
#include "load_Jacobi.h"
#include "copy_Prol.h"
#include "gather_B_QR.h"
#include "gather_B_dump.h"
#include "gather_f.h"
#include "cpt_Trace_Acc.h"
#include "print_Q.h"
#include "Prol_add_Cnodes.h"
#include "DEFL_PCG_matfree.h"

int EMIN_matfree(const int np, const int itmax, const double en_tol, const double condmax,
                 const double maxwgt, const int min_lfil, const int max_lfil, const int D_lfil,
                 const int nn, const int nn_C, const int ntv, const int nt_A, const int nt_P,
                 const int nt_patt, const int *fcnode, const int *iat_A, const int *ja_A,
                 const double *coef_A, const int *iat_Pin, const int *ja_Pin,
                 const double *coef_Pin, const int *iat_patt, const int *ja_patt,
                 const double *const *TV, int *&iat_Pout, int *&ja_Pout, double *&coef_Pout,
                 double *info){

   int ierr = 0;

   std::chrono::time_point<std::chrono::system_clock> start, end;
   std::chrono::duration<double> elaps_sec;
   printf("step 0\n");

   int *iat_Tpatt = (int*) malloc( (nn_C+1)*sizeof(int) );
   int *ja_Tpatt  = (int*) malloc( nt_patt*sizeof(int) );
   int *perm  = (int*) malloc( nt_patt*sizeof(int) );
   int *iperm  = (int*) malloc( nt_patt*sizeof(int) );
   if (iat_Tpatt == nullptr || ja_Tpatt == nullptr || perm == nullptr ||
       iperm == nullptr)
      return ierr = 1;
   ierr = Transp_Patt(np,nn,nn_C,nt_patt,iat_patt,ja_patt,iat_Tpatt,ja_Tpatt,perm,
                      iperm);
   if (ierr != 0) return ierr = 1;
   printf("step 1\n");

   int nnz_J = 0;
   double *D_inv = nullptr;
   double time_prec_K = 0.0;

   start = std::chrono::system_clock::now();
   D_inv = (double*) malloc( nt_patt*sizeof(double) );
   double *scr = (double*) malloc( nn*sizeof(double) );
   if (D_inv == nullptr || scr == nullptr) return ierr = 3;
   load_Jacobi(np,nn,nt_patt,iat_Tpatt,ja_Tpatt,iat_A,ja_A,coef_A,scr,D_inv);
   free(scr);
   nnz_J = nt_patt;
   end = std::chrono::system_clock::now();
   elaps_sec = end - start;
   time_prec_K = elaps_sec.count();
   if (DUMP_PREC){
      FILE *dfile = fopen("D_mat","w");
      for (int k = 0; k < nnz_J; k++) fprintf(dfile,"%20.11e\n",D_inv[k]);
      fflush(dfile);
      fclose(dfile);
   }
   printf("step 2\n");


   if (dump) std::cout << "---- COPY PROL ----" << std::endl << std::endl;
   double *coef_P0 = (double*) calloc( nt_patt , sizeof(double) );
   if (coef_P0 == nullptr) return ierr = 1;
   copy_Prol(np,nn,fcnode,iat_Pin,ja_Pin,coef_Pin,iat_patt,ja_patt,coef_P0);
   if (DUMP_PREC){
      FILE *origPfile = fopen("origProl.csr","w");
      wrCSRmat(origPfile,false,nn,iat_Pin,ja_Pin,coef_Pin);
      fclose(origPfile);
      FILE *extPfile = fopen("extProl.csr","w");
      wrCSRmat(extPfile,false,nn,iat_patt,ja_patt,coef_P0);
      fclose(extPfile);
   }
   printf("step 3\n");

   if (dump) std::cout << "---- gather_B_Qr ----" << std::endl << std::endl;
   start = std::chrono::system_clock::now();
   double *mat_Q = nullptr;
   double *vec_f = nullptr;
   if (DUMP_PREC){
      double *mat_B = nullptr;
      ierr = gather_B_dump(np,nn,nn_C,ntv,fcnode,iat_patt,ja_patt,TV,mat_B);
      FILE *Bfile = fopen("mat_B","w");
      print_Q(Bfile,nn,ntv,iat_patt,mat_B);
      fflush(Bfile);
      fclose(Bfile);
      free(mat_B);
   }
   ierr = gather_B_QR(np,condmax,maxwgt,nn,nn_C,ntv,fcnode,iat_patt,ja_patt,TV,
                      mat_Q,coef_P0);
   if (ierr != 0) return ierr = 4;
   int nnz_Q = ntv*iat_patt[nn];
   printf("step 4\n");

   int nrows_Q = iat_patt[nn];
   vec_f = (double*) malloc( nrows_Q*sizeof(double) );
   int *c2glo = (int*) malloc( nn_C*sizeof(int) );
   if (vec_f == nullptr || c2glo == nullptr) return ierr = 1;

   #pragma omp parallel for num_threads(np)
   for (int i = 0; i < nn; i++){
      int k = fcnode[i];
      if (k >= 0) c2glo[k] = i;
   }

   #pragma omp parallel for num_threads(np)
   for (int icol = 0; icol < nn_C; icol++){
       int istart = iat_Tpatt[icol];
       int iend = iat_Tpatt[icol+1];
       int irow = c2glo[icol];
       int istart_A = iat_A[irow];
       int iend_A = iat_A[irow+1];
       gather_f(iend-istart,&(ja_Tpatt[istart]),iend_A-istart_A,&(ja_A[istart_A]),
                &(coef_A[istart_A]),&(vec_f[istart]));
   }
   printf("step 5\n");
   end = std::chrono::system_clock::now();
   elaps_sec = end - start;
   double time_gath_B = elaps_sec.count();
   if (DUMP_PREC){
      FILE *corrPfile = fopen("corrProl.csr","w");
      wrCSRmat(corrPfile,false,nn,iat_patt,ja_patt,coef_P0);
      fclose(corrPfile);
      FILE *Qfile = fopen("mat_Q2","w");
      print_Q(Qfile,nn,ntv,iat_patt,mat_Q);
      fflush(Qfile);
      fclose(Qfile);
      FILE *ffile = fopen("vec_f","w");
      for (int i = 0; i < nt_patt; i++) fprintf(ffile,"%15.6e\n",vec_f[i]);
      fflush(ffile);
      fclose(ffile);
   }
   double Tr_A = 0;

   #if COMP_ENRG
   Tr_A = cpt_Trace_Acc(np,nn,fcnode,iat_A,ja_A,coef_A);
   #endif

   int iter;

   if (dump) std::cout << "---- DEFL_PCG_matfree ----" << std::endl << std::endl;
   start = std::chrono::system_clock::now();

   double *DP = (double*) malloc( nt_patt*sizeof(double) );
   if (DP == nullptr) return ierr = 1;

   ierr = DEFL_PCG_matfree(np,nn,nn_C,nt_patt,ntv,perm,iperm,D_inv,iat_A,ja_A,coef_A,Tr_A,
                           iat_patt,ja_patt,mat_Q,coef_P0,vec_f,itmax,en_tol,iter,DP);
   if (ierr != 0) return ierr = 5;
   end = std::chrono::system_clock::now();
   elaps_sec = end - start;
   double time_PCG = elaps_sec.count();

   #pragma omp parallel for num_threads(np)
   for (int i = 0; i < nt_patt; i++) coef_P0[i] -= DP[i];

   free(DP);
   free(vec_f);
   free(mat_Q);
   free(D_inv);
   free(perm);
   free(iperm);
   free(iat_Tpatt);
   free(ja_Tpatt);

   ierr = Prol_add_Cnodes(np,nn,nn_C,fcnode,iat_patt,ja_patt,coef_P0,iat_Pout,ja_Pout,
                          coef_Pout);
   if (ierr != 0) ierr = 6;
   free(coef_P0);

   info[0]  = 0.0;
   info[1]  = time_prec_K;
   info[2]  = time_gath_B;
   info[3]  = time_PCG;
   info[6]  = 0.0;
   info[7]  = static_cast<double>(iter);
   info[8]  = 0.0;
   info[9]  = 0.0;
   info[10] = static_cast<double>(nnz_J);
   info[11] = static_cast<double>(nnz_Q);

   return ierr;
}
