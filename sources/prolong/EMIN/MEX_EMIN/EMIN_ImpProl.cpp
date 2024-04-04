#include <iostream>
#include <chrono>

#if defined PRINT
#define dump true
#else
#define dump false
#endif

#include "parm_EMIN.h"
#include "wrCSRmat.h"
#include "Transp_Patt.h"
#include "EMIN_matfree.h"
#include "copy_Prol.h"
#include "gather_B_QR.h"
#include "gather_f.h"
#include "cpt_Trace_Acc.h"
#include "print_Q.h"
#include "print_Z.h"

int EMIN_ImpProl(const int np, const int itmax, const double en_tol, const double condmax,
                 const double maxwgt, const int prec_type, const int sol_type,
                 const int min_lfil, const int max_lfil, const int D_lfil, const int nn,
                 const int nn_C, const int ntv, const int nt_A, const int nt_P,
                 const int nt_patt, const int *fcnode, const int *iat_A, const int *ja_A,
                 const double *coef_A, const int *iat_Pin, const int *ja_Pin,
                 const double *coef_Pin, const int *iat_patt, const int *ja_patt,
                 const double *const *TV, int *&iat_Pout, int *&ja_Pout, double *&coef_Pout,
                 double *info)
{


   int ierr = 0;

   int iter = 0;
   int nnz_K = 0, nnz_PK = 0, nnz_ZQ = 0;
   double avg_lfil = 0.0;
   double relres = 0.0;

   std::chrono::time_point<std::chrono::system_clock> start, end, glob_start, glob_end;
   std::chrono::duration<double> elaps_sec;

   glob_start = std::chrono::system_clock::now();

   ierr = EMIN_matfree(np,itmax,en_tol,condmax,maxwgt,min_lfil,max_lfil,D_lfil,
                       nn,nn_C,ntv,nt_A,nt_P,nt_patt,fcnode,iat_A,ja_A,coef_A,
                       iat_Pin,ja_Pin,coef_Pin,iat_patt,ja_patt,TV,iat_Pout,
                       ja_Pout,coef_Pout,info);

   glob_end = std::chrono::system_clock::now();
   elaps_sec = glob_end - glob_start;
   double time_glob =  elaps_sec.count();

   info[4]  = 0.0;
   info[5]  = time_glob;

   return ierr;

}
