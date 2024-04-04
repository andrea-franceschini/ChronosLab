


#include "cpt_nsy_rfsai.h"

int cpt_nsy_rfsai(const int nstep, const int step_size, const double eps, const int nn_A,
                  const int nt_A, const double *diag_A, const int *iat_A, const int *ja_A,
                  const double *coef_A, const double *coef_AT, int *&iat_FL, int *&ja_FL,
                  double *&coef_FL, double *&coef_FUT){


   Open_DebugLog();


   int mmax = nstep*step_size;
   int ntmax_F = nn_A*(mmax+1) + nn_A;
   iat_FL   = (int*) malloc((nn_A+1) * sizeof(int));
   ja_FL    = (int*) malloc((ntmax_F) * sizeof(int));
   coef_FL  = (double*) malloc((ntmax_F) * sizeof(double));
   coef_FUT = (double*) malloc((ntmax_F) * sizeof(double));
   if ( iat_FL == nullptr ||  ja_FL == nullptr ||
        coef_FL == nullptr || coef_FUT == nullptr ) return 1;



   int *JWN = (int*) malloc(nn_A * sizeof(int));
   double *WR_L = (double*) malloc(nn_A * sizeof(double));
   double *WR_U = (double*) malloc(nn_A * sizeof(double));
   if ( JWN == nullptr || WR_L == nullptr || WR_U == nullptr ) return 2;

   double *full_A = (double*) malloc( (mmax*mmax) * sizeof(double));
   double *rhs_L = (double*) malloc( (mmax+1) * sizeof(double));
   double *rhs_U = (double*) malloc( (mmax+1) * sizeof(double));
   lapack_int *ipvt = (lapack_int*) malloc( mmax * sizeof(lapack_int));
   if (full_A == nullptr || rhs_L == nullptr || rhs_U == nullptr || ipvt == nullptr)
      return 2;
   double *rhs_L_sav = (double*) malloc( (mmax+1) * sizeof(double));
   double *rhs_U_sav = (double*) malloc( (mmax+1) * sizeof(double));
   if ( rhs_L_sav == nullptr || rhs_U_sav == nullptr ) return 2;


   std::fill_n(JWN,nn_A,0);


   int ind_FL = 0;
   iat_FL[0] = ind_FL;


   for( int irow = 0; irow < nn_A; irow++){


      if (DEBUG){
         fprintf(dbfile,"-------------------------------\n");
         fprintf(dbfile,"IROW: %d\n",irow);
      }


      int mrow = 0;


      int istep = 0;
      double DKap_old = 0.0;
      bool Refine = (nstep >= 1);
      while (Refine){

         istep++;

         if (DEBUG) fprintf(dbfile,"istep %6d mroww %6d\n",istep,mrow);



         int mrow_old = mrow;
         KapGrad_NSY(istep,irow,mrow,irow,step_size,iat_A,ja_A,coef_A,coef_AT,rhs_L,rhs_U,
                     &(ja_FL[ind_FL]),JWN,WR_L,WR_U);


         if (mrow > mrow_old){


            bool null_L;
            bool null_U;
            gather_fullsys(irow,mrow,&(ja_FL[ind_FL]),nn_A,iat_A,ja_A,coef_A,full_A,
                           rhs_L,rhs_U,null_L,null_U);

            if (DEBUG){
               fprintf(dbfile,"full_A:\n");
               for (int i = 0; i < mrow; i++){
                  for (int j = 0; j < mrow; j++) fprintf(dbfile," %15.6e",full_A[j*mrow+i]);
                  fprintf(dbfile,"\n");
               }
               fprintf(dbfile,"JCOLS: ");
               for (int i = 0; i < mrow; i++) fprintf(dbfile," %15d",ja_FL[ind_FL+i]);
               fprintf(dbfile,"\n");
               fprintf(dbfile,"RHS_L: ");
               for (int i = 0; i < mrow; i++) fprintf(dbfile," %15.6e",rhs_L[i]);
               fprintf(dbfile,"\n");
               fprintf(dbfile,"RHS_U: ");
               for (int i = 0; i < mrow; i++) fprintf(dbfile," %15.6e",rhs_U[i]);
               fprintf(dbfile,"\n");
            }



            if (!null_L || !null_U){
               lapack_int l_mrow = static_cast<lapack_int>( mrow );
               lapack_int info = LAPACKE_dgetrf(LAPACK_ROW_MAJOR,l_mrow,l_mrow,full_A,l_mrow,
                                                ipvt);
               if (info != 0) return 3;





            }



            for (int k = 0; k < mrow; k++) rhs_L_sav[k] = rhs_L[k];
            for (int k = 0; k < mrow; k++) rhs_U_sav[k] = rhs_U[k];


            if (!null_L){
               lapack_int l_mrow = static_cast<lapack_int>( mrow );
               lapack_int l_one = static_cast<lapack_int>( 1 );
               lapack_int info = LAPACKE_dgetrs(LAPACK_ROW_MAJOR,'T',l_mrow,1,full_A,l_mrow,
                                                ipvt,rhs_L,l_one);
               if (info != 0) return 3;
            }
            if (!null_U){
               lapack_int l_mrow = static_cast<lapack_int>( mrow );
               lapack_int l_one = static_cast<lapack_int>( 1 );
               lapack_int info = LAPACKE_dgetrs(LAPACK_ROW_MAJOR,'N',l_mrow,1,full_A,l_mrow,
                                                ipvt,rhs_U,l_one);
               if (info != 0) return 3;
            }

            if (DEBUG){
               fprintf(dbfile,"SOL_L: ");
               for (int i = 0; i < mrow; i++) fprintf(dbfile," %15.6e",rhs_L[i]);
               fprintf(dbfile,"\n");
               fprintf(dbfile,"SOL_U: ");
               for (int i = 0; i < mrow; i++) fprintf(dbfile," %15.6e",rhs_U[i]);
               fprintf(dbfile,"\n");
            }



            double DKap_L_new = inl_ddot(mrow,rhs_L,1,rhs_U_sav,1);
            double DKap_U_new = inl_ddot(mrow,rhs_U,1,rhs_L_sav,1);
            double DKap_new = fabs(DKap_L_new + DKap_U_new);


            if (istep == nstep){
               Refine = false;
            } else {
               Refine = (fabs(DKap_new-DKap_old) >= eps*DKap_old) && (DKap_new != 0.);
               DKap_old = fabs(DKap_new);
            }

         } else {


            Refine = false;

         }

      }


      exit_Refinement: ;


      double diag_entry = diag_A[irow];
      double scal_fac = diag_entry - 0.5*( inl_ddot(mrow,rhs_U,1,rhs_L_sav,1) +
                                           inl_ddot(mrow,rhs_L,1,rhs_U_sav,1) );

      if (DEBUG){
         fprintf(dbfile,"mrow: %d\n",mrow);
         fprintf(dbfile,"RHS_L: ");
         for (int i = 0; i < mrow; i++) fprintf(dbfile," %15.6e",rhs_L[i]);
         fprintf(dbfile,"\n");
         fprintf(dbfile,"RHS_U: ");
         for (int i = 0; i < mrow; i++) fprintf(dbfile," %15.6e",rhs_U[i]);
         fprintf(dbfile,"\n");
         fprintf(dbfile,"RHS_L_SAV: ");
         for (int i = 0; i < mrow; i++) fprintf(dbfile," %15.6e",rhs_L_sav[i]);
         fprintf(dbfile,"\n");
         fprintf(dbfile,"RHS_U_SAV: ");
         for (int i = 0; i < mrow; i++) fprintf(dbfile," %15.6e",rhs_U_sav[i]);
         fprintf(dbfile,"\n");
         fprintf(dbfile,"diag_entry: %15.6e\n",diag_entry);
         fprintf(dbfile,"f*b: %15.6e\n",inl_ddot(mrow,rhs_U,1,rhs_L_sav,1));
         fprintf(dbfile,"c*g: %15.6e\n",inl_ddot(mrow,rhs_L,1,rhs_U_sav,1));
         fprintf(dbfile,"SCAL FACTOR: %15.6e\n",scal_fac);
      }



      double check_val = fabs(scal_fac / diag_entry);
      //if ( check_val < 1.0e-10 ){
      //   cout << "SMALL DIAGONAL = " << check_val << " IN ROW: " << irow << endl;
      //}


      double fac = 1.0 / sqrt(fabs(scal_fac));
      for (int k = 0; k < mrow; k++) coef_FL[ind_FL+k] = fac*rhs_L[k];

      coef_FL[ind_FL+mrow] = fac;
      ja_FL[ind_FL+mrow] = irow;


      if (scal_fac < 0.0) fac = -fac;
      for (int k = 0; k < mrow; k++) coef_FUT[ind_FL+k] = fac*rhs_U[k];

      coef_FUT[ind_FL+mrow] = fac;


      for (int k = 0; k < mrow; k++) JWN[ja_FL[ind_FL+k]] = 0;


      ind_FL += mrow + 1;
      iat_FL[irow+1] = ind_FL;

   }


   free(JWN);
   free(WR_L);
   free(WR_U);
   free(full_A);
   free(rhs_L);
   free(rhs_U);
   free(ipvt);
   free(rhs_L_sav);
   free(rhs_U_sav);


   int nt_F = ind_FL;
   ja_FL    = (int*) realloc( ja_FL, (nt_F) * sizeof(int));
   coef_FL  = (double*) realloc( coef_FL , (nt_F) * sizeof(double));
   coef_FUT = (double*) realloc( coef_FUT , (nt_F) * sizeof(double));
   if ( iat_FL == nullptr ||  ja_FL == nullptr ||
        coef_FL == nullptr || coef_FUT == nullptr ) return 1;


   Close_DebugLog();

   return 0;
}
