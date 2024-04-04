#include "Compute_nsy_rfsai.h"

int Compute_nsy_rfsai(const int nstep, const int step_size, const double eps,
                      const int nn_A, const int *iat_A, const int *ja_A,
                      const double *coef_A, int *&iat_FL, int *&ja_FL,
                      double *&coef_FL, int *&iat_FU, int *&ja_FU, double *&coef_FU){


   int ierr = 0;


   if (false){
      std::cout << "Printing input Matrix" << std::endl;
      FILE *ofile = fopen("mat_input.csr","w"); if (!ofile) exit(1);
      for (int i = 0; i < nn_A; i++){
          for (int j = iat_A[i]; j < iat_A[i+1]; j++){
             fprintf(ofile,"%10d %10d %25.15e\n",i+1,ja_A[j]+1,coef_A[j]);
          }
      }
      fclose(ofile);
   }



   double *diag_A = nullptr;
   extract_diag(nn_A,iat_A,ja_A,coef_A,diag_A);



   int nt_A = iat_A[nn_A];
   int *iat_Asym     = (int*) malloc((nn_A+1) * sizeof(int));
   int *ja_Asym      = (int*) malloc(nt_A * sizeof(int));
   double *coef_Asym = (double*) malloc(nt_A * sizeof(double));
   if (iat_Asym == nullptr || ja_Asym == nullptr || coef_Asym == nullptr){
      std::cout << "Error in allocating symmetrize pattern A" << std::endl;
      return 1;
   }
   memcpy(iat_Asym,iat_A,(nn_A+1)*sizeof(int));
   memcpy(ja_Asym,ja_A,nt_A*sizeof(int));
   memcpy(coef_Asym,coef_A,nt_A*sizeof(double));

   double *coef_AT  = nullptr;
   ierr = SymmetrizePattern(nn_A,iat_Asym,ja_Asym,coef_Asym,coef_AT);
   if (ierr != 0){
      std::cout << "Error in SymmetrizePattern" << std::endl;
      return 2;
   }


   if (false){
      std::cout << "Printing transposed input Matrix" << std::endl;
      FILE *ofile = fopen("mat_transposed.csr","w"); if (!ofile) exit(1);
      for (int i = 0; i < nn_A; i++){
          for (int j = iat_Asym[i]; j < iat_Asym[i+1]; j++){
             fprintf(ofile,"%10d %10d %25.15e\n",i+1,ja_Asym[j]+1,coef_AT[j]);
          }
      }
      fclose(ofile);
      std::cout << "Printing symmetrized input Matrix" << std::endl;
      ofile = fopen("mat_symmetrized.csr","w"); if (!ofile) exit(1);
      for (int i = 0; i < nn_A; i++){
          for (int j = iat_Asym[i]; j < iat_Asym[i+1]; j++){
             fprintf(ofile,"%10d %10d %25.15e\n",i+1,ja_Asym[j]+1,coef_Asym[j]);
          }
      }
      fclose(ofile);
   }



   ierr = cpt_nsy_rfsai(nstep,step_size,eps,nn_A,nt_A,diag_A,iat_Asym,ja_Asym,coef_Asym,
                        coef_AT,iat_FL,ja_FL,coef_FL,coef_FU);
   if (ierr != 0){
      std::cout << "Error in cpt_nsy_rfsai" << std::endl;
      return 3;
   }

   if (false){
      FILE *pfile = fopen("mat_FL.csr","w"); if (!pfile) exit(1);
      wrCSRmat(pfile,false,nn_A,iat_FL,ja_FL,coef_FL);
      fclose(pfile);
      pfile = fopen("mat_FU.csr","w"); if (!pfile) exit(1);
      wrCSRmat(pfile,false,nn_A,iat_FL,ja_FL,coef_FU);
      fclose(pfile);
   }



   int nt_FL;
   int nt_FU;
   ierr = compress_nsy_fsai(nn_A,nt_FL,nt_FU,iat_FL,ja_FL,iat_FU,ja_FU,coef_FL,coef_FU);
   if (ierr != 0){
      std::cout << "Error in compress_nsy_fsai" << std::endl;
      return 4;
   }


   free(diag_A);
   free(iat_Asym);
   free(ja_Asym);
   free(coef_Asym);
   free(coef_AT);

   return 0;

}
