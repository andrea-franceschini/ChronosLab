#include "compress_nsy_fsai.h"

int compress_nsy_fsai(const int nn, int &nt_FL, int &nt_FU,
                      int *&iat_FL, int *&ja_FL, int *&iat_FU, int *&ja_FU,
                      double *&coef_FL, double *&coef_FU){


   nt_FL = iat_FL[nn];
   nt_FU = nt_FL;



   int *iat_FL_new = (int*) malloc((nn+1) * sizeof(int));
   int *ja_FL_new = (int*) malloc(nt_FL * sizeof(int));
   double *coef_new = (double*) malloc(nt_FL * sizeof(double));
   iat_FU = (int*) malloc((nn+1) * sizeof(int));
   if (iat_FL_new == nullptr || ja_FL_new == nullptr || coef_new == nullptr ||
       iat_FU == nullptr){

      std::cout << "Allocation Error in compress_nsy_fsai" << std::endl;
      return 1;
   }


   std::fill_n(iat_FU,nn+1,0);


   nt_FL = 0;
   nt_FU = 0;
   iat_FL_new[0] = 0;
   for (int i = 0; i < nn; i++){
      for (int j = iat_FL[i]; j < iat_FL[i+1]; j++){
         if (coef_FL[j] != 0.0){
            coef_new[nt_FL] = coef_FL[j];
            ja_FL_new[nt_FL] = ja_FL[j];
            nt_FL++;
         }
         if (coef_FU[j] != 0.0) {
            iat_FU[ja_FL[j]]++;
            nt_FU++;
         }
      }
      iat_FL_new[i+1] = nt_FL;
   }
   nt_FL = nt_FL;
   nt_FU = nt_FU;






   free(coef_FL);
   coef_FL = coef_new;


   ja_FL_new = (int*) realloc( ja_FL_new, (nt_FL) * sizeof(int));
   coef_FL  = (double*) realloc( coef_FL , (nt_FL) * sizeof(double));
   if ( ja_FL == nullptr || coef_FL == nullptr ) return 1;


   ja_FU = (int*) malloc(nt_FU * sizeof(int));
   coef_new  = (double*) malloc( nt_FU * sizeof(double));
   int *ISCR = (int*) malloc((nn+1) * sizeof(int));
   if ( ja_FU == nullptr || coef_new == nullptr || ISCR == nullptr) return 1;


   ISCR[0] = 0;
   for ( int i = 1; i < nn+1; i++ ) ISCR[i] = ISCR[i-1] + iat_FU[i-1];
   for ( int i = 0; i < nn+1; i++ ) iat_FU[i] = ISCR[i];




   for ( int i = 0; i < nn; i++ ){
      for ( int j = iat_FL[i]; j < iat_FL[i+1]; j++ ){
         if (coef_FU[j] != 0.0){
            int ind  = ISCR[ja_FL[j]];
            ja_FU[ind] = i;
            coef_new[ind] = coef_FU[j];
            ISCR[ja_FL[j]] = ind+1;
         }
      }
   }


   free(coef_FU);
   coef_FU = coef_new;


   free(ISCR);


   free(iat_FL);
   free(ja_FL);
   iat_FL = iat_FL_new;
   ja_FL = ja_FL_new;

   return 0;

}
