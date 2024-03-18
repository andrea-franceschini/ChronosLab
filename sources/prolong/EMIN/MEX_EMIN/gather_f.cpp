void gather_f(const int nn, const int *indices, const int len, const int *jcol,
              const double *coef, double *vec_f){

   int ii = 0;
   int jj = 0;
   if (len > 0){
      while (ii < nn){

         while (jcol[jj] < indices[ii]){
            jj++;
            if (jj == len) goto exit_loop;
         }
         if (jcol[jj] == indices[ii]){
            vec_f[ii] = -coef[jj];
         } else {
            vec_f[ii] = 0.0;
         }
         ii++;
      }
   }
   exit_loop: ;

   for (int k = ii; k < nn; k++) vec_f[k] = 0.0;



}
