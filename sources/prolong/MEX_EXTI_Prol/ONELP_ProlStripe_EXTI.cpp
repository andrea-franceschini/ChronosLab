#include <iostream>
#include <iomanip>
#define DBNODE  7
#define ENDNODE 7

#include <stdlib.h>
#include <cmath>
#include <algorithm>
#include <limits>

#include "ir_heapsort.h"

const double ONE  = 1.0;
const double ZERO = 0.0;
const double EPS = std::numeric_limits<double>::epsilon();

int ONELP_ProlStripe_EXTI(const int firstrow, const int lastrow, const int nn_loc,
                          const int nn_A, const int nt_A, const int ntmax_P,
                          const int *const iat_A, const int *const ja_A,
                          const double *const coef_A, const int *const coef_S,
                          const int *const iat_C, const int *const ja_C,
                          const double *const coef_C, const int *const fcnodes,
                          int &nn_P, int &nt_P, int *iat_P, int *ja_P,
                          double *coef_P){

   double avg_nnz = 20.0*static_cast<double>(nt_A) / static_cast<double>(lastrow-firstrow);
   int size_scr = static_cast<int>(avg_nnz);
   int *WI         = new int [nn_A](); if (WI == NULL) return 1;
   std::fill_n(WI,nn_A,0);
   int *ja_FS      = new int [size_scr](); if (ja_FS == NULL) return 1;
   double *coef_FS = new double [size_scr](); if (coef_FS == NULL) return 1;
   int *list_weak  = new int [size_scr](); if (list_weak == NULL) return 1;


   nn_P = 0;
   int ind_P = 0;
   iat_P[0] = ind_P;


   int shift = firstrow - 1;
   for (int inod = firstrow; inod < lastrow ; inod ++){

      int inod_coarse = fcnodes[inod];

      if (inod_coarse >= 0){
         nn_P++;


         ja_P[ind_P] = inod_coarse;
         coef_P[ind_P] = ONE;
         ind_P++;

      } else if (inod_coarse == -1) {
         nn_P++;

         int n_int = 0;
         int n_weak = 0;
         double denom = ZERO;
         double a_ii;
         for (int j = iat_A[inod]; j < iat_A[inod+1]; j++){
            int jcol = ja_A[j];
            if (jcol == inod){

               a_ii = coef_A[j];
               denom += a_ii;
            } else {
               if (fcnodes[jcol] >= 0){

                  if (coef_S[j] > 0){

                     ja_P[ind_P+n_int] = jcol;
                     coef_P[ind_P+n_int] = coef_A[j];
                     WI[jcol] = n_int+1;
                     n_int++;

                  } else {

                     denom += coef_A[j];
                     WI[jcol] = -j;
                     list_weak[n_weak] = jcol;
                     n_weak++;
                  }
               } else {

                  if (coef_S[j] < 0){

                     denom += coef_A[j];
                  }
               }
            }
         }

         for (int j = iat_C[inod]; j < iat_C[inod+1]; j++){
            int jcol = ja_C[j];
            if (jcol != inod){

               if (fcnodes[jcol] < 0){
                  double a_ik = coef_C[j];
                  double a_ki_bar;
                  double ext_sum = 0.0;
                  int n_FS = 0;

                  for (int k = iat_C[jcol]; k < iat_C[jcol+1]; k++){
                     int kcol = ja_C[k];
                     if (kcol == inod){

                        a_ki_bar = std::min(ZERO,coef_C[k]);
                        ext_sum += a_ki_bar;

                     } else {

                        if (fcnodes[kcol] >= 0){

                           double a_kj_bar = std::min(ZERO,coef_C[k]);
                           ja_FS[n_FS] = kcol;
                           coef_FS[n_FS] = a_kj_bar;
                           n_FS++;
                           ext_sum += a_kj_bar;

                           if ( WI[kcol] <= 0 ){
                              double value = ZERO;
                              if (WI[kcol] < 0){

                                 value = coef_A[-WI[kcol]];

                                 denom -= value;
                              }

                              ja_P[ind_P+n_int] = kcol;
                              coef_P[ind_P+n_int] = value;
                              WI[kcol] = n_int+1;
                              n_int++;
                           }
                        }
                     }
                  }

                  if ( std::abs(ext_sum) > EPS*a_ii ){

                     for (int k = 0; k < n_FS; k++){
                        int pos = WI[ja_FS[k]]-1;
                        coef_P[ind_P+pos] += (a_ik*coef_FS[k]) / ext_sum;
                     }
                     denom += (a_ik * a_ki_bar) / ext_sum;

                  } else {

                     double fac = ONE / (static_cast<double>(n_int+1));
                     denom += a_ik*fac;
                     for (int k = 0; k < n_FS; k++){
                        int pos = WI[ja_FS[k]]-1;
                        coef_P[ind_P+pos] += a_ik*fac;
                     }

                  }
               }
            }
         }

         for (int i = ind_P; i < ind_P+n_int; i++){

            WI[ja_P[i]] = 0;
            ja_P[i] = fcnodes[ja_P[i]];
            coef_P[i] = -coef_P[i] / denom;
         }
         ir_heapsort(&(ja_P[ind_P]),&(coef_P[ind_P]),n_int);

         ind_P += n_int;
         for (int i = 0; i < n_weak; i++) WI[list_weak[i]] = 0;
      }

      iat_P[nn_P] = ind_P;
   }

   nt_P = ind_P;

   delete [] WI;
   delete [] ja_FS;
   delete [] coef_FS;
   delete [] list_weak;

   return 0;
}
