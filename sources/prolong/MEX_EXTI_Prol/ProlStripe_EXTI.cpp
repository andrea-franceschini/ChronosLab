
#include <iostream>
#include <iomanip>
#define DBNODE  000
#define ENDNODE 200

#include <stdlib.h>   
#include <cmath>      
#include <algorithm>  
#include <limits>     

#include "ir_heapsort.h"

const double ONE  = 1.0;
const double ZERO = 0.0;
const double EPS = numeric_limits<double>::epsilon();

int ProlStripe_EXTI(const int firstrow, const int lastrow, const int nn_loc,
                    const int nn_A, const int nt_A, const int ntmax_P,
                    const int *const iat_A, const int *const ja_A,
                    const double *const coef_A, const int *const coef_S,
                    const int *const fcnodes, int &nn_P, int &nt_P, int *iat_P, int *ja_P,
                    double *coef_P){

   double avg_nnz = 1.2*static_cast<double>(nt_A) / static_cast<double>(lastrow-firstrow);
   int size_scr = static_cast<int>(avg_nnz)*static_cast<int>(avg_nnz);
   int *WI         = new int [nn_A](); if (WI == NULL) return 1;
   fill_n(WI,nn_A,0);
   int *ja_CC      = new int [size_scr](); if (ja_CC == NULL) return 1;
   int *ja_FS      = new int [size_scr](); if (ja_FS == NULL) return 1;
   int *pos_kj     = new int [size_scr](); if (pos_kj == NULL) return 1;
   double *coef_CC = new double [size_scr](); if (coef_CC == NULL) return 1;
   double *coef_FS = new double [size_scr](); if (coef_FS == NULL) return 1;
   double *a_kj    = new double [size_scr](); if (a_kj == NULL) return 1;

   nn_P=0;
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
         
         int n_FS = 0;
         int n_CC = 0;
         double a_ii;
         double denom = ZERO;
         for (int j = iat_A[inod]; j < iat_A[inod+1]; j++){ 
            int jcol = ja_A[j];
            if (jcol == inod){
               
               a_ii = coef_A[j];
               denom += a_ii;
            } else {
               if (coef_S[j] > 0){
                  
                  int jnod = fcnodes[jcol];
                  if (jnod >= 0){
                     
                     if (WI[jcol] > 0){
                        
                        
                        int pos = WI[jcol]-1;
                        coef_CC[pos] = coef_A[j];
                     } else if (WI[jcol] == 0){
                        
                        ja_CC[n_CC] = jcol;
                        coef_CC[n_CC] = coef_A[j];
                        WI[jcol] = n_CC+1;
                        n_CC++;
                     } else {
                        
                        
                        WI[jcol] = -WI[jcol];
                     }
                  } else {
                     
                     
                     ja_FS[n_FS] = jcol;
                     coef_FS[n_FS] = coef_A[j];
                     n_FS++;
                     
                     for (int k = iat_A[jcol]; k < iat_A[jcol+1]; k++){
                        if (coef_S[k] > 0){
                           int kcol = ja_A[k];
                           if (fcnodes[kcol] >= 0){
                              
                              if (WI[kcol] == 0){
                                 
                                 ja_CC[n_CC] = kcol;
                                 coef_CC[n_CC] = 0.0;
                                 
                                 WI[kcol] = n_CC+1;
                                 n_CC++;
                              } else if (WI[kcol] < 0){
                                 
                                 int pos = -(WI[kcol]+1);
                                 
                                 denom -= coef_CC[pos];
                                 WI[kcol] = -WI[kcol];
                              }
                           }
                        }
                     }
                  }
               } else {
                  
                  
                  if (WI[jcol] == 0){
                     denom += coef_A[j];
                     
                     if (fcnodes[jcol] >= 0){
                        
                        ja_CC[n_CC] = jcol;
                        
                        coef_CC[n_CC] = coef_A[j];
                        WI[jcol] = -(n_CC+1);
                        n_CC++;
                     }
                  } else {
                     
                     int pos = WI[jcol] - 1;
                     coef_CC[pos] = coef_A[j];
                  }
               }
            }
         } 
         
         int n_int = 0;
         for (int i = 0; i < n_CC; i++){
            int jcol = ja_CC[i];
            if (WI[jcol] > 0){
               
               ja_P[ind_P+n_int] = jcol;
               coef_P[ind_P+n_int] = coef_CC[i];
               WI[jcol] = n_int+1;
               n_int++;
               
            } else {
               
               WI[jcol] = 0;
            }
         }
         
         for (int k = 0; k < n_FS; k++){
            int knod = ja_FS[k];
            double a_ik = coef_FS[k];
            
            double ext_sum = min(ZERO,a_ik);
            
            
            
            
            int ind = 0;
            for (int l = iat_A[knod]; l < iat_A[knod+1]; l++){
               int lcol = ja_A[l];
               if ( WI[lcol] > 0){
                  
                  double a_kl = min(ZERO,coef_A[l]);
                  
                  ext_sum += a_kl;
                  a_kj[ind] = a_kl;
                  pos_kj[ind] = WI[lcol] - 1;
                  ind++;
               }
            }
            
            if ( abs(ext_sum) > EPS*a_ii ){
               
               denom += a_ik*min(ZERO,a_ik) / ext_sum;
               
               for (int jj = 0; jj < ind; jj++){
                  int jpos = pos_kj[jj];
                  coef_P[ind_P+jpos] += a_ik*a_kj[jj] / ext_sum;
               }
               
            } else {
               
               double fac = ONE / (static_cast<double>(n_int+1));
               
               denom += a_ik*fac;
               
               for (int jj = 0; jj < ind; jj++){
                  int jpos = pos_kj[jj];
                  coef_P[ind_P+jpos] += a_ik*fac;
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

      }

      iat_P[nn_P] = ind_P;

   } 
   
   nt_P = ind_P;
   
   delete [] WI;
   delete [] ja_CC;
   delete [] ja_FS;
   delete [] pos_kj;
   delete [] coef_CC;
   delete [] coef_FS;
   delete [] a_kj;

   return 0;
}
