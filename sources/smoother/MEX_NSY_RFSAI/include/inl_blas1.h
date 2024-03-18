

#pragma once

#include <cmath>



inline void inl_dcopy(const int n, const double *const v1, const int k1,
                      double *const v2, const int k2){
   int j = 0;
   for (int i = 0; i < n; i += k1){
      v2[i] = v1[j];
      j += k2;
   }
   return;
}


inline double inl_ddot(const int n, const double *const v1, const int k1,
                     const double *const v2, const int k2){
   double ddot = 0.0;
   int j = 0;
   for (int i = 0; i < n; i += k1){
      ddot += v1[i]*v2[j];
      j += k2;
   }
   return ddot;
}


inline double inl_dnrm1(const int n, const double *const v1, const int k1){
   double dnrm1 = 0.0;
   for (int i = 0; i < n; i += k1) dnrm1 += std::abs(v1[i]);
   return dnrm1;
}


inline double inl_dnrm2(const int n, const double *const v1, const int k1){
   double dnrm2 = 0.0;
   for (int i = 0; i < n; i += k1) dnrm2 += (v1[i])*(v1[i]);
   return sqrt(dnrm2);
}


inline void inl_daxpy(const int n, const double alpha, const double *const x,
                      const int k1, double *const y, const int k2){
   int j = 0;
   for (int i = 0; i < n; i += k1){
      y[j] = y[j] + alpha*x[i];
      j += k2;
   }
   return;
}


