int DEFL_PCG_matfree(const int np, const int nn, const int nn_C, const int nn_K,
                     const int ntv, const int *perm, const int *iperm,
                     const double *D_inv, const int *iat_A, const int *ja_A,
                     const double *coef_A, const double Tr_A, const int *iat_patt,
                     const int *ja_patt, const double *mat_Q, const double *vec_P0,
                     const double *vec_f, const int itmax, const double energy_tol,
                     int &iter, double *vec_DP);
