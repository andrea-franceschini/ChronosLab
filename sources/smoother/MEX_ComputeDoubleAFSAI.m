function [mat_G1,mat_G2,nnzrs,times] = MEX_ComputeDoubleAFSAI(nthreads,PRINT_MAT,...
                                       nstep,step_size,eps,avg_nnzr,mat_A)

% Unpack the input matrix
nn_A = size(mat_A,1);
nt_A = nnz(mat_A);
[iat_A,ja_A,coef_A] = unpack_csr(mat_A);

% Switch from column to row
iat_A     = iat_A';
ja_A      = ja_A';
coef_A    = coef_A';

% Switch from matlab type to c++ type
nthreads  = int32(nthreads);
PRINT_MAT = int32(PRINT_MAT);
nstep     = int32(nstep);
step_size = int32(step_size);
iat_A     = int64(iat_A) - 1;
ja_A      = int32(ja_A) - 1;

% Compute DoubleFSAI preconditioner ------------------------------------------------------
[iat_G1,ja_G1,coef_G1,iat_G2,ja_G2,coef_G2,nnzrs,times] = ...
        Entry_ComputeDoubleFSAI(nthreads,PRINT_MAT,nstep,step_size,eps,avg_nnzr,...
        nn_A,iat_A,ja_A,coef_A);

% Create the sparse FSAI factors
% G1 --->
nt_G1 = size(ja_G1,2);
irow_G1 = zeros(nt_G1,1);
iend   = iat_G1(1)-1;
for i = 1:nn_A
   istart = iend + 1;
   iend = iat_G1(i+1)-1;
   irow_G1(istart:iend) = i;
end
ja_G1 = double(ja_G1);
mat_G1 = sparse(irow_G1,ja_G1,coef_G1,nn_A,nn_A);
% G2 --->
nt_G2 = size(ja_G2,2);
irow_G2 = zeros(nt_G2,1);
iend   = iat_G2(1)-1;
for i = 1:nn_A
   istart = iend + 1;
   iend = iat_G2(i+1)-1;
   irow_G2(istart:iend) = i;
end
ja_G2 = double(ja_G2);
mat_G2 = sparse(irow_G2,ja_G2,coef_G2,nn_A,nn_A);

end
