function [AMG_hrc] = cpt_aspAMG(param,A,TV0)
%-----------------------------------------------------------------------------------------
%
% Function representing the entry point for the recursive computation of the AMG hierarchy
%
% Input
%
% param:     structure containing all the parameters for AMG construction
% A:         system matrix
% TV0:       initial testspace which is eventually padded to the desired size
%
% Output
%
% AMG_hrc:  AMG hierachy, a linked list of all the AMG components at every level
%
%-----------------------------------------------------------------------------------------

% Init the level counter
level = 0;

VERBOSE = param.VERBOSE;

% Padd the test space with random vectors if TV0 is not large enough
if VERBOSE > 2
   fprintf('BEGIN: Initializing initial test space\n');
end
n = size(A,1);
ntv = param.tspace.ntv;
ntv0 = size(TV0,2);
TV = zeros(n,ntv);
TV(:,1:ntv0) = TV0;
TV(:,ntv0+1:ntv) = rand(n,ntv-ntv0);

% Orthonormalize the test space
[TV,~] = qr(TV,0);
if VERBOSE > 2
   fprintf('END: Initializing initial test space\n');
end

AMG_hrc = REC_cpt_aspAMG(level,param,A,TV);

end

%-----------------------------------------------------------------------------------------

function [AMG_hrc] = REC_cpt_aspAMG(level,param,A,TV0)

global DEBINFO;
VERBOSE = param.VERBOSE;
% Remove some DEBUG file
%if DEBINFO.flag
%   if level == 1
%      ! rm Logfile.*
%      ! rm IBDEBUG
%      ! rm ROWDEBUG
%      ! rm LOG_CORRECT
%   end
%end

amg_times

% Increase the level counter
level = level + 1;
AMG_hrc.level = level;
if VERBOSE >= 1
   fprintf('*****   LEVEL: %5d   ***** \n',level);
end

% Get the current size
n = size(A,1);
% Compute operator max eigenvalue
lmax_A = eigs(A,1,'lm','Tolerance',1.e-4);
AMG_hrc.lmax_A = lmax_A;

% Print operator info
if VERBOSE > 2
   fprintf('Operator # of rows:               %10d\n',n);
   fprintf('Operator # of non-zeroes:         %10d\n',nnz(A));
   fprintf('Operator # of non-zeroes per row: %10.2f\n',nnz(A)/n);
   fprintf('Operator max lambda:              %15.6e\n\n',lmax_A);
end

% Compute the preconditioner for this level
if n <= param.amg.maxCoarseSZ || level >= param.amg.nLevMax

   %***** This level is the last one *****

   % Compute the Factorization and store it
   if param.symm
      [L,p,S] = chol(A,'lower');
      if p ~= 0
         fprintf('ERROR IN FACTORIZING LAST LEVEL\n');
      end
      AMG_hrc.L = L;
      AMG_hrc.S = S;
   else
      if VERBOSE > 2
         fprintf('Using UNSYMMETRIC factorization\n');
      end
      [L,U,P,Q] = lu(A);
      % P*A*Q = L*U ---> A*x = b ---> P'*LU*Q'*x = b ---> LU*Q'*x = P*b ---> x = Q*(U\(L\(P*b))
      AMG_hrc.L = L;
      AMG_hrc.U = U;
      AMG_hrc.P = P;
      AMG_hrc.Q = Q;
   end
   AMG_hrc.symm = param.symm;
   AMG_hrc.A = A;
   AMG_hrc.next = 0;

else

   %***** Compute this level using adaptive AMG *****

   % Compute the smoother
   tic;
   if VERBOSE > 2
      fprintf('BEGIN: Computing the smoother\n');
   end
   smootherOp = smoother(A, param.smoother, VERBOSE);
   if VERBOSE > 2
      fprintf('END: Computing the smoother\n\n');
   end
   T_smoo = T_smoo + toc;

   %--------------------------------------------------------------------------------------

   % Compute the test space
   tic;
   if VERBOSE > 2
      fprintf('BEGIN: Computing the test space\n');
   end
   [TV, lambda, res] = tspace(TV0, A, smootherOp, param.tspace, VERBOSE);
   if VERBOSE > 2
      fprintf('   i          lambda     res_lam     res_vec\n');
   end
   nl = numel(lambda);
   if VERBOSE > 2
      fprintf('%4d %15.6e %11.2e %11.2e\n',[(1:nl)' lambda(end:-1:1) res(end:-1:1,:)]');
      fprintf('END: Computing the test space\n\n');
   end
   T_tspa = T_tspa + toc;

   %--------------------------------------------------------------------------------------

   % Check if aggressive coarsening is required
   use_aggressive = level <= param.coarsen.nl_agg;

   % Compute coarsening
   tic;
   if VERBOSE > 2
      fprintf('BEGIN: Computing the coarse nodes\n');
   end
   if ~use_aggressive
      [fcnode,clist,flist,S,S_patt] = coarsen(param.coarsen,A,smootherOp,TV,VERBOSE);
      num_f = numel(flist);
      num_c = numel(clist);
      fclist = 0;
      fflist = 0;
   else
      [fcnode,clist,fclist,fflist,S,S_patt] = agg_coarsen(param.coarsen,A,TV,VERBOSE);
      num_f = numel(fclist) + numel(fflist);
      num_c = numel(clist);
   end
   if VERBOSE > 2
      fprintf('END: Computing the coarse nodes\n\n');
   end
   T_coar = T_coar + toc;
   if VERBOSE > 2
      fprintf('Number of FINE nodes:   %10d\n',num_f);
      fprintf('Number of COARSE nodes: %10d\n',num_c);
      fprintf('Coarse node percentage: %10.2f\n\n',100*num_c/numel(fcnode));
   end
   if num_f + num_c ~= numel(fcnode)
      err_msg = 'Number of F/C node is not consistent with level size';
      error(err_msg);
   end

   %--------------------------------------------------------------------------------------

   % Perform compatible relaxations to estimate coarse grid quality
   tic;
   if param.coarsen.CR_iters > 0
      if VERBOSE > 2
         fprintf('BEGIN: Compatible relaxation\n');
      end
      rho = comp_relax(param.coarsen.CR_iters,fcnode,smootherOp,A);
      if VERBOSE > 2
         fprintf('END: Compatible relaxation\n\n');
      end
   end
   T_CR = T_CR + toc;
   if VERBOSE > 2
      fprintf('Number of FINE nodes:   %10d\n',num_f);
   end

   %--------------------------------------------------------------------------------------

   % Compute prolongation
   tprol = tic;
   if VERBOSE > 2
      fprintf('BEGIN: Computing the prolongation\n');
   end
   [P,clist,fcnode,emin_info] = prolong(use_aggressive,level,param,clist,fclist,fflist,...
                                        fcnode,S,S_patt,smootherOp,A,TV);
   if VERBOSE > 2
      fprintf('END: Computing the prolongation\n\n');
   end
   prol_tot_time = toc(tprol);
   if VERBOSE > 2
      fprintf('Prolongation time %10.3f [s]\n',prol_tot_time);
   end
   T_prol = T_prol + prol_tot_time;

   % Store EMIN detailed information
   if ( strcmpi(param.prolong.prol_emin,'EMIN') )
      emin_info.time_prolTot  = prol_tot_time;
   end
   AMG_hrc.emin_info = emin_info;

   % Create next level test space (by simple injection of coarse nodes)
   nc = numel(clist);
   TVnext(1:nc,:) = TV(clist,:);

   %--------------------------------------------------------------------------------------

   % Store the components of the current level
   if level > 1
      AMG_hrc.A = A;
   end
   AMG_hrc.P = P;
   AMG_hrc.S = S;
   AMG_hrc.S_patt = S_patt;
   AMG_hrc.TV = TV;
   AMG_hrc.fcnode = fcnode;
   AMG_hrc.omega =smootherOp.omega;
   AMG_hrc.Snnz = nnz(smootherOp.right) + nnz(smootherOp.left);
   AMG_hrc.nupre = param.smoother.nupre;
   AMG_hrc.nupost = param.smoother.nupost;
   AMG_hrc.smt_lft = smootherOp.left;
   AMG_hrc.smt_rgt = smootherOp.right;
   if smootherOp.LS_deg > 0
      AMG_hrc.Minv = smootherOp.polyPrec;
   else
      if (numel(smootherOp.left_out) + numel(smootherOp.right_out)) == 0
         % Simple smoother
         AMG_hrc.Minv = @(x) smootherOp.omega*(smootherOp.right*(smootherOp.left*x));
      else
         AMG_hrc.Snnz = AMG_hrc.Snnz + nnz(smootherOp.left_out) + nnz(smootherOp.right_out);
         % Nested smoother
         AMG_hrc.Minv = @(x) smootherOp.omega*(smootherOp.right*(smootherOp.right_out*...
                             (smootherOp.left_out*(smootherOp.left*x))));
      end
   end

   % Create next level operator
   Anext = P'*(A*P);
   % Enforce symmetry if a symmetric accelerator is required
   if param.symm
      Anext = 0.5*(Anext+Anext');
   end

   %--------------------------------------------------------------------------------------

   % Compute next level of the hierarchy
   AMG_hrc.next = REC_cpt_aspAMG(level,param,Anext,TVnext);

end

end
