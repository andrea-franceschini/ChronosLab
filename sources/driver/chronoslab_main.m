% Adaptive Smoother and Prolongation AMG main driver
clear;
clc;
close all;

rand('state',0);

fprintf('EXECUTION BEGIN\n\n');

% Read names of the input files
fileIN = fopen('chronoslab.fnames','r');
C = textscan(fgetl(fileIN),'%s'); D = C{1}; file_XML    = D{1};
C = textscan(fgetl(fileIN),'%s'); D = C{1}; file_SOLVER = D{1};
C = textscan(fgetl(fileIN),'%s'); D = C{1}; file_BIN    = D{1};
fclose(fileIN);

% Read general parameters
[solv_method,itmax,tol,restart] = read_solver(file_SOLVER);

% Input from .mat file
load(file_BIN);

% Treat Dir nodes
A = treat_dirBC(A);

if exist('TV0')
   ntv0 = size(TV0,2);
else
   ntv0 = 1;
end
param = parse_xml(file_XML,ntv0);
param.symm = issymmetric(A);
param.VERBOSE = 10;

fprintf('END INPUT\n\n');

% Compute right-hand side
rhs = ones(size(A,1),1);

% Set timings (if needed)
amg_times

% Compute the AMG hierarchy
fprintf('BEGIN: Preconditioner computation\n');
AMG_prec = cpt_aspAMG(param,A,TV0);
fprintf('END: Preconditioner computation\n\n');

%-----------------------------------------------------------------------------------------

tic;
fprintf('BEGIN: System solution\n');
switch lower(solv_method)
    case 'stat_amg'

      % Solve the system by AMG stationary iteration
      fprintf('BEGIN: System solution by stationary AMG\n');
      [sol,flag,relres,iter,resvec] = stat_AMG(A,rhs,tol,itmax,AMG_prec);
      fprintf('END: System solution by stationary AMG\n');

    case 'pcg'

       % Solve the system by PCG
       fprintf('BEGIN: System solution by PCG\n');
       Mfun = @(r) AMG_Vcycle(AMG_prec,A,r);
       [sol,flag,relres,iter,resvec] = pcg(A,rhs,tol,itmax,Mfun);
       relres = resvec(end)/resvec(1);
       fprintf('END: System solution by PCG\n');

    case 'bicgstab'

      % Solve the system by BiCGstab
      fprintf('BEGIN: System solution by BiCGstab\n');
      Mfun = @(r) AMG_Vcycle(AMG_prec,A,r);
      [sol,flag,relres,iter,resvec] = bicgstab(A,rhs,tol,itmax,Mfun);
      fprintf('END: System solution by BiCGstab\n');

    case 'gmres'

      % Solve the system by BiCGstab
      fprintf('BEGIN: System solution by GMRES\n');
      Mfun = @(r) AMG_Vcycle(AMG_prec,A,r);
      [sol,flag,relres,iter,resvec] = gmres(A,rhs,restart,tol,itmax,Mfun);
      fprintf('END: System solution by GMRES\n');

end
fprintf('END: System solution\n');
T_iter = toc;

%-----------------------------------------------------------------------------------------

fprintf('\n');
fprintf('------------------------------------------------------------------------------------------');
fprintf('\n');

%% Get AMG hierarchy information
AMG_info = get_AMG_info(AMG_prec,A);

%% Print AMG hierarchy information
%print_AMG_info(AMG_info);
%fprintf('\n');

fprintf('Grid Complexity:     %10.4f\n',AMG_info(1).grid_comp);
fprintf('Operator Complexity: %10.4f\n',AMG_info(1).oper_comp);

% Print results
if flag
   fprintf('Convergence not achieved\n');
end
fprintf('\n');
fprintf('# of iterations:   %15d\n',iter);
fprintf('Relative residual: %15.6e\n',relres);
resreal = norm(rhs-A*sol)/norm(rhs);
fprintf('Real residual:     %15.6e\n',resreal);

fprintf('\n');
fprintf('AMG time:          %10.2f\n',T_smoo+T_tspa+T_coar+T_CR+T_prol);
fprintf('Iteration time:    %10.2f\n',T_iter);
