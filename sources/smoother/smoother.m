function smootherOp = smoother(A, param, VERBOSE)

%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% Parameters for polynomial acceleration (temporarily hard coded...)
fac_max = 1.01;
LS_deg = 0;
disp = 0;
%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

% Mandatory parameters
nthread      = param.nthread;
nstep        = param.nstep;
step_size    = param.step_size;
epsilon      = param.epsilon;
avg_nnzr     = param.avg_nnzr;
method       = param.method;

% Init the outer smoother to void
smootherOp.left_out = [];
smootherOp.right_out = [];

switch lower(method)

    case 'afsai_cpp'
        % Set-up AFSAI (afsai with mex-cpp code)
        F = afsai_cpp(A,nthread,nstep,step_size,epsilon);
        % Correct NaN for very ill-conditioned problems
        irow = find(isnan(diag(F)));
        if VERBOSE > 2 && numel(irow) > 0
           fprintf('WARNING: Correcting %d diagonals\n',numel(irow));
        end
        for i = 1:numel(irow)
           ii = irow(i);
           F(ii,:) = 0;
           F(ii,ii) = 1 / sqrt(A(ii,ii));
           A(ii,ii)
        end
        % Compute damping parameter
        FAFT = @(x) F*(A*(F'*x));
        opts.issym = 1;
        opts.disp = disp;
        opts.tol = 5.e-4;
        lambda = eigs(FAFT,size(A,1),1,'la',opts);
        if VERBOSE > 2
           fprintf('Max Lambda: %10.4f\n',lambda);
        end
        omega = min(1,1.9 / lambda);
        % Append the smoother
        smootherOp.left = F;
        smootherOp.right = F';
        smootherOp.omega = omega;
        smootherOp.lambda = lambda;

    case 'afsai_nsy'
        % Set-up AFSAI_NSY (afsai for nsy systems with mex-cpp code)
        [FL,FU] = NSY_rfsai_cpp(nstep,step_size,epsilon,A);
        % Compute damping parameter
        FAFT = @(x) FL*(A*(FU*x));
        opts.issym = 0;
        opts.disp = 1;
        opts.tol = 5.e-4;
        lambda = eigs(FAFT,size(A,1),1,'lm',opts);
        if VERBOSE > 2
           fprintf('Max Lambda: %10.4f\n',lambda);
        end
        omega = min(1,1.9 / lambda);
        % Append the smoother
        smootherOp.left = FL;
        smootherOp.right = FU;
        smootherOp.omega = omega;
        smootherOp.lambda = lambda;

    case 'jacobi'
        % Compute Diagonal
        F = 1 ./ sqrt(full(diag(A)));
	    F = diag(sparse(F));
        % Compute damping parameter
        FAFT = @(x) F*(A*(F'*x));
        opts.issym = 1;
        if true
           lambda = eigs(FAFT,size(A,1),1,'lm','IsFunctionSymmetric',1,'Tolerance',1.e-2,'Display',disp,'FailureTreatment','keep');
           if VERBOSE > 2
              fprintf('Max Lambda: %10.4f\n',lambda);
           end
           omega = min(1,1.9 / lambda);
        else
           lambda = 2.0;
           omega = 1.0;
        end
        % Append the smoother
        smootherOp.left = F;
        smootherOp.right = F';
        smootherOp.omega = omega;
        smootherOp.lambda = lambda;

    case 'dbafsai_cpp'
        % Set-up DBAFSAI (double AFSAI with mex-cpp code)
        PRINT_MAT = 0;
        [G1,G2,nnzrs,times] = MEX_ComputeDoubleAFSAI(nthread,PRINT_MAT,...
                              nstep,step_size,eps,avg_nnzr,A);
        % Compute damping parameter
        FAFT = @(x) G2*(G1*(A*(G1'*(G2'*x))));
        opts.issym = 1;
        opts.disp = disp;
        opts.tol = 5.e-4;
        if true
           lambda = eigs(FAFT,size(A,1),1,'lm',opts);
           if VERBOSE > 2
              fprintf('Max Lambda: %10.4f\n',lambda);
           end
           omega = min(1,1.9 / lambda);
        else
           lambda = 2.0;
           omega = 1.0;
        end
        % Append the smoother
        smootherOp.left = G1;
        smootherOp.right = G1';
        smootherOp.left_out = G2;
        smootherOp.right_out = G2';
        smootherOp.omega = omega;
        smootherOp.lambda = lambda;

    otherwise
        error('Not existing method');
end

% Compute polynomial acceleration if needed
if LS_deg > 0
   nn = size(A,1);
   % Perturbe eigenvalues to prevent bad estimates
   lmin = 0.0;
   lmax = lambda*fac_max;
   % Compute polynomial acceleration coefficients
   [alpha, beta, gamma] = cpt_LSpol([lmin lmax],LS_deg);
   % Check damping factor
   prec_inner = @(x) smootherOp.left*(A*(smootherOp.right*x));
   prec = @(x) smootherOp.right*(apply_LSpol(prec_inner,(smootherOp.left*x),LS_deg,alpha, beta, gamma));
   MA = @(x) prec(A*x);
   lmax_poly = eigs(MA,nn,1,'lr','Display',disp,'Tolerance',0.001);
   if VERBOSE > 2
      fprintf('Post-PolyAcc - Max lambda: %e\n',lmax_poly);
   end
   omega_poly = 1.9 / lmax_poly;
   omega_poly = min(omega_poly,1);
   if VERBOSE > 2
      fprintf('omega: %e\n',omega_poly);
   end
   smootherOp.polyPrec = prec;
end

% Store polynomial degree
smootherOp.LS_deg = LS_deg;

end
