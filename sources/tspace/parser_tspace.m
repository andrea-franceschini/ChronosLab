function param = parser_tspace(node,ntv)

itmax       = 0;
tol         = 1.e-8;
method      = 'none';
init_approx = 0;

if isfield(node,'tspace_ntv')
   ntv = str2double(getfield(node,'tspace_ntv'));
end
if isfield(node,'tspace_iter')
   itmax = str2double(getfield(node,'tspace_iter'));
end

if itmax > 0
   method = 'srqcg';
end

% Check method and init_approx
switch lower(method)
   case 'none'
   case 'smoothing'
   case 'srqcg'
   case 'lanczos'
   case 'ng-srqcg'
   case 'ng-lanczos'
   case 'arnoldi'
   otherwise
      error('Not existing method');
end

param.ntv         = ntv;
param.itmax       = itmax;
param.tol         = tol;
param.method      = method;
param.init_approx = init_approx;

end
