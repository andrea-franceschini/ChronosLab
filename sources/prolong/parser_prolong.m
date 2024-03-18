function param = parser_prolong(node,prob_type)

maxNumCompThreads('automatic');
np = maxNumCompThreads('automatic');

if strcmpi(prob_type, 'mech')
   proltype = 'MEX_BAMG';
elseif strcmpi(prob_type, 'lapl')
   proltype = 'MEX_EXTI';
end
prol_emin    = 'EMIN';
itmax_Vol    = 100;
tol_Vol      = 0.01;
dist_min     = 1;
dist_max     = 15;
maxcond      = 1.e+13;
maxrownrm    = 5.0;
eps_prol     = 1.e-10;
updateCF     = 0;
patt_pow     = 2;
patt_tau     = 0.00;
nnzr_max     = 30;
itmax_emin   = 10;
energ_tol    = 0.001;
condmax_emin = 1.e+10;
maxwgt_emin  = 1.e+10;
prec_emin    = 1;
solv_emin    = 4;
min_lfil     = 0;
max_lfil     = 20;
D_lfil       = 5;

if isfield(node,'prolo_smooth')
   prol_emin = getfield(node,'prolo_smooth');
end

% Check that flags are ok
if ( ~( strcmpi(proltype,'MEX_BAMG') || strcmpi(proltype,'MEX_EXTI') || ...
        strcmpi(proltype,'MEX_CLAS') || strcmpi(proltype,'MEX_HYBC') ) )
   error('Wrong proltype flag: %s\n',proltype);
end
if ( ~( strcmpi(prol_emin,'EMIN') || strcmpi(prol_emin,'SMOOTH') || ...
        strcmpi(prol_emin,'NONE') ) )
   error('Wrong proltype flag: %s\n',prol_emin);
end

param.proltype     = proltype;
param.prol_emin    = prol_emin;
param.np           = np;
param.itmax_Vol    = itmax_Vol;
param.tol_Vol      = tol_Vol;
param.dist_min     = dist_min;
param.dist_max     = dist_max;
param.maxcond      = maxcond;
param.maxrownrm    = maxrownrm;
param.eps_prol     = eps_prol;
param.updateCF     = updateCF;
param.patt_pow     = patt_pow;
param.patt_tau     = patt_tau;
param.nnzr_max     = nnzr_max;
param.itmax_emin   = itmax_emin;
param.energ_tol    = energ_tol;
param.condmax_emin = condmax_emin;
param.maxwgt_emin  = maxwgt_emin;
param.prec_emin    = prec_emin;
param.solv_emin    = solv_emin;
param.min_lfil     = min_lfil;
param.max_lfil     = max_lfil;
param.D_lfil       = D_lfil;

end
