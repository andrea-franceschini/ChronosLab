function param = parser_coarsen(node,prob_type)

if strcmpi(prob_type, 'mech')
   SoCtype = 'DOM';
   tau = 0.01;
elseif strcmpi(prob_type, 'lapl')
   SoCtype = 'CLA';
   tau = 0.25;
end

if isfield(node,'coarsen_tau')
   tau = str2double(getfield(node,'coarsen_tau'));
end

SoC_perc = 50;
nl_agg   = 0;
CR_iters = 1;

param.SoC_type = SoCtype;
param.tau      = tau;
param.SoC_perc = SoC_perc;
param.nl_agg   = nl_agg;
param.CR_iters = CR_iters;

end
