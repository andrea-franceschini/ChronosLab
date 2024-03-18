function param = parser_smoother(node,smooth_type)

maxNumCompThreads('automatic');
nthread = maxNumCompThreads('automatic');

nupre    = 1;
nupost   = 1;
avg_nnzr = 50;

switch lower(smooth_type)
   case 'jacobi'
      method = 'jacobi';
   case 'lightfsai'
      method    = 'afsai_cpp';
      nstep     = 15;
      step_size = 1;
      epsilon   = 1.e-3;
   case 'mediumfsai'
      method    = 'afsai_cpp';
      nstep     = 25;
      step_size = 1;
      epsilon   = 1.e-3;
   case 'heavyfsai'
      method    = 'afsai_cpp';
      nstep     = 25;
      step_size = 2;
      epsilon   = 1.e-3;
   case 'nestedfsai'
      method    = 'dbafsai_cpp';
      nstep     = 10;
      step_size = 1;
      epsilon   = 1.e-3;
      avg_nnzr  = 50;
   otherwise
      error('smooth_type not supported yet');
end

param.nthread   = nthread;
param.nstep     = nstep;
param.step_size = step_size;
param.epsilon   = epsilon;
param.avg_nnzr  = avg_nnzr;
param.method    = method;
param.nupre     = nupre;
param.nupost    = nupost;

end
