function [ppol] = initpol(intv) 
%% initializes the polynomial ppol
%% first interval 
ninterv = length(intv)-1;
for k=1:ninterv
  ppol(1,k) = 1.; 
  ppol(2,k) = 0.; 
end
