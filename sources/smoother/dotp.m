   function r = dotp(p, q) 
%% function r = dotp(p, q)
%% returns the inner product of two polynomials
%% whose expansions are given by p and q. 
%% dotprod = int1 + scal*int2
%% where int1 and int2 are integrals using chebyshev weights
%% in (0,a) and (a,b) respectively..
%% scal is set in the function itself -- 
scal =  [0.1, 1, 1, 1, 1.0] ; %% works for up to 5 intervals.

m = min(size(q,1),size(p,1));  
nintv = size(p,2);
r = 0.0;
if (m ~= 0) 
   for ii = 1:nintv 
       r = r + scal(ii)*q(1,ii)*p(1,ii) ; 
   end
end;
 
for j=1:m
   for ii =1:nintv
      r = r + scal(ii)*q(j,ii)*p(j,ii) ;
   end
end;

r = r*pi/2; 
