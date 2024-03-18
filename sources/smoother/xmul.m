 function q=xmul(p,intv) 
%  function q=xmul(p, intv) 
%  returns the coefficents of the polynomial expansion
%  t*  sum  p_i * T_i ( alp* t - bet) in the basis
%  { T_i( alp * t - bet) } 
%  p = input polynomial -- intv = intervals.. 

 q = 0;
nintv = length(intv)-1;

for ii = 1:nintv 
   c = 0.5* (intv(ii) + intv(ii+1));
   h = 0.5* (intv(ii+1) - intv(ii));
   m = size(p,1) ; 
%%
   if (m == 0)
      return
   end	
%%
   q(1:m,ii) = c*p(1:m,ii); 
   q(m+1,ii) = 0;
%%  contribution of   t * T_0(alp t - bet) = t  = [ s+bet]/alp
   q(2,ii) = q(2,ii) + h*p(1,ii) ;
%%
   h2 = 0.5*h;
   for j=2:m 
     q(j+1,ii) = q(j+1,ii) + h2*p(j,ii);
     q(j-1,ii) = q(j-1,ii) + h2*p(j,ii);
   end;
end
