function p2=polsum (p,gam,q)
%  function p2=poldif(p,gam,q) 
%  returns the poynomial q2 = p + gam * q; 
%  
nintv = size(p,2);
n = size(p,1);
m = size(q,1);
 
p2=zeros(max(m,n),nintv);
for ii=1:nintv 
 for j=1:n
   p2(j,ii)=p(j,ii);
 end
 for j=1:m
   p2(j,ii) = p2(j,ii) + gam*q(j,ii) ; 
 end;
end


