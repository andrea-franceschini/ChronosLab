function x = apply_LSpol(prod_A,b,deg,alpha,beta,gamma)

n = size(b,1);

% Use zero initial solution
x0 = zeros(n,1);
% Compute initial residual
r0 = b - prod_A(x0);

% Prepare for the iteration
v1 = r0 / beta(1);
x  = x0 + gamma(1)*v1;
v0 = zeros(n,1);
bet = 0.0;

% Perform iteration  [CG-like] 
for k=1:deg
    v  = prod_A(v1);
    v  = v - alpha(k)*v1 - bet*v0;
    v0 = v1;
    v1 = v / beta(k+1);
    x  = x + gamma(k+1)*v1;
    bet = beta(k+1);
end

end
