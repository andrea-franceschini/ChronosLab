function rho = comp_relax(CR_iters,fcnode,smoother,A)

ind_F = find(fcnode < 0);
As = A(ind_F,ind_F);

%S'*(L*L')*S;
% CR = I - inv(S'*M*S)*(S'*A*S)

rho = 1;

end
