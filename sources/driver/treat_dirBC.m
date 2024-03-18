function A = treat_dirBC(A)
lmax = eigs(A,1,'lm','FailureTreatment','keep','Display',0,'Tolerance',0.001,'MaxIterations',3);
A(A==1) = lmax/10;
end
