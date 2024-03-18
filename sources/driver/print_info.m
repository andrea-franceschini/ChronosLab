function print_info(AMG_prec,A,iter,relres,resreal)

%-----------------------------------------------------------------------------------------

fprintf('\n');
fprintf('------------------------------------------------------------------------------------------');
fprintf('\n');

% Get AMG hierarchy information
AMG_info = get_AMG_info(AMG_prec,A);

% Print AMG hierarchy information
print_AMG_info(AMG_info);

fprintf('\n');
fprintf('Grid Complexity:     %10.4f\n',AMG_info(1).grid_comp);
fprintf('Operator Complexity: %10.4f\n',AMG_info(1).oper_comp);

% Print results
if flag
   fprintf('Convergence not achieved\n');
end
fprintf('\n');
fprintf('\n');
fprintf('# of iterations:   %15d\n',iter);
fprintf('Relative residual: %15.6e\n',relres);
fprintf('Real residual:     %15.6e\n',resreal);

fprintf('\n');
fprintf('Smoother time:     %10.2f\n',T_smoo);
fprintf('Tspace time:       %10.2f\n',T_tspa);
fprintf('Coarsening time:   %10.2f\n',T_coar);
fprintf('Comp. Relax. time: %10.2f\n',T_CR);
fprintf('Prolongation time: %10.2f\n',T_prol);
fprintf('Iteration time:    %10.2f\n',T_iter);

end
