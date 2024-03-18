function [fcnode, clist, flist] = cpt_MIS(S, fcnode, VERBOSE)

%  0 --> not assigned
% -1 --> FINE
%  1 --> COARSE

nn = size(S,1);
list = (1:nn)';

tic;
for i = 1 : nn

    % Pick up a random entry from the list
    k = i + floor(rand(1)*(nn+1-i));
    inod = list(k);
    list(k) = list(i);

    % Check whether it has been assigned yet
    if (fcnode(inod) == 0)
        % It is free, mark its neighbours as fine
        [jlist,~,~] = find(S(:,inod));
        % Dirichlet nodes
        ID = (fcnode(jlist) == -2);
        % Mark them as fine
        fcnode(jlist) = -1;
        % Restore Dirichlet nodes
        fcnode(jlist(ID)) = -2;
        % Mark it as coarse
        fcnode(inod) = 1;
    end
end
time = toc;
if VERBOSE > 2
   fprintf('Time for MIS: %f\n',time);
end

% Count and number coarse nodes
nc = 0;
nf = 0;
clist = zeros(nn,1);
flist = zeros(nn,1);
for i = 1:nn
   if fcnode(i) > 0
      nc = nc + 1;
      fcnode(i) = nc;
      clist(nc) = i;
   elseif fcnode(i) < 0
      nf = nf + 1;
      fcnode(i) = -nf;
      flist(nf) = i;
   end
end
clist = clist(1:nc);
flist = flist(1:nf);

end
