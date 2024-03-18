function [S,S_patt] = cpt_SoC_aff(SoC_perc,A,S_patt,TV,VERBOSE)
%-----------------------------------------------------------------------------------------
%
% Function to compute affinity-based strength of connection. This version also computes
% negative affinities to check that they are not taken into account.
%
% Input
%
% SoC_perc   desired percentage of non-zeroes to be retained in S relative to the
%            original non-zeroes of A
% A          system matrix
% S_patt     sparsity pattern containing every possible connection
% TV         test space used
%
% Output
%
% S       SoC matrix
% S_patt  SoC matrix with 1 indicating strong connections and zero elsewhere
%
%-----------------------------------------------------------------------------------------

global DEBINFO

S_patt = tril(S_patt,-1);
nn = size(S_patt,1);
nt = nnz(S_patt);
TV = TV';
% Compute row affinity
[ii,jj,aa] = find(S_patt);
ss = zeros(nt+1,1);
for ind = 1:nt

   if (VERBOSE > )2 && (mod(ind,floor(nt/10)) == 0)
      fprintf('Processed %9i entries out of %9i\n',ind,nt)
   end

   v = TV(:,ii(ind));
   w = TV(:,jj(ind));
   ss(ind) = cpt_sign_aff(v,w);

end
ii(nt+1) = nn; % This is for the diagonal entry
jj(nt+1) = nn; % This is for the diagonal entry
S = sparse(ii,jj,ss);

% Complete the matrix with a unitary diagonal
S = S + S' + speye(nn);

% Filter SoC sig retaining a given percentage of A off-diagonals
% Since S contains the diagonal, we retain nn more
if VERBOSE > 2
   fprintf('Filtering SoC sig\n')
end
nnz_eff = nnz(A) - nn;
ncut = min(nnz(S),floor(nnz_eff*SoC_perc/100)+nn);
[ii,jj,ss] = find(S);
[ss,perm] = sort(ss,'descend');
ii_f = ii(perm(1:ncut));
jj_f = jj(perm(1:ncut));
sf = ss(1:ncut);
S = sparse(ii_f,jj_f,sf);
if VERBOSE > 2
   fprintf('End Filtering SoC sig\n')
   fprintf('# of S off-diagonal entries: %10d\n',nnz(S)-nn);
   fprintf('\nS density over A: # %10.4f\n',nnz(S)/nnz(A));
end

%Create the pattern for SoC
ii = ii(perm);
jj = jj(perm);
ss(1:ncut) = 1;
ss(ncut+1:end) = -1;
S_patt = sparse(ii,jj,ss);
%save('SSS','S','S_patt');

if DEBINFO.coarsen.draw_dist
   % Draw SoC distribution
   figure; plot(ss(end:-1:1),'-b'); title('SoC distribution');
   figure; plot(sf(end:-1:1),'-r'); title('Filtered SoC distribution');
   drawnow
end

end
