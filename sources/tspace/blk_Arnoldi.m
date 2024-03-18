function [vec, val, iter, resid] = blk_Arnoldi(V0,ProdMat,blk_size,m,ntv,dual_orth)

n = size(V0, 1);
V = zeros(n,blk_size*m);
H = zeros(blk_size*m);
W = zeros(n,blk_size);
FAC = zeros(blk_size);
[V0,~] = qr(V0,0);
V(:,1:blk_size) = V0(:,1:blk_size);

istart_n = 1;
iend_n = blk_size;
for i = 1:m

   istart = istart_n;
   iend = iend_n;
   W = ProdMat(V(:,istart:iend));
   % First orthogonalization
   jend = 0;
   for j = 1:i
      jstart = jend + 1;
      jend = jstart + blk_size - 1;
      H(jstart:jend,istart:iend) = V(:,jstart:jend)'*W;
      W = W - V(:,jstart:jend)*H(jstart:jend,istart:iend);
   end
   % Second orthogonalization
   if dual_orth
      jend = 0;
      for j = 1:i
         jstart = jend + 1;
         jend = jstart + blk_size - 1;
         FAC = V(:,jstart:jend)'*W;
         H(jstart:jend,istart:iend) = H(jstart:jend,istart:iend) + FAC;
         W = W - V(:,jstart:jend)*FAC;
      end
   end
   if i == m; break; end
   istart_n = iend + 1;
   iend_n = iend + blk_size;
   [V(:,istart_n:iend_n),H(istart_n:iend_n,istart:iend)] = qr(W,0);

end

% Compute eigenpairs of H
[Vr,lambda] = eig(H);
val = diag(lambda);
vec = V*Vr;

% Estimate residual
res = ProdMat(vec)-vec*diag(val);
lres = zeros(m,1);
for i = 1:m
   lres(i) = norm(res(:,i));
end
lres_rel = lres./val;
resid = lres_rel;

% Sort according to minimal val
[val, perm] = sort(val, 'ascend');
vec = vec(:, perm);
resid = resid(perm);

vec = vec(:, 1:ntv);
val = val(1:ntv);
resid = [nan(ntv,1), resid(1:ntv)];

% Just for compatibility
iter = m;

end
