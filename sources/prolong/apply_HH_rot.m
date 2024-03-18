function x = apply_HH_rot(w_rot,y)
%------------------------------------------------------------------------
%                                                                       %
% Function to apply an Householder rotation to a vector: x = (I-2w*w')y %
%                                                                       %
%------------------------------------------------------------------------

x = y - 2*w_rot*(w_rot'*y);

return
