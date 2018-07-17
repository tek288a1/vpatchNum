function jt = jtval(nj, rho, zeta)
%% j(nu) calculation
% input:
%   nj = number of terms to be summed
%   rho
%   zeta (not eta)
% output:
%   jt
%
% This m-file numerically evaluates j(\nu) = J( \zeta(\nu) ) at uniformly
% spaced-out \nu-points on [0, 2\pi).
%
%   J(\zeta) = 2 \zeta \sum_{k=1}^{\infty} \rho^{2k-1} 
%       ( 1/( \zeta + \rho^{2k-1} )^2 + 1/( 1 + \rho^{2k-1} \zeta )^2 )
%
% This function is related to Jacobi-Theta function. 
  jt = zeros(size(zeta));
  rhotmp = rho;
  for k = 1:nj
      jt = jt + rhotmp*(1./(zeta+rhotmp).^2 + 1./(1+rhotmp*zeta).^2);
      rhotmp = rhotmp*rho^2;
  end
  jt = 2*zeta.*jt;
end
