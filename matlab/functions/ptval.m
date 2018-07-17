function [z, zd, zeta, dth_dnu] = ptval(rho, beta, a, eta)
%PTVAL	Point evaluation of F, F_{nu}, and zeta(eta) at given eta points
%
%   [z, zd, zeta] = ptval(rho, beta, a, eta) calculates, using
%   Horner's method
%
%	z = 1i*[ Fm_hat(eta) + alpha_hat(eta) - alpha_hat(eta_hat) ]
%
%   and its derivative with respect to the angular variable nu
%
%	zd = -eta*[ Fm_hat'(eta) + alpha_hat'(eta) 
%       	     - alpha_hat'(eta_hat)*eta_hat'(eta) ]
%
%   and 
%
%	zeta = (eta - beta)/(1 - beta*eta) .
%
%%%% Useful expressions
%
%  eta_hat = eta(rho^2/zeta) = (rho^2+beta*zeta)/(zeta+beta*rho^2)
%
%  dth_dnu = (1-beta^2)*eta/( (eta-beta)*(1-beta*eta) )
%
%  dnuh_dth = -( rho^2*(1-beta^2)*zeta )
%            /( (zeta+beta*rho)*(rho^2+beta*zeta) )
%
%  Fm_hat(eta) = (zeta-rho)/(zeta+rho)
%
%  eta*Fm_hat(eta) = 2*rho*zeta/(zeta+rho)^2*dth_dnu
%
%  alpha_hat(eta) = sum_{k=1}^{N} a_k * eta^k
%
%  eta*alpha_hat'(eta) = sum_{k=1}^{N} k * a_k * eta^k
%
%  alpha_hat(eta_hat) = sum_{k=1}^{N} a_k * eta_hat^k
%
%  eta*alpha_hat'(eta_hat)*eta_hat'(eta) 
%    = ( sum_{k=1}^{N} k * a_k * eta_hat^k ) * dnuh_dth * dth_dnu
%
%  INPUTs: a = [U, a_1, a_2, ... , a_{N-1}]' is a column vector of 
%             length N consisting of translational velocity and
%             series coefficients
%          rho = conformal parameter
%          beta = Mobius parameter
%          eta = equispaced points on the circle |eta| = r.
%
% Written by Tae Eun Kim, 12/03/2017
    n = length(a);
    nm1 = n-1;
    sh = 1;
    zeta = (eta - beta)./(1 - beta*eta);
    etah = (rho^2 + beta*zeta)./(zeta + beta*rho^2);
    dth_dnu = (1-beta^2)*eta./( (eta-beta).*(1-beta*eta) );
    dnuh_dth = -rho^2*(1-beta^2)*zeta...
        ./( (zeta+beta*rho^2).*(beta*zeta+rho^2) );
    fm = (zeta - rho)./(zeta + rho);    % Mobius
    fmd = 2*rho*zeta./(zeta + rho).^2.*dth_dnu;
    al = a(nm1+sh);                     % alpha(eta)
    ah = a(nm1+sh);                     % alpha_hat(eta)
    ald = nm1*al;
    ahd = nm1*ah;
    a(0+sh) = 0;
    for j = nm1-1:-1:0
        al = al.*eta + a(j+sh);
        ah = ah.*etah + a(j+sh);
        ald = ald.*eta + j*a(j+sh);
        ahd = ahd.*etah + j*a(j+sh);
    end
    z = 1i*(fm + al - ah);
    zd = -(fmd + ald - ahd.*dnuh_dth.*dth_dnu);
end
