function c11g = c11gval(m, rho, beta, A0, A1g, omega0)
% C11GVAL - c11g = c1gval(m, A0, A1g, omega0)
%   evaluates the functional c_1[g_m] which is defined to be the average
%   
%   < g'*A0 > + < omega0*A1[g] >
%   
    n = length(A0);
    nu = linspace(0, 2*pi, n+1)';
    nu(end) = [];
    eta = exp(1i*nu);
    zeta = (eta - beta)./(1 - beta*eta);
    etah = (rho^2 + beta*zeta)./(zeta + beta*rho^2);
    dzeta_deta = (1 - beta^2)./(1 - beta*eta).^2;
    % detah_dzeta = (beta*(zeta+beta*rho^2)-(rho^2+beta*zeta))./(zeta+beta*rho^2).^2;
    detah_dzeta = rho^2*(beta^2 - 1)./(zeta + beta*rho^2).^2;
    g = 1i*( eta.^m - etah.^m );
    gd = -m*( eta.^m - eta.*etah.^(m-1).*detah_dzeta.*dzeta_deta );    
    % g = 1i*exp(1i*m*nu);
    % gd = -m*exp(1i*m*nu);
    c11g = 1/n*sum( A0.*gd + omega0.*A1g );
end
