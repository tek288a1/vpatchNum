function c1g = c1gval(m, A0, A1g, omega0)
% C1GVAL - c1g = c1gval(m, A0, A1g, omega0)
%   evaluates the functional c_1[g_m] which is defined to be the average
%   
%   < g'*A0 > + < omega0*A1[g] >
%   
    n = length(A0);
    nu = linspace(0, 2*pi, n+1)';
    nu(end) = [];
    g = 1i*exp(1i*m*nu);
    gd = -m*exp(1i*m*nu);
    c1g = 1/n*sum( A0.*gd + omega0.*A1g );
end
