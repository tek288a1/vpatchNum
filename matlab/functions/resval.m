function res = resval(rho, beta, a)
%RESVAL res = resval(rho, beta, a)
%
%  Evaluation of residual Im[ (A0(nu_j) + U0)*f0'(nu_j) ]
%  at n evenly-spaced out off-axis points in the upper-half
%  eta-plane. 
%
%  The integral A0(nu) is calculated at n collocation points
%
%	nu_j = pi/(2*n) + (j-1)*pi/n, j=1..n
%
%  with 2*n quadrature points
%
%	nu'_k = 2*pi*(k-1)/(2*n), k=1..2*n
%
%  by alternating-point trapezoidal method. 
%
%   Written by Tae Eun Kim, 12/03/2017
    a = a(:);
    u = a(1);
    res = zeros(size(a));
    n = length(a);
    n2 = 2*n;
    n4 = 4*n;
    ic = 2:2:n2;                        % collocation indices
    iq = 1:2:n4;                        % quadrature indices
    nu = linspace(0, 2*pi, n4+1)';
    nu(end) = [];
    eta = exp(1i*nu);
    [z, zd] = ptval(rho, beta, a, eta);
    a0 = a0val(z(ic), z(iq), zd(iq), u);
    res = imag( (a0 + u).*zd(ic) );
end
