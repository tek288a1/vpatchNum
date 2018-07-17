rho = 0.4;
% beta = (1 - sqrt(1-rho^2))/rho;
beta = 0;
nj = 150;
n = 128;
nu = linspace(0, 2*pi, n+1)';
nu(end) = [];
eta = exp(1i*nu);
zeta = (eta - beta)./(1 - beta*eta);
jt = real( jtval(nj, rho, zeta) );

jt_maple = load('~/Dropbox/kim/jacobitheta.dat');

norm(jt - jt_maple, 'inf')
