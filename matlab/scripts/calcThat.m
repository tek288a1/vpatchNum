%% loading data
clear;
addpath('../fnc')
n_in = 128;
rho_in = 0.40;
beta_opt = 0;
ddir = '~/Google Drive/Vpatch1/matlab1/data.fortran/';
dfile = sprintf('vp_b%1d_n%02d_r%4d.dat',...
                   beta_opt, log2(n_in), rho_in*1e4);
outdir = '~/Google Drive/Vpatch1/matlab1/data.matlab/';
outfile = sprintf('amat_b%1d_n%02d_r%4d.dat',...
                   beta_opt, log2(n_in), rho_in*1e4);

fort_dat = load( strcat(ddir, dfile) );

% inputs:
% n = size of quasi-solution; number of coefficients + 1
% rho
% beta
% a = quasi-solution
n = fort_dat(1);
rho = fort_dat(2);
beta = fort_dat(3);
u = fort_dat(4);
a = fort_dat(4:end);

% body of function begins here:
n4 = 4*n;
n2 = 2*n;
npt = n4;
nj = 100;
nu = linspace(0, 2*pi, npt+1)';
nu(end) = [];

eta = exp(1i*nu);
[z, zd, zeta, dthdnu] = ptval(rho, beta, a, eta);
[A0, q0] = A0valnew(z, zd, u);
jt = jtval(nj, rho, zeta);
c0 = c0val(A0, zd);
omega0 = q0 .* dthdnu .* (c0 - u*jt);
[~, A1] = j1val(eta, z, zd, omega0);
[~, A2] = j2val(eta, z, zd, omega0);
[~, A3] = j3val(eta, z, zd, omega0);
That = thval(A1, A2, A3);              % time consuming

