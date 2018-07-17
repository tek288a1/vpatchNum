%% loading data
clear;
added_path = ['../functions']; 
addpath(added_path);

% n_in = input(' Input N = ');
% rho_in = input(' Input rho = ');
% beta_opt = input(' Input beta_opt = ');
% kk = input(' Input K = ');
n_in = 256;
rho_in = 0.92;
beta_opt = 1;
ddir = '../../data/data.fortran/';
dfile = sprintf('vp_b%1d_n%02d_r%4d.dat', ...
                beta_opt, log2(n_in), rho_in*1e4);
data = load( strcat(ddir, dfile) );
n = data(1);
rho = data(2);
beta = data(3);
u = data(4);
a = data(4:end);
n4 = 4*n;
n2 = 2*n;
npt = n4;
nj = 100;
nu = linspace(0, 2*pi, npt+1)';
nu(end) = [];

% inner and outer 
r_in = rho^(1/8);
r_out = 1/r_in;
eta = exp(1i*nu);
eta_in = r_in*eta;
eta_out = r_out*eta;

%% calculation of point values
M = trmatval(n, rho, beta);

% ind_trunc = 1:100;
n_coef = length(a)-1;
ind_trunc = 1:length(a);
[z, zd, zeta, dthdnu] = ptval(rho, beta, a(ind_trunc), eta);
[z_in, zd_in, zeta_in, dthdnu_in] = ptval(rho, beta, a(ind_trunc), eta_in);
[z_out, zd_out, zeta_out, dthdnu_out] = ptval(rho, beta, a(ind_trunc), eta_out);

jt = jtval(nj, rho, zeta);
jt_in = jtval(nj, rho, zeta_in);
jt_out = jtval(nj, rho, zeta_out);

[A0, q0] = a0valnew(z, zd, u);
[A0_in, q0_in] = a0val(z_in, z, zd, u); % careful with quadrature code
[A0_out, q0_out] = a0val(z_out, z, zd, u); % same here

c0 = c0val(A0, zd);
c0_in = c0val(A0_in, zd_in);
c0_out = c0val(A0_out, zd_out);

omega0 = q0 .* dthdnu .* (c0 - u*jt);
omega0_in = q0_in .* dthdnu_in .* (c0_in - u*jt_in);
omega0_out = q0_out .* dthdnu_out .* (c0_out - u*jt_out);

Lu = q0 .* (omega0 + jt.*dthdnu);
Lu_in = q0_in .* (omega0_in + jt_in.*dthdnu_in);
Lu_out = q0_out .* (omega0_out + jt_out.*dthdnu_out);

%% check: image of the eta-annulus be contained in the zeta-annulus
clf
hold off
hold on
% plot(real(z), imag(z), real(z_in), imag(z_in), real(z_out), ...
%      imag(z_out))
plot(real(zeta), imag(zeta), 'b')
plot(real(zeta_in), imag(zeta_in), 'r')
plot(real(zeta_out), imag(zeta_out), 'r')
plot(real(rho*zeta), imag(rho*zeta), 'b')
plot(real(1/rho*zeta), imag(1/rho*zeta), 'b')
axis equal

%% upper bounds on divided sums/differences
% Upper bounds of four functions Ql(nu, theta), l=1..4 are examined.
% In the following loops, 4 "boundaries" of S^2 are examined.
% k = 1 : (i,j) = (1,1) : |eta| = r_in, |zeta| = r_in
% k = 2 : (i,j) = (1,2) : |eta| = r_in, |zeta| = r_out
% k = 3 : (i,j) = (1,1) : |eta| = r_out, |zeta| = r_in
% k = 4 : (i,j) = (1,1) : |eta| = r_out, |zeta| = r_out
% 
% The maximal values are recorded on Q_tmp, rows for different Ql and
% columns for different k.
Q = zeros(n, n, 4);
Q_tmp = zeros(4, 4);
k = 0;
for i = 1:2
    if i == 1
        eiNu = r_in*repmat(eta(1:2:n2), 1, n);
        Znu = repmat(z_in(1:2:n2), 1, n);
    else
        eiNu = r_out*repmat(eta(1:2:n2), 1, n);
        Znu = repmat(z_out(1:2:n2), 1, n);
    end
    for j = 1:2
        if j == 1
            eiTh = r_in*repmat(eta(2:2:n2).', n, 1);
            Zth = repmat(z_in(2:2:n2).', n, 1);
        else
            eiTh = r_out*repmat(eta(2:2:n2).', n, 1);
            Zth = repmat(z_out(2:2:n2).', n, 1);
        end
        k = k + 1;
        Q(:,:,1) = (eiNu - eiTh)./(Znu.^2 - Zth.^2);
        Q(:,:,2) = conj(eiNu - eiTh)./(Znu.^2 - Zth.^2);
        Q(:,:,3) = 1./(Znu + Zth);
        Q(:,:,4) = conj(Znu-Zth)./(Znu-Zth) + ...
            conj(Znu+Zth)./(Znu+ Zth);
        Q_tmp(:, k) = reshape( max( max((abs(Q))) ), 4, 1 );
    end
end
Q_max = max(Q_tmp, [], 2)

%% 2*q0*exp(i*nu)
2*norm([r_in*q0_in; r_out*q0_out], inf)

%% q0*omega0
norm([q0_in.*omega0_in; q0_out.*omega0_out], inf)
% norm([q0_in; q0_out], inf)
% norm([omega0_in; omega0_out], inf)

%% q0*omega0*f0/(2*exp(i*nu))
1/2*norm([q0_in.*omega0_in.*z_in/r_in,...
          q0_out.*omega0_out.*z_out/r_out], inf)

%% omega0*f0'
norm([omega0_in.*zd_in; omega0.*zd_out], inf)

%% f0'
norm([z_in; z_out], inf)

%%
norm([q0_in; q0_out], inf)*norm([omega0_in; omega0_out], inf)
norm([q0_in.*omega0_in; q0_out.*omega0_out], inf)

%% Remove path (at end of script/script clean-up)
rmpath(added_path);
