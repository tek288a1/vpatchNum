%% script m-file: double fft routine
clear, clc
vp = load('vp_r05to90_b0_n128.mat');
rho = 0.5;
id = find(vp.RHO == rho);
a = vp.A(:, id);
beta = vp.BETA(id);
n = vp.nr_term
%-------------
nrpt = 4*n;
n2 = 2*n;
nu = linspace(0, 2*pi, nrpt+1)';
nu(end) = [];
eta = exp(1i*nu);
[z, zd, zeta] = ptval(rho, beta, a, eta);
etah = (rho^2 + beta*zeta)./(zeta + beta*rho^2);

%%
inu = 1:2:nrpt;
ith = 2:2:nrpt;
f0nu = z(inu);
f0th = z(ith).';                         % transpose, not hermitian
f0pth = zd(ith).';
kern = f0pth.*( abs(f0nu).^2 - abs(f0th).^2 )./( f0nu.^2 - f0th.^2 );

%%
f0nu_comp = load('~/Google Drive/VP_fortran/src/output.dat');
f0 = [f0nu_comp(:,1) + 1i*f0nu_comp(:,2)]
norm(real(f0nu)-f0nu_comp(:,1))
norm(imag(f0nu)-f0nu_comp(:,2))
norm(f0nu-f0)


%% 
data = load('~/Google Drive/VP_fortran/data/vp_b0_n07_r4000.dat');
adata = data(4:end);
[z1, zd1] = ptval(rho, beta, adata, eta);
f0odd = z(inu);


%% shifted fft
fc = fft(f0nu)/n;
fcs = fft(f0th)/n .* exp(-1i*pi*(0:n-1)/n); % this is wrong; see below
norm(fc-fcs.')

%% double fft
% fft across rows (fixed nu, varying theta): required phase-shift
ft_row = fft(kern, n, 2)/n;
phase_factor = exp(-1i*pi*(0:n-1)/n);
ft_row = ft_row.*phase_factor;
% fft across columns (fixed theta, varying nu)
ft = fft(ft_row, n, 1)/n;

%% alternately, we can do fft across columns first then rows.
ft_col = fft(kern, n, 1)/n;
ft_comp = fft(ft_col, n, 2)/n;
ft_comp = ft_comp.*phase_factor;
norm(ft-ft_comp)

%% check with a0
a0ft = n*ifft(ft(:, 1));
a0int = a0val(f0nu, f0th, f0pth, a(1));
norm(a0ft-a0int)

%% import solution and basic pt evaluation
clear, clc
sh = 1;
vp = load('vp_r05to90_b1_n128.mat');
rho = 0.4;
id = find(vp.RHO == rho);
a = vp.A(:, id);
beta = vp.BETA(id);
u = a(1);
n = vp.nr_term;
n4 = 4*n;
n2 = 2*n;
npt = n2;
nj = 100;
nu = linspace(0, 2*pi, n2+1)';
nu(end) = [];
eta = exp(1i*nu);
[z, zd, zeta, dthdnu] = ptval(rho, beta, a, eta);
M = trmatval(n, rho, beta);
jt = jtval(nj, rho, zeta);
[a0, q0] = a0valnew(z, zd, u);
c0 = c0val(a0, zd);
omega0 = q0 .* dthdnu .* (c0 - u*jt);
[H, B] = h4val(z, zd);
% [H, B] = h4valnew(z, zd);
[J1, A1] = j1val(eta, z, zd, omega0);

%% check
%%%% loading comparison data
base = '~/Google Drive/VP_fortran/tmp/';
% M1 = load('~/Google Drive/VP_fortran/tmp/mdata.dat');
% a1 = load('~/Google Drive/VP_fortran/tmp/adata.dat');
% jt1 = load('~/Google Drive/VP_fortran/tmp/jval.dat');
% jt1 = complex(jt1(:,1), jt1(:,2));
% zeta1 = load('~/Google Drive/VP_fortran/tmp/zeta.dat');
% zeta1 = complex( zeta1(:,1), zeta1(:,2) );
% a2 = load('~/Google Drive/VP_fortran/tmp/adataST.dat');
% a2(1:3) = [];
% H1 = load([base, 'h4.dat']);
% H1 = complex(H1(:,1), H1(:,2));
% H1 = reshape(H1, npt/2, npt/2).';
B1 = load([base, 'b4.dat']);
B1 = complex(B1(:,1), B1(:,2));
B1 = reshape(B1, npt/2, npt/2).'; 
B2 = load([base, 'b41.dat']);
B2 = complex(B2(:,1), B2(:,2));
B2 = reshape(B2, npt/2, npt/2).'; 
% a01 = n*ifft(B(:,1));

%%%% comparison
% norm(M-M1)
% norm(a-a1)
% norm(a1-a2)
% norm(jt-jt1)
% norm(zeta-zeta1)
% norm(H-H1)
% norm(B-B1)
norm(B1-B2)
% norm(a0-a01)

%% shifted fft
clear, clc
n = 32;
th = linspace(0, 2*pi, n+1)';
th(end) = [];
z = exp(1i*th);
x = sum( (z/2).^(1:n/2-1), 2 ) + 1 + sum( (1./(2*z)).^(1:n/2-1), 2 );
y = fft(x)/n;

thsh = pi/n + th;
zsh = exp(1i*thsh);
xsh = sum( (zsh/2).^(1:n/2-1), 2 ) + 1 + sum( (1./(2*zsh)).^(1:n/2-1), 2 );
% ph_adj = exp(-1i*pi*(0:n-1)/n).';
% ph_adj(n/2+2:end) = -ph_adj(n/2+2:end);
% ph_adj(n/2+1) = 0;
v = fftshift( -n/2+(0:n-1)' )/n;
% v = (0:n-1)'/n; % this possibly is the mistake/error in Tanveer's code
ph_adj = exp(-1i*pi*v);
ys = fft(xsh)/n .* ph_adj;
ys(n/2+1) = 0;
norm(y - ys)

% recovery
norm( x - n*ifft(y) )
norm( x - n*ifft(ys) )
norm( xsh - n*ifft(ys./ ph_adj) )

%% fftshift
clear, clc
n = 8;
block = ones(n/2, n/2);
X = [ block, 2*block; 3*block, 4*block ]
Xsh = fftshift(X)
