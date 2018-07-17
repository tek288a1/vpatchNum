%% loading quasi-solution data: fortran and matlab routines
clear, clc
fort_dat = load('testdata.fortran.dat');
mat_dat = load('testdata.matlab.dat');
n = fort_dat(1);
rho = fort_dat(2);
beta = fort_dat(3);
u = fort_dat(4);
a = fort_dat(4:end);
u1 = mat_dat(4);
a1 = mat_dat(4:end);
n4 = 4*n;
n2 = 2*n;
npt = n4;
nj = 100;
nu = linspace(0, 2*pi, npt+1)';
nu(end) = [];
eta = exp(1i*nu);

%% loading data
clear; clc
n_in = 128;
rho_in = 0.3;
beta_opt = 1;
ddir = '~/Google Drive/Vpatch1/matlab1/data.fortran/';
dfile = sprintf('vp_b%1d_n%02d_r%4d.dat',...
                   beta_opt, log2(n_in), rho_in*1e4);
fort_dat = load( strcat(ddir, dfile) );
n = fort_dat(1);
rho = fort_dat(2);
beta = fort_dat(3);
u = fort_dat(4);
a = fort_dat(4:end);
n4 = 4*n;
n2 = 2*n;
npt = n4;
nj = 100;
nu = linspace(0, 2*pi, npt+1)';
nu(end) = [];
eta = exp(1i*nu);

%% calculation: matlab routines using fortran data
[z, zd, zeta, dthdnu] = ptval(rho, beta, a, eta);
M = trmatval(n, rho, beta);
jt = jtval(nj, rho, zeta);
[A0, q0] = A0valnew(z, zd, u);
c0 = c0val(A0, zd);
omega0 = q0 .* dthdnu .* (c0 - u*jt);
Lu = q0 .* (omega0 + jt.*dthdnu);
[J1, A1] = j1val(eta, z, zd, omega0);
[J2, A2] = j2val(eta, z, zd, omega0);
[J3, A3] = j3val(eta, z, zd, omega0);
[H4, B4] = h4val(z, zd);
T_hat = thval(A1, A2, A3);
T1 = t1val(T_hat, B4, dthdnu);
A1jm = A1jmval(q0, T1, M);
Ajm = Ajmval(q0, T1, M, Lu);

%% forming matrix by directly fft-ing point-value data
i_rng = (1:n/2);
A = real( Ajm(i_rng, i_rng) );
Ach = zeros(size(A));
for m = i_rng
    A1g = A1gval(m, rho, beta, z, zd);
    c11g = c11gval(m, rho, beta, A0, A1g, omega0);
    [T1g, T10g] = T1gval(m, q0, omega0, dthdnu, A1g, c11g, Lu);
    Tg = Tgval(T1g, Lu);
    Tgm = real(fft(Tg)/n4);
    Ach(:,m) = Tgm(1+i_rng);
end
%% For comparison, use modes from 1 to n/2
fprintf([' Max. error btw a_{j,m} computed via convolution and via ' ...
         'direct fft: %12.4e\n'], norm(A-Ach, 'inf'))

%% inspecting (super-/sub-)diagonal terms
imagesc(A);
axis image, colorbar
fprintf('\n');
fprintf(' mean of k-th diagonal terms:\n');
fprintf(' the 2nd super-diagonal terms are dominant.\n');
fprintf('\n');
fprintf(' %5s %13s \n', 'k', 'mean');
fprintf(' %s\n', repmat('-', 1, 19));
for k = 4+(-10:10)
    fprintf(' %5d %13.4g\n', k, mean(diag(A, k)));
end

%% L-inverse and an approximate inverse
nb2 = length(i_rng);
Linv = inv(-diag(i_rng)+A);
sh = 1;
d(0+sh) = mean(diag(A, 0));
d(1+sh) = mean(diag(A, 1));
d(2+sh) = mean(diag(A, 2));

I = eye(nb2);
D = diag( -i_rng+d(0+sh) );
Dinv = inv(D);
N1 = d(1+sh)*Dinv*diag(ones(nb2-1,1), 1);
N2 = d(2+sh)*Dinv*diag(ones(nb2-2,1), 2);
% Lsinv = (I - N1 - N2 + N1^2 + N2^2 + N1*N2 + N2*N1)*Dinv;
Lsinv = (I - N1 - N2)*Dinv;

norm(Linv - Lsinv)



%% manipulation
% imagesc(real(Ajm));
% axis image, colorbar
imagesc(Linv(1:10,1:10));
axis image, colorbar
for d = -5:5
    fprintf('mean of k-th diagonal terms: k = %3d, mean = %11.4g\n', ...
            d, mean(diag(Linv, d)));
end



%% check if coefficients are real
jt_fc = fft(jt)/n4;
A0_fc = fft(A0)/n4;
q0_fc = fft(q0)/n4;
omega0_fc = fft(omega0)/n4;
Lu_fc = fft(Lu)/n4;
dthdnu_fc = fft(dthdnu)/n4;

fprintf('\n')
fprintf('%s\n', repmat('-', 1, 80))
fprintf('%s\n', ... 
        'Checking real-ness of coefficients and matrix entries computed.')
fprintf('%s\n', ...
        ' Following are max-norm of imaginary parts of the results of FFT.')
fprintf('%s %4d \n', ' : N used for FFT: ', n4)
fprintf('%s %4d%s%4d \n', ...
        ' : Dimensions of A1, A2, A3, B4: ', n2, 'x', n2)
fprintf('%s %4d%s%4d \n', ...
        ' : Dimensions of T_hat and T1: ', n2, 'x', n2)
fprintf('%s %4d%s%4d \n', ...
        ' : Dimensions of A1jm: ', n, 'x', n)
fprintf('%s %4d%s%4d \n', ...
        ' : Dimensions of Ajm: ', n-1, 'x', n-1)
fprintf('\n')
fprintf('%16s %12.4e\n', 'j(nu):', norm( imag(jt_fc), 'inf' ))
fprintf('%16s %12.4e\n', 'A0(nu):', norm( imag(A0_fc), 'inf' ))
fprintf('%16s %12.4e\n', 'q0(nu):', norm( imag(q0_fc), 'inf' ))
fprintf('%16s %12.4e\n', 'omega0(nu):', norm( imag(omega0_fc), 'inf' ))
fprintf('%16s %12.4e\n', 'Lu(nu):', norm( imag(Lu_fc), 'inf' ))
fprintf('%16s %12.4e\n', 'A1_{j,m}:', norm( imag(A1(:)), 'inf' ))
fprintf('%16s %12.4e\n', 'A2_{j,m}:', norm( imag(A2(:)), 'inf' ))
fprintf('%16s %12.4e\n', 'A3_{j,m}:', norm( imag(A3(:)), 'inf' ))
fprintf('%16s %12.4e\n', 'B4_{j,m}:', norm( imag(B4(:)), 'inf' ))
fprintf('%16s %12.4e\n', 't_hat_{j,m}:', norm( imag(T_hat(:)), 'inf' ))
fprintf('%16s %12.4e\n', 't1_{j,m}:', norm( imag(T1(:)), 'inf' ))
fprintf('%16s %12.4e\n', 'a1_{j,m}:', norm( imag(A1jm(:)), 'inf' ))
fprintf('%16s %12.4e\n', 'a_{j,m}:', norm( imag(Ajm(:)), 'inf' ))

%% calculation: matlab routines using matlab data
[z1, zd1, zeta1, dthdnu1] = ptval(rho, beta, a1, eta);
jt1 = jtval(nj, rho, zeta1);
[A01, q01] = A0valnew(z1, zd1, u1);
c01 = c0val(A01, zd1);
omega01 = q01 .* dthdnu1 .* (c01 - u1*jt1);
Lu1 = q01 .* (omega01 + jt1.*dthdnu1);
[J11, A11] = j1val(eta, z1, zd1, omega01);
[J21, A21] = j2val(eta, z1, zd1, omega01);
[J31, A31] = j3val(eta, z1, zd1, omega01);
[H41, B41] = h4val(z1, zd1);
T_hat1 = thval(A11, A21, A31);
T11 = t1val(T_hat1, B41, dthdnu1);
A1jm1 = A1jmval(q01, T11, M);
Ajm1 = Ajmval(q01, T11, M, Lu1);

%% loading data for comparison: same quasi-solutions, different routines
base = '~/Google Drive/VP_fortran/tmp/';
z2 = cplxload([base, 'z.dat']);
zd2 = cplxload([base, 'zd.dat']);
zeta2 = cplxload([base, 'zeta.dat']);
dthdnu2 = cplxload([base, 'dthdnu.dat']);
jt2 = cplxload([base, 'jt.dat']);
A02 = cplxload([base, 'A0.dat']);
q02 = cplxload([base, 'q0.dat']);
c02 = cplxload([base, 'c0.dat']);
omega02 = cplxload([base, 'omega0.dat']);
J12 = cmatload([base, 'J1.dat']);
J22 = cmatload([base, 'J2.dat']);
J32 = cmatload([base, 'J3.dat']);
H42 = cmatload([base, 'H4.dat']);
A12 = cmatload([base, 'A1.dat']);
A22 = cmatload([base, 'A2.dat']);
A32 = cmatload([base, 'A3.dat']);
B42 = cmatload([base, 'B4.dat']);
T_hat2 = cmatload([base, 'T_hat.dat']);
T12 = cmatload([base, 'T1.dat']);
A1jm2 = cmatload([base, 'A1jm.dat']);
Ajm2 = cmatload([base, 'Ajm.dat']);

%% compare: different quasi-solutions used in matlab routines
% _ : fortran data; matlab routines
% _1: matlab data; matlab routines
% _2: fortran data; fortran routines
fprintf('\n')
fprintf('%s\n', 'The difference between values calculated using')
fprintf('%s\n', '  - different quasi-solution data (matlab vs fortran)')
fprintf('%s\n', '  - same matlab routines')
fprintf('\n')
fprintf('%16s %12.4e\n', 'quasi-solution:', norm(a-a1, 'inf'))
fprintf('%16s %12.4e\n', 'f0(nu):', norm(z-z1, 'inf'))
fprintf('%16s %12.4e\n', 'f0p(nu):', norm(zd-zd1, 'inf'))
fprintf('%16s %12.4e\n', 'zeta(nu):', norm(zeta-zeta1, 'inf'))
fprintf('%16s %12.4e\n', 'dth_dnu:', norm(dthdnu-dthdnu1, 'inf'))
fprintf('%16s %12.4e\n', 'j(nu):', norm(jt-jt1, 'inf'))
fprintf('%16s %12.4e\n', 'A0(nu):', norm(A0-A01, 'inf'))
fprintf('%16s %12.4e\n', 'q0(nu):', norm(q0-q01, 'inf'))
fprintf('%16s %12.4e\n', 'c0:', norm(c0-c01, 'inf'))
fprintf('%16s %12.4e\n', 'omega0(nu):', norm(omega0-omega01, 'inf'))
fprintf('%16s %12.4e\n', 'j1(nu,th):', norm(J1-J11, 'inf')) 
fprintf('%16s %12.4e\n', 'j2(nu,th):', norm(J2-J21, 'inf')) 
fprintf('%16s %12.4e\n', 'j3(nu,th):', norm(J3-J31, 'inf')) 
fprintf('%16s %12.4e\n', 'h4(nu,th):', norm(H4-H41, 'inf'))
fprintf('%16s %12.4e\n', 'A1_{j,k}:', norm(A1-A11, 'inf'))
fprintf('%16s %12.4e\n', 'A2_{j,k}:', norm(A2-A21, 'inf'))
fprintf('%16s %12.4e\n', 'A3_{j,k}:', norm(A3-A31, 'inf'))
fprintf('%16s %12.4e\n', 'B4_{j,k}:', norm(B4-B41, 'inf'))
fprintf('%16s %12.4e\n', 't_hat_{j,k}:', norm(T_hat-T_hat1, 'inf'))
fprintf('%16s %12.4e\n', 't1_{j,k}:', norm(T1-T11, 'inf'))
fprintf('%16s %12.4e\n', 'a1_{j,k}:', norm(A1jm-A1jm1, 'inf'))
fprintf('%16s %12.4e\n', 'a_{j,k}:', norm(Ajm-Ajm1, 'inf'))

%%
fprintf('\n')
fprintf('%s\n', 'The difference between values calculated using ')
fprintf('%s\n', '  - same quasi-solution data                   ')
fprintf('%s\n', '  - different routines (matlab vs fortran)     ')
fprintf('\n')
fprintf('%16s %12s\n', 'quasi-solution:', 'N/A')
fprintf('%16s %12.4e\n', 'f0(nu):', norm(z-z2, 'inf'))
fprintf('%16s %12.4e\n', 'f0p(nu):', norm(zd-zd2, 'inf'))
fprintf('%16s %12.4e\n', 'zeta(nu):', norm(zeta-zeta2, 'inf'))
fprintf('%16s %12.4e\n', 'dth_dnu:', norm(dthdnu-dthdnu2, 'inf'))
fprintf('%16s %12.4e\n', 'j(nu):', norm(jt-jt2, 'inf'))
fprintf('%16s %12.4e\n', 'A0(nu):', norm(A0-A02, 'inf'))
fprintf('%16s %12.4e\n', 'q0(nu):', norm(q0-q02, 'inf'))
fprintf('%16s %12.4e\n', 'c0:', norm(c0-c02, 'inf'))
fprintf('%16s %12.4e\n', 'omega0(nu):', norm(omega0-omega02, 'inf'))
fprintf('%16s %12.4e\n', 'j1(nu,th):', norm(J1-J12, 'inf'))
fprintf('%16s %12.4e\n', 'j2(nu,th):', norm(J2-J22, 'inf')) 
fprintf('%16s %12.4e\n', 'j3(nu,th):', norm(J3-J32, 'inf')) 
fprintf('%16s %12.4e\n', 'h4(nu,th):', norm(H4-H42, 'inf'))
fprintf('%16s %12.4e\n', 'A1_{j,k}:', norm(A1-A12, 'inf'))
fprintf('%16s %12.4e\n', 'A2_{j,k}:', norm(A2-A22, 'inf'))
fprintf('%16s %12.4e\n', 'A3_{j,k}:', norm(A3-A32, 'inf'))
fprintf('%16s %12.4e\n', 'B4_{j,k}:', norm(B4-B42, 'inf'))
fprintf('%16s %12.4e\n', 't_hat_{j,k}:', norm(T_hat-T_hat2, 'inf'))
fprintf('%16s %12.4e\n', 't1_{j,k}:', norm(T1-T12, 'inf'))
fprintf('%16s %12.4e\n', 'a1_{j,k}:', norm(A1jm-A1jm2, 'inf'))
fprintf('%16s %12.4e\n', 'a_{j,k}:', norm(Ajm-Ajm2, 'inf'))

%%
addpath('../fnc')
a_tanveer = load('../data.tanveer/tmp.txt');
a_matlab = load('../data.matlab/vp_b1_n07_r8500.matlab.dat');
a_fortran = load('../data.fortran/vp_b1_n07_r8500.dat');
a_tanveer(129:end) = [];
a_matlab(1:3) = [];
a_fortran(1:3) = [];

norm(a_matlab-a_fortran)
norm(a_matlab-a_tanveer)
norm(a_fortran-a_tanveer)
rmpath('../fnc')

%% Matrix A comparison: fortran vs. matlab
% 07/04/18

clear
addpath('../functions');
n = 256;

A = cmatload('../../fortran90/tmp/Ajm.dat'); % fortran code
                                             % produces 255 x 255
A = real(A);
Afor = A(1:n/2, 1:n/2);

fname = '../../data/data.matlab/amat_b1_n08_r8500.bin';
fid = fopen(fname);
Amat = fread(fid, [n/2, n/2], 'double');
fclose(fid);


norm(Afor-Amat, 'inf')

%% Matrix A comparison: fortran vs. matlab
% 07/05/18

clear
addpath('../functions');
n = 1024;

A = ...                                 % fortran 
    cmatload('../../data/data.fortran/amat_b1_n10_r9000.dat');
%A = real(A);
%Afor = A(1:n/2, 1:n/2);

%%
fname = '../../data/data.matlab/amat_b1_n08_r8500.bin';
fid = fopen(fname);
Amat = fread(fid, [n/2, n/2], 'double');
fclose(fid);


norm(Afor-Amat, 'inf')
