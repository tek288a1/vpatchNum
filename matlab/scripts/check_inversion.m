%% loading data
clear;
% n_in = input(' Input N = ');
% rho_in = input(' Input rho = ');
% beta_opt = input(' Input beta_opt = ');
% kk = input(' Input K = ');
n_in = 256;
rho_in = 0.91;
beta_opt = 1;
kk = 5;
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
T_hat = thval(A1, A2, A3);              % time consuming
T1 = t1val(T_hat, B4, dthdnu);
A1jm = A1jmval(q0, T1, M);              % time consuming
Ajm = Ajmval(q0, T1, M, Lu);            % time consuming


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
hline = repmat('-', 1, 70);
fprintf('\n');
fprintf('%s \n', hline);
fprintf(' Checking validity of approximate inverse\n');
fprintf('%s \n', hline);
fprintf(' * %8s %6.4f\n', 'rho = ', rho);
fprintf(' * %8s %6.4f\n', 'beta = ', beta);
fprintf('\n');
fprintf(' * Max. error btw a_{j,m} computed via\n');
fprintf('   convolution vs. direct fft:          ');
fprintf('%15.4e\n', norm(A-Ach, 'inf'))

%% inspecting (super-/sub-)diagonal terms
imagesc(A);
axis image, colorbar
fprintf('\n');
fprintf(' * Mean of k-th diagonal terms:\n');
fprintf('\n');
fprintf('    The 2nd super-diagonal terms are dominant.\n');
fprintf('    (Check MatLab figure window for colored \n');
fprintf('     image of the matrix A.)\n');
fprintf('\n');
fprintf(' %5s %13s \n', 'k', 'mean');
fprintf(' %19s\n', repmat('-', 1, 17));
for k = 2+(-3:3)
    fprintf(' %5d %13.4g\n', k, mean(diag(A, k)));
end

%% Lsinv, approximate inverse, without block structure
nb2 = length(i_rng);
I = eye(nb2);
D = diag( -i_rng );
Dinv = inv(D);

L = D + A;
Linv = inv(L);                          % not a recommended practice
sh = 1;
d(0+sh) = mean(diag(A, 0));
d(1+sh) = mean(diag(A, 1));
d(2+sh) = mean(diag(A, 2));
d(3+sh) = mean(diag(A, 3));

N0 = d(0+sh)*Dinv*diag(ones(nb2-0,1), 0); % diagonal
N1 = d(1+sh)*Dinv*diag(ones(nb2-1,1), 1); % 1-sup. diagonal; nilpotent
N2 = d(2+sh)*Dinv*diag(ones(nb2-2,1), 2); % 2-sup. diagonal; nilpotent
Lsinv = (I - N0 - N1 - N2)*Dinv;        % Neumann series, linear
fprintf('\n');
fprintf(' * Calculating the norm || I - Ls^{I}*L ||_2: \n');
fprintf('   ( dim(L) = %d x %d )\n', size(I));
fprintf('\n')
fprintf('   ** Without block structure: %8.4f\n', norm(I - Lsinv*L));
fprintf('      Here Ls^{I} = (I - As)*D^{-1} is a banded matrix\n');
fprintf('      of width 3 - a diagonal and 2 super-diagonals.\n');

%% Lsinv, approximate inverse, with block structure
% kk = 5;
ind_f = 1:kk;
ind_c = kk+1:nb2;
nf = length(ind_f);
nc = length(ind_c);

A11 = A(ind_f, ind_f);
A12 = A(ind_f, ind_c);
A21 = A(ind_c, ind_f);
A22 = A(ind_c, ind_c);

% Forming As with block structure: combining A12 and A22
for k = 0:4
    dd(k+sh) = mean(diag([A12; A22], -kk+k));
end
A11s = triu(A11);                       % upper triangular
A12s = zeros(size(A12));                % lower triangular
A12s(end, 1) = dd(1+sh);
A12s(end, 2) = dd(2+sh);
A12s(end-1, 1) = dd(2+sh);
A12s(end, 3) = dd(3+sh);
A12s(end-1,2) = dd(3+sh);
A12s(end-2, 1) = dd(3+sh);
A12s(end, 4) = dd(4+sh);
A12s(end-1, 3) = dd(4+sh);
A12s(end-2, 2) = dd(4+sh);
A12s(end-3, 1) = dd(4+sh);
A21s = zeros(size(A21));                % zero matrix
A22s = dd(0+sh)*diag(ones(nc-0,1), 0);  % banded matrix
for k = 1:4
    A22s = A22s + dd(k+sh)*diag(ones(nc-k,1), k);
end


As_blk = [A11s A12s; A21s A22s];
Lsinv_blk = (I - Dinv*As_blk)*Dinv;
fprintf('\n');
fprintf('   ** With block structure:    %8.4f\n', norm(I - Lsinv_blk*L));
fprintf('      Here, Ls^{I} has a 2-by-2, upper-triangular\n');
fprintf('      block structure, with K = %d.\n', kk);


% % simplified
% for k = 0:nc-1                          % mean diagonals of A22
%     dd(k+sh) = mean(diag(A22, k));      % get stablized as n gets large
% end
% A11s = triu(A11);                       % upper triangular
% A12s = tril(A12);                       % lower triangular
% A21s = zeros(size(A21));                % zero matrix
% A22s = dd(0+sh)*diag(ones(nc-0,1), 0);
% for k = 1:nc-1
%     A22s = A22s + dd(k+sh)*diag(ones(nc-k,1), k);
% end

%%
figure(1)
subplot(221), imagesc(A(ind_f, ind_f)), colorbar
subplot(222), imagesc(A(ind_f, ind_c)), colorbar
subplot(223), imagesc(A(ind_c, ind_f)), colorbar
subplot(224), imagesc(A(ind_c, ind_c)), colorbar

figure(2)
subplot(221), imagesc(A11s), colorbar
subplot(222), imagesc(A12s), colorbar
subplot(223), imagesc(A21s), colorbar
subplot(224), imagesc(A22s), colorbar

%% save to a binary file
fid = fopen('data.matlab/amat_b1_n08_r9100.bin', 'w+');
fwrite(fid, A, 'double');
fclose(fid);
fid = fopen('data.matlab/amat_b1_n08_r9100.bin');
Aread = fread(fid, [nb2,nb2], 'double');

%% 
fid = fopen('data.matlab/amat_b1_n08_r9100.dat', 'w+');
fprintf(fid, '%26.16e\n', A);
fclose(fid);
Araw = load('data.matlab/amat_b1_n08_r9100.dat');
Atxt = reshape(Araw, nb2, nb2);
