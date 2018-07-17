%% loading data
clear;
% Add path (at beginning of script)
added_path = ['../functions']; 
addpath(added_path);
% n_in = input(' Input N = ');
rho_in = input(' Input rho = ');
% beta_opt = input(' Input beta_opt = ');
% kk = input(' Input K = ');
n_in = 256;
% rho_in = 0.93;
beta_opt = 1;
kk = 5;
ddir = '../../data/data.fortran/';
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
% ind_trunc = 1:100;
n_coef = length(a)-1;
ind_trunc = 1:length(a);
[z, zd, zeta, dthdnu] = ptval(rho, beta, a(ind_trunc), eta);
M = trmatval(n, rho, beta);
jt = jtval(nj, rho, zeta);
[A0, q0] = a0valnew(z, zd, u);
c0 = c0val(A0, zd);
omega0 = q0 .* dthdnu .* (c0 - u*jt);
Lu = q0 .* (omega0 + jt.*dthdnu);
% [J1, A1] = j1val(eta, z, zd, omega0);
% [J2, A2] = j2val(eta, z, zd, omega0);
% [J3, A3] = j3val(eta, z, zd, omega0);
% [H4, B4] = h4val(z, zd);
% T_hat = thval(A1, A2, A3);              % time consuming
% T1 = t1val(T_hat, B4, dthdnu);
% A1jm = A1jmval(q0, T1, M);              % time consuming
% Ajm = Ajmval(q0, T1, M, Lu);            % time consuming

% fft's
f0ft = fft(z)/n4;
q0ft = fft(q0)/n4;
q0tnft = fft(q0.*dthdnu)/n4;            % q0*dthdnu
r0tn = zd - omega0;                     % r0*dthdnu
r0tnft = fft(r0tn)/n4; 

k_vec = fftshift( [0:n4-1]-n2 )';       % multiplication by k
k_vec1 = k_vec;                         % preserve 0th mode
k_vec1(1) = 1;

%% f0p, f0m, and q1m
f0p_inf = norm( k_vec .* f0ft, 1 );
f0p_h0 = norm( k_vec1 .* f0ft, 2 );
f0p_h1 = norm( k_vec1.^2 .* f0ft, 2 );
f0m = abs(z(1));
q1m = min(abs(zd));
% plot( nu, real(zd), nu, imag(zd), nu, abs(zd) )
% subplot(121); plot( nu, abs(zd), nu, q1m*ones(size(nu)) )
% subplot(122); plot( real(zd(1:n2)), imag(zd(1:n2)) )

%% check: confirmed
% Note that we are only taking values from 0 to pi since we want
% |nu-theta|<pi,
Nu = repmat(nu(1:2:n2), 1, n);        % nu, column
Th = repmat(nu(2:2:n2).', n, 1);      % theta, row
Znu = repmat(z(1:2:n2), 1, n);        % nu, column
Zth = repmat(z(2:2:n2).', n, 1);      % theta, row
Q1 = (Znu - Zth)./(Nu - Th);
min(min(abs(Q1)));
clf;
imagesc(abs(Q1)); colorbar, axis equal


%% q0
q0_inf = norm(q0ft, 1);
q0_h1 = norm( k_vec1 .* q0ft, 2 );

%% q0*dthdnu
q0tn_h1 = norm( k_vec1 .* q0tnft, 2 );


%% r0*dthdnu
% subplot(121); plot(nu, abs(r0tn))
% subplot(122); semilogy(1:n4, abs(r0tnft))
r0tn_inf = norm(r0tnft, 1);
r0tn_h1 = norm( k_vec1 .* r0tnft, 2 );

%% transfer matrix M and K_inf, K_1, and K_2
% recall that M is n x n matrix
K_inf = norm( (1 + sum(M, 1)) ./ ((1:n).^2), 2);
l_diag = diag(0:n-1); 
l_diag(1) = 1;
k_diag = diag((1:n).^(-2));
K_1 = norm( l_diag * M * k_diag, 'fro' );
K_2 = norm( l_diag^2 * M * k_diag, 'fro' );

%% beta's
beta1 = K_inf/q1m;
beta2 = 2*K_2/(3*q1m) + K_inf*f0p_h1/(3*q1m^2);
beta3 = K_inf/f0m;
beta4 = K_1/(2*f0m) + K_inf*f0p_inf/(2*f0m^2);
beta5 = (beta1 + beta3)*f0p_h0;
beta6 = (beta2 + beta4)*f0p_h0;

%% write to a file
result_file = sprintf('result_r%4d.dat', rho_in*1e4);
fid = fopen(result_file, 'w+');
fprintf(fid, 'rho       %8.2f\n', rho);
fprintf(fid, 'beta      %8.2f\n', beta);
fprintf(fid, 'n_coeff   %8d  \n', n_coef);
fprintf(fid, 'n_fft     %8d  \n', n4);
fprintf(fid, 'f0p_inf   %8.2f\n', f0p_inf);
fprintf(fid, 'f0p_h0    %8.2f\n', f0p_h0);
fprintf(fid, 'f0p_h1    %8.2f\n', f0p_h1);
fprintf(fid, 'f0m       %8.2f\n', f0m);
fprintf(fid, 'q1m       %8.2f\n', q1m);
fprintf(fid, 'q0_inf    %8.2f\n', q0_inf);
fprintf(fid, 'q0_h1     %8.2f\n', q0_h1);
fprintf(fid, 'q0tn_h1   %8.2f\n', q0tn_h1);
fprintf(fid, 'r0tn_inf  %8.2e\n', r0tn_inf);
fprintf(fid, 'r0tn_h1   %8.2e\n', r0tn_h1);
fprintf(fid, 'K_inf     %8.2f\n', K_inf);
fprintf(fid, 'K_1       %8.2f\n', K_1);
fprintf(fid, 'K_2       %8.2f\n', K_2);
fprintf(fid, 'beta1     %8.2f\n', beta1);
fprintf(fid, 'beta2     %8.2f\n', beta2);
fprintf(fid, 'beta3     %8.2f\n', beta3);
fprintf(fid, 'beta4     %8.2f\n', beta4);
fprintf(fid, 'beta5     %8.2f\n', beta5);
fprintf(fid, 'beta6     %8.2f\n', beta6);
fclose(fid);







%% the following block checks.
% n_coef = length(a)-1;
% c_pos = 1i*a(2:end);
% c_neg = -M(1:end-1,1:end-1)*c_pos;
% cm_0 = 1i*(1+beta*rho)/(1-beta*rho);
% cm_tmp = 1i*2*rho*(1-beta^2)/(rho-beta)/(1-beta*rho)*((beta-rho)/(1-beta*rho)).^[1:n_coef-1]';
% cm_neg = [cm_0; cm_tmp];
% ft_pos = ft(2:n_coef+1);
% ft_neg = [ft(1); ft(end:-1:end-n_coef+2)];
% norm(ft_pos-c_pos, inf)
% norm(ft_neg-(c_neg+cm_neg), inf)
% 
% % norms of f0'
% zd_h1 = norm([ft_pos.*([1:n_coef]').^2; ft_neg.*([0:n_coef-1].').^2],2)
% zd_inf = norm([ft_pos.*([1:n_coef]'); ft_neg.*([0:n_coef-1].')], 1)
% % sqrt(sum( (abs(ft_pos).^2) .* ([1:n_coef]').^4 ) + ...
% %     sum( (abs(ft_neg).^2) .* ([0:n_coef-1].').^4 )) % check



%% Remove path (at end of script/script clean-up)
rmpath(added_path);
