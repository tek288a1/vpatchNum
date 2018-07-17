%% script m-file: saveit.m
clear, clc

which_beta = 1;
init = load('vp_init_info.dat');
rho_init = init(1);
rho_term = 0.9;
a = init(2:end);

rho_vec = [rho_init:0.05:rho_term];
nr_step = 2;
nr_term = 128;
RHO = rho_vec(2:end);
BETA = zeros(size(RHO));
A = zeros(nr_term, length(RHO));
for j = 1:length(rho_vec)-1
    [a, ~, beta] = ...
        vsolver(rho_vec(j), rho_vec(j+1), nr_step, which_beta, a, nr_term);
    BETA(j) = beta;
    A(:,j) = a;
end

vp.RHO = RHO;
vp.BETA = BETA;
vp.A = A;
vp.nr_term = nr_term;
vp.nr_data = length(RHO);
% save('vp_r05to90_b1_n128.mat', '-struct', 'vp')

%%
clc
for rho = 0.1:0.05:0.9
    id = find(vp.RHO==rho);
    a = vp.A(:, id);
    dat_name = ...
        sprintf('%s%4.4d%s',...
        '../../VP_fortran/data/vp_b1_n07_r', round(rho*1e4), '.dat');
    dat_fort = load(dat_name);
    a_fort = dat_fort(4:end);
    fprintf('%6.4f \t %12.5e\n', rho, norm(a - a_fort))
end

%% as a data file; for compatibility with
fid = fopen('test.dat', 'wt');
fprintf(fid, '%13s %10.4f\n', 'rho = ', vp.rho);
fprintf(fid, '%13s %10.4f\n', 'beta = ', vp.beta);
fprintf(fid, '%13s %10d\n\n', 'n = ', vp.n);
fprintf(fid, '%24.16e \n', vp.a);
fclose(fid);

%% from matlab data file to text file; betaopt=0
clear
vp = load('vp_r05to90_b0_n128.mat');
n = 128;
nexp = log2(n);
betaopt = 0;
for id = 1:vp.nr_data
    rho = vp.RHO(id);
    a = vp.A(:, id);
    beta = vp.BETA(id);
    filename = sprintf('data/vp_b%1d_n%02d_r%4d.matlab.dat', ...
        betaopt, nexp, round(rho*1e4))
    fid = fopen(filename, 'wt');
    fprintf( fid, '%24d\n', n );
    fprintf( fid, '%24.4f\n', rho );
    fprintf( fid, '%24.16f\n', beta );
    fprintf( fid, '%24.16e\n', a );
    fclose(fid);
end

%% from matlab data file to text file; betaopt=1
clear
vp = load('vp_r05to90_b1_n128.mat');
n = 128;
nexp = log2(n);
betaopt = 1;
for id = 1:vp.nr_data
    rho = vp.RHO(id);
    a = vp.A(:, id);
    beta = vp.BETA(id);
    filename = sprintf('data/vp_b%1d_n%02d_r%4d.matlab.dat', ...
        betaopt, nexp, round(rho*1e4))
    fid = fopen(filename, 'wt');
    fprintf( fid, '%24d\n', n );
    fprintf( fid, '%24.4f\n', rho );
    fprintf( fid, '%24.16f\n', beta );
    fprintf( fid, '%24.16e\n', a );
    fclose(fid);
end
