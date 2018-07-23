%% Matrix A calculation
% 07/22/18

clear
addpath('../functions');
fname = '../../data/data.matlab/amat_b1_n10_r8500.bin';
fid = fopen(fname);
n = 1024;
A = fread(fid, [n/2, n/2], 'double');
fclose(fid);

norm(A)
size(A)

%%
outfilename = 'atest.bin';
fid = fopen(outfilename, 'w+');
fwrite(fid, A, 'double');
fclose(fid);

%% testing
n_in = 256;
rho_in = 0.9;
beta_opt = 1;
K = 1024;
indir = '../../data/data.fortran/';
infile = ...
    sprintf('vp_b%1d_n%02d_r%04d.dat', ...
    beta_opt, log2(n_in), round(rho_in*1e4));
outfile = ...
    sprintf('amat_b%1d_n%02d_r%04d.bin',...
    beta_opt, log2(K), round(rho_in*1e4));
infilename = strcat(indir, infile);
outfilename = strcat(outfile);
raw = load( infilename );
%     n = raw(1);
rho = raw(2);
beta = raw(3);
a = raw(4:end);

%%
fprintf('** Calculating A for rho = %7.4f and K = %5d **\n', rho, K);
% tic
% A = calcA(K, rho, beta, a);
% toc

fid = fopen(outfilename, 'w+');
fwrite(fid, A, 'double');
fclose(fid);
rmpath('../functions')

