%% evenly spaced-out rho data
clear;
addpath('../functions')
n_in = 256;
K = 512;
beta_opt = 1;
indir = '~/Google Drive/Vpatch1/matlab1/data.fortran/';
outdir = '~/Google Drive/Vpatch1/matlab1/data.matlab/';

tic
%for rho_in = 0.90:0.01:0.95
rho_in = 0.9;
    infile = sprintf('vp_b%1d_n%02d_r%04d.dat',...
                     beta_opt, log2(n_in), round(rho_in*1e4));
    infilename = strcat(indir, infile);
    outfile = sprintf('amat_b%1d_n%02d_r%04d.bin',...
                      beta_opt, log2(K), round(rho_in*1e4));
    outfilename = strcat(outdir, outfile);

    raw = load( infilename );
    %n = raw(1);
    rho = raw(2);
    beta = raw(3);
    a = raw(4:end);
    
    fprintf('** Calculating A for rho = %7.4f **\n', rho);
    A = calcA(K, rho, beta, a);
    
    fid = fopen(outfilename, 'w+');
    fwrite(fid, A, 'double');
    fclose(fid);
%end
toc
rmpath('../functions')

%% working with a binary file
% writing A into a binary file
%   fid = fopen(filename, 'w+');
%   fwrite(fid, A, 'double');
%   fclose(fid);
%   
% reading a binary file
%   fid = fopen(filename);
%   A = fread(fid, [m, n], 'double');

