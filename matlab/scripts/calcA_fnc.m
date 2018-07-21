function calcA_fnc(n_in, rho_in, beta_opt, K)
%% evenly spaced-out rho data
addpath('../functions')
%n_in = 256;
%K = 32;
%beta_opt = 1;
%rho_in = 0.85;
indir = '../../data/data.fortran/';
outdir = '../../data/data.matlab/';

tic
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
    
    fprintf('** Calculating A for rho = %7.4f and K = %5d **\n', rho, K);
    A = calcA(K, rho, beta, a);
    
    fid = fopen(outfilename, 'w+');
    fwrite(fid, A, 'double');
    fclose(fid);
toc
rmpath('../functions')
end
