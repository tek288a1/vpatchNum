function calcA_fnc(n_in, rho_in, beta_opt, K)
% calcA_fnc is a function file calculating A matrix.
% inputs:  n_in      number of coefficients + velocity in the initial data
%          rho_in    rho
%          beta_opt  0, 1, 2
%          K         dimension of A
    addpath('../functions')
    indir = '../../data/data.fortran/';
    outdir = '../../data/data.matlab/';
    infile = ...
        sprintf('vp_b%1d_n%02d_r%04d.dat', ...
                beta_opt, log2(n_in), round(rho_in*1e4));
    outfile = ...
        sprintf('amat_b%1d_n%02d_r%04d.bin',...
                beta_opt, log2(K), round(rho_in*1e4));
    infilename = strcat(indir, infile);
    outfilename = strcat(outdir, outfile);
    raw = load( infilename );
%     n = raw(1);
    rho = raw(2);
    beta = raw(3);
    a = raw(4:end);
    
    fprintf('** Calculating A for rho = %7.4f and K = %5d **\n', rho, K);
    tic
%     A = calcA(K, rho, beta, a);
    A = rand(K);
    toc
    
    fid = fopen(outfilename, 'w+');
    fwrite(fid, A, 'double');
    fclose(fid);
    rmpath('../functions')
end
