function T_hat = thval(A1, A2, A3)
    %% calculation of \hat{t}_{j,m} for -N <= j,m < N
    %
    % \hat{t}_{j,m} is defined as
    %   = - A1(j-m,0) + A1(j,-m) - A2(j+m,0) + A2(j,m) - A3(j,-m)
    %     - 2*sgn(m)* \sum_{k=0}^{m}' (A1(j-k,k-m) + A2(j-k,k+m))
    %     
    % where \sum_{k=0}^{m}' a(k)
    %   = a(0)/2 + \sum_{k=1}^{m-1} a(k) + a(m)/2 .
    %   
    % 1. FFTSHIFT 2Nx2N matrices (A1, A2, A3) and pad with zeros to
    %    form 4Nx4N matrices. 
    % 2. Form T which is a 2Nx2N matrix.
        
    n2 = size(A1, 1);                   % better be a multiple of 2
    n4 = 2*n2;
    n = n2/2;

    % 1. 4Nx4N matrices with zeroth-modes in the middle of spectrum
    A1c = zeros(n4, n4);
    A2c = zeros(n4, n4);
    A3c = zeros(n4, n4);
    i_ctr = (n+1:3*n);                  % indices for center
    A1c(i_ctr, i_ctr) = fftshift(A1);
    A2c(i_ctr, i_ctr) = fftshift(A2);
    A3c(i_ctr, i_ctr) = fftshift(A3);

    % 2. In order to enhance readability, introduce index
    % shifts. The codes will resemble mathematical expressions. 
    s1 = n+1;                           % index shift for 2N x 2N
    s2 = n2+1;                          % index shift for 4N x 4N
    T_hat = zeros(n2, n2);
    for j = -n:n-1
        for m = -n:n-1
            % tmp vector below corresponds to the slices of matrices
            % A1 and A2 with indices given by
            % 	A1(j-k1,k1-m) and A2(j-k2,k2+m)
            % with k1 between 0 and m and
            % with k2 between 0 and -m for given j and m.            
            kk1 = (min(0,m):max(0,m));
            kk2 = (min(0,-m):max(0,-m));
            tmp = diag(...
                A1c(j-kk1+s2,kk1-m+s2) + A2c(j-kk2+s2,kk2+m+s2) );            
            T_hat(j+s1,m+s1) = ...
                - A1c(j-m+s2,0+s2) + A1c(j+s2,-m+s2) ...
                - A2c(j+m+s2,0+s2) + A2c(j+s2,m+s2) ...
                - A3c(j+s2,-m+s2) ...
                - sign(m)*sum( tmp(1:end-1) + tmp(2:end) );
        end
    end
end
