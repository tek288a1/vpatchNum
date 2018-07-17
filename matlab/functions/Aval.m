function A = Aval(q0, T1, M)
    %% calculation of a_{j,m} for 0 <= j <= N-1 and 1 <= m <= N
    %
    % a_{j,m} is calculated by 
    %   = \sum_{k=-N}^{N} 
    %     Q(k)*( T1(j-k,m) - \sum_{l=1}^{N} T1(j-k,1-l)*M(l,m) )
    %
    % 1. Calculate Fourier coefficients Q(k) of q0 in fftshifted form.
    % 2. t^{(1)}_{j,m} for -N <= j,m < N have been calculated. We need
    %    values for N <= j < 2N. So pad T1 (2Nx2N) with zeros to form
    %    T1c (4Nx4N). 
    % 3. Form A according to the formula.

    n4 = length(q0);                    % better be multiple of 4
    n2 = n4/2;                          % T1 is 2Nx2N
    n = n2/2;                           % M and A are NxN
    
    % 1. fft on q0 followed by fftshift
    Q = fftshift( fft(q0)/length(q0) );

    % 2. Recall that the zeroth modes of T1 are already in the middle of
    %    spectrum, so no need for fftshift. We call the enlarged matrix
    %    T1c.
    T1c = zeros(n4, n4);
    i_ctr = (n+1:3*n);
    T1c(i_ctr,i_ctr) = T1;

    % 3. Forming A. Though A is a NxN square matrix, keep in mind that
    % the indices represented are not symmetric, i.e., 
    %	0 <= j <= N-1  while  1 <= m <= N
    A = zeros(n, n);
    sh = 1;                                 % j-index of A needs a shift
    s2 = n2+1;                              % T1c and Q requires this shift
    kk = (-n:n-1);                          % indices for k-summation
    ll = (1:n);                             % indices for l-summation

    %%
    for j = 0:n-1
        for m = 1:n
            % Forming summand. Note how l-summation is handled via
            % matrix-vector multiplication (inner product).
            smd = Q(kk+s2) .* ...
                  ( T1c(j-kk+s2,m+s2) - T1c(j-kk+s2,1-ll+s2)*M(ll,m) );
            A(j+sh,m) = sum( smd );
        end
    end
end
