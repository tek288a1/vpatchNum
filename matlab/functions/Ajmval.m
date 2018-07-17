function Ajm = Ajmval(q0, T1, M, Lu)
%% Calculation of a_{j,m} for 1 <= j,m <= N-1
% This will be done in two stages 
    
    n4 = length(q0);
    n2 = n4/2;
    n = n2/2;

    %% calculation of a^{(1)}_{j,m} for 0 <= j <= N-1 and 1 <= m <= N
    %
    % a^{(1)}_{j,m} is calculated by 
    %   = \sum_{k=-N}^{N} 
    %     Q(k)*( T1(j-k,m) - \sum_{l=1}^{N} T1(j-k,1-l)*M(l,m) )
    %
    % 1. Calculate Fourier coefficients Q(k) of q0 in fftshifted form.
    % 2. t^{(1)}_{j,m} for -N <= j,m < N have been calculated. We need
    %    values for N <= j < 2N. So pad T1 (2Nx2N) with zeros to form
    %    T1c (4Nx4N). 
    % 3. Form A1jm according to the formula.

    % 1. fft on q0 followed by fftshift
    Q = fftshift( fft(q0)/length(q0) );

    % 2. Recall that the zeroth modes of T1 are already in the middle of
    % spectrum, so no need for fftshift. We call the enlarged matrix
    % T1c.
    T1c = zeros(n4, n4);
    i_ctr = (n+1:3*n);
    T1c(i_ctr,i_ctr) = T1;

    % 3. Forming A1jm. Though A1jm is a NxN square matrix, keep in mind
    % that the indices represented are not symmetric, i.e.,
    %	0 <= j <= N-1  while  1 <= m <= N
    A1jm = zeros(n, n);
    sh = 1;                            % j-index of A1jm needs a shift
    s2 = n2+1;                         % T1c and Q requires this shift
    kk = (-n:n-1);                      % indices for k-summation
    ll = (1:n);                         % indices for l-summation

    %%
    for j = 0:n-1
        for m = 1:n
            % Forming summand. Note how l-summation is handled via
            % matrix-vector multiplication (inner product).
            summand = Q(kk+s2) .* ...
                  ( T1c(j-kk+s2,m+s2) - T1c(j-kk+s2,1-ll+s2)*M(ll,m) );
            A1jm(j+sh,m) = sum( summand );
        end
    end


    %% calculation of a_{j,m} for 1 <= j,m <= N-1
    %
    % a_{j,m} is calculated by 
    %   = a^{(1)}_{j,m} + a^{(1)}_{0,m}*\tilde{l}_j
    %
    % where \tilde{l}_j = - l_j / l_0 and l_j's are Fourier
    % coefficients of Lu defined by
    % 
    % Lu(nu) = q0(nu) * ( omega0(nu) + j(nu)*dthdnu(nu) )
    % 
    % Note that A1jm = ( a^{(1)}_{j,m} ) is already calculated for 
    % 	0 <= j <= N-1  and  1 <= k <= N .
    % 
    % 1. Calculate l_j via FFT on point-values of Lu. 
    % 2. Form Ajm which is (N-1)x(N-1).

    % 1. FFT on Lu. Note that length(Lu) = n4. Since we only need
    % positive modes of Lu, fftshift is not necessary. 
    l = fft(Lu)/length(Lu);
    lt = -l/l(1);

    % 2. Forming Ajm.
    mm = (1:n-1);                       % range of m-index
    jj = (2:n);                         % range of j-index
    A10m = repmat(A1jm(0+sh, mm), n-1, 1);
    ltj = repmat(lt(jj), 1, n-1);
    Ajm = A1jm(jj,mm) + A10m .* ltj;
    % %% doing with for-loops
    % Ajmch = zeros(n-1,n-1);
    % for j = 1:n-1
    %     for m = 1:n-1
    %         Ajmch(j,m) = A1jm(j+sh,m) + A1jm(0+sh,m)*lt(j+sh);
    %     end
    % end
end
