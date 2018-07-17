function T1 = t1val(T_hat, B4, dthdnu)
    %% calculation of t^{(1)}_{j,m} for -N <= j,m, < N
    %
    % t^{(1)}_{j,m} is defined as
    %   = \hat{t}_{j,m} - ( \hat{t}_{0,m} - m*B^{(4)}_{-m,0} ) v_j 
    %   
    % 1. Grab v_j for -2N <= j < 2N. (FFT on dthdnu)
    % 2. FFTSHIFT B4 (2Nx2N) and pad with zeros; call it B4c (4Nx4N).
    % 3. Take \hat{t}_{j,m} for -N <= j,m < N.
    % 4. Form t^{(1)}_{j,m}

    n4 = length(dthdnu);                % better be a multiple of 4
    n2 = n4/2;
    n = n2/2;
    
    % 1. We have point values of dthdnu at 4N equi-spaced points:
    v = fft(dthdnu)/length(dthdnu);
    v = fftshift(v);

    % 2. Do fftshift on B4 and put in the center of 4Nx4N matrix.
    B4c = zeros(n4, n4);
    i_ctr = (n+1:3*n);                  % indices for center
    B4c(i_ctr,i_ctr) = fftshift(B4);

    % 3. T = \hat{t}_{j,m} was already calculated using the
    % function tval.
    %    * dimensions of T = 2N x 2N
    %    * indices represented: -N <= j, m < N

    % 4. Let's calculate t^{(1)}_{j,m} for -N <= j,m < N
    s1 = n+1;                           % index shift for 2N x 2N
    s2 = n2+1;                          % index shift for 4N x 4N
    T1 = zeros(size(T_hat));
    for j = -n:n-1
        for m = -n:n-1
            T1(j+s1,m+s1) = T_hat(j+s1,m+s1) - ...
                (T_hat(0+s1,m+s1) - m*B4c(-m+s2,0+s2))*v(j+s2);
        end
    end
    % % The block above can be vectorized.
    % % Though the codes below are shorter, may not be clearer to
    % % readers. 
    % T0m = repmat(T_hat(0+s1,:), n2, 1);
    % mB4 = repmat( (-n:n-1).*B4c(n+s2:-1:-n+1+s2,0+s2).', n2, 1);
    % vj = repmat( v(-n+s2:n-1+s2), 1, n2);
    % T1 = T - (T0m - mB4) .* vj;
end
