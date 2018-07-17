function T = tval(A1, A2, A3, B4, dthdnu)
%% calculation of t(j,m) = \tilde{t}_{j,m} + m*B^{(4)}_{-m,0} v_j
% 1. 2-D fftshift on 2N x 2N matrices
% 2. pad with zeros to have enlarged 4N x 4N matrices
% 3. calculate T that is 2N x 2N
% 
%%%% THIS IS OUTDATED %%%%

    n2 = size(A1, 1);
    n4 = 2*n2;
    n = n2/2;
    
    A1c = zeros(n4, n4);
    A2c = zeros(n4, n4);
    A3c = zeros(n4, n4);
    B4c = zeros(n4, n4);
    T = zeros(n2, n2);
    sh = n2+1;                              % index shift for 4N x 4N
    sh1 = n+1;                              % index shift for 2N x 2N
    i_ctr = (n+1:3*n);                      % indices for center

    A1c(i_ctr, i_ctr) = fftshift(A1);
    A2c(i_ctr, i_ctr) = fftshift(A2);
    A3c(i_ctr, i_ctr) = fftshift(A3);
    B4c(i_ctr, i_ctr) = fftshift(B4);

    for j = -n:n-1
        for m = -n:n-1
            % tmp vector below corresponds to the slices of matrices
            % A1 and A2 with indices given by
            % 	A1(j-k1,k1-m) and A2(j-k2,k2+m)
            % with k1 between 0 and m and
            % with k2 between 0 and -m for given j and m.            
            
            % tmp = diag(...
            %     A1c(j+sh:-sign(m):j-m+sh, 0-m+sh:sign(m):m-m+sh) + ...
            %     A2c(j+sh:sign(m):j+m+sh, 0+m+sh:-sign(m):-m+m+sh) );
            tmp = zeros(abs(m)+1,1);
            kk1 = (min(0,m):max(0,m));
            kk2 = (min(0,-m):max(0,-m));
            tmp = diag(...
                A1c(j-kk1+sh,kk1-m+sh) + A2c(j-kk2+sh,kk2+m+sh) );            
            
            T(j+sh1,m+sh1) = ...
                - A1c(j-m+sh,0+sh) + A1c(j+sh,-m+sh) ...
                - A2c(j+m+sh,0+sh) + A2c(j+sh,m+sh) ...
                - A3c(j+sh,-m+sh) ...
                - sign(m)*sum( tmp(1:end-1) + tmp(2:end) );
            % T(j+sh1,m+sh1) = ...
            %     - A1c(j-m+sh,0+sh) + A1c(j+sh,-m+sh) ...
            %     - A2c(j+m+sh,0+sh) + A2c(j+sh,m+sh) ...
            %     - A3c(j+sh,-m+sh);            
        end
    end
end

% $$$ %% calculation of t0(j,m) = \tilde{t}_{j,m}
% $$$ % 1. 2-D fftshift on 2N x 2N matrices
% $$$ % 2. pad with zeros to have enlarged 4N x 4N matrices
% $$$ % 3. calculate t0 that is 2N x 2N
% $$$ 
% $$$ A1c = zeros(n4, n4);
% $$$ A2c = zeros(n4, n4);
% $$$ A3c = zeros(n4, n4);
% $$$ B4c = zeros(n4, n4);
% $$$ T01 = zeros(n2, n2);
% $$$ sh = n2+1;                              % index shift for 4N x 4N
% $$$ sh1 = n+1;                              % index shift for 2N x 2N
% $$$ i_ctr = (n+1:3*n);                      % indices for center
% $$$ A1c(i_ctr, i_ctr) = fftshift(A1);
% $$$ A2c(i_ctr, i_ctr) = fftshift(A2);
% $$$ A3c(i_ctr, i_ctr) = fftshift(A3);
% $$$ B4c(i_ctr, i_ctr) = fftshift(B4);
% $$$ for j = -n:n-1
% $$$     for m = -n:n-1
% $$$         % tmp vector below corresponds to the slices of matrices
% $$$         % A1 and A2 with indices given by
% $$$         % 	A1(j-k1,k1-m) and A2(j-k2,k2+m)
% $$$         % with k1 between 0 and m and
% $$$         % with k2 between 0 and -m for given j and m.
% $$$         tmp = zeros(abs(m)+1,1);
% $$$         kk1 = (min(0,m):max(0,m));
% $$$         kk2 = (min(0,-m):max(0,-m));
% $$$         tmp = diag( A1c(j-kk1+sh,kk1-m+sh) + A2c(j-kk2+sh,kk2+m+sh) );
% $$$         T01(j+sh1,m+sh1) = ...
% $$$             - A1c(j-m+sh,0+sh) + A1c(j+sh,-m+sh) ...
% $$$             - A2c(j+m+sh,0+sh) + A2c(j+sh,m+sh) ...
% $$$             - A3c(j+sh,-m+sh) ...
% $$$             - sign(m)*sum( tmp(1:end-1) + tmp(2:end) );
% $$$     end
% $$$ end
