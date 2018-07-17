function A1g = A1gval(m, rho, beta,  z, zd)
% A1GVAL - A1g = A1gval(m, z, zd) 
%   evaluates A_1[g_m](nu) at N equispaced points between 0 and 2*pi
%   (including 0, excluding 2*pi) where A_1[g] is defined as
%
%   1/(4*pi) \int_{0}^{2*pi} d th ...
%      zdth ( conj(gnu+gth)/(znu+zth) + conj(gnu-gth)/(znu-zth) )
%      - conj(zdth) ( (gnu+gth)/(znu+zth) + (gnu-gth)/(znu-zth) )
%   
%   and g_m(nu) = i*exp(i*m*nu)
%
    n2 = length(z);                     % better be an even number
    n = n2/2;
    A1g = zeros(size(z));
    io = (1:2:n2);                      % odd indices
    ie = (2:2:n2);                      % even indices
    nu = linspace(0, 2*pi, n2+1)';
    nu(end) = [];
    eta = exp(1i*nu);
    zeta = (eta - beta)./(1 - beta*eta);
    etah = (rho^2./zeta + beta)./(1 + beta*rho^2./zeta);
    g = 1i*(eta.^m - etah.^m);
%     g = 1i*exp(1i*m*nu);
    
    gnu = repmat( g([io, ie]), 1, n );
    znu = repmat( z([io, ie]), 1, n );
    gth = [repmat( g(ie).', n, 1 ); repmat( g(io).', n, 1 )];
    zth = [repmat( z(ie).', n, 1 ); repmat( z(io).', n, 1 )];
    zdth = [repmat( zd(ie).', n, 1 ); repmat( zd(io).', n, 1 )];
    
    integ = zdth.*( conj( gnu+gth )./( znu+zth ) + ...
                    conj( gnu-gth )./( znu-zth ) ) - ...
            conj(zdth).*( ( gnu+gth )./( znu+zth ) + ...
                          ( gnu-gth )./( znu-zth ) );
    A1g = 1/(2*n)*sum( integ, 2 );
    A1g([io, ie]) = A1g;
end
