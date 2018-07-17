function [J3, A3] = j3val(eta, z, zd, omega)
    npt = length(z);	% number of points need to be even
    nt = npt/2;		% number of transform modes
    io = 1:2:npt;
    ie = 2:2:npt;
    Znu = repmat( z(io), 1, nt );
    Zth = repmat( z(ie).', nt, 1 );
    Zdth = repmat( zd(ie).', nt, 1);
    Omnu = repmat( omega(io), 1, nt );
    Einu = repmat( eta(io), 1, nt );
    Eith = repmat( eta(ie).', nt, 1 ); 
    
    J3 = 1i*Omnu .* conj(Zdth) ./ (Znu + Zth);
    A3 = fft2(J3)/nt^2;
    v = fftshift((-nt/2:nt/2-1))/nt;
    psh = exp(-1i*pi*v);
    A3 = A3 .* psh;
    A3(nt/2+1, :) = 0;
    A3(:, nt/2+1) = 0;
end 
