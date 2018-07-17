function [a0, q0] = a0valnew(z, zd, u)
%A0VALNEW [a0, q0] = a0valnew(z, zd, u)
%
%  Calculation of A0(nu) and q0(nu) by alternating point
%  trapezoidal method and spectral differentiation. -u-a0 is the complex
%  velocity.
%
%  INPUTS:
%    z = f0(nu'), length(z) needs to be an even number
%    zd = f0'(nu')
%
%   Written by Tae Eun Kim, 02/23/2018
    n = length(z);
    nb2 = n/2;
    a0 = zeros(size(z));
    zo = z(1:2:n);
    ze = z(2:2:n);
    zdo = zd(1:2:n);
    zde = zd(2:2:n);
    for i = 1:n
        if mod(i, 2) == 1               % odd
            zc = zo(ceil(i/2));
            zq = ze;
            zdq = zde;
        else                            % even
            zc = ze(i/2);
            zq = zo;
            zdq = zdo;
        end
        a0(i) = 1/nb2* sum( ...
            (abs(zc)^2 - abs(zq).^2)./(zc^2 - zq.^2).*zdq ...
            );
    end
    q0 = 1./(u + a0);
    % % spectral differentiation of a0
    % a0ft = real( fft(a0)/n );               % real coeffs
    % dfac = 1i*fftshift([0, -n/2+1:n/2-1]'); % column vector
    % a0d = ifft( dfac.*a0ft )*n;
end
