function [a0, q0] = a0val(zz, z, zd, u)
%A0VAL [a0, q0] = a0val(zz, z, zd, u)
%
%  Calculation of A0(nu) and q0(nu) by alternating point
%  trapezoidal method and spectral differentiation. -u-a0 is the complex
%  velocity.
%
%  INPUTS:
%    zz = f0(nu) at collocation points
%    z = f0(nu') at quadrature points
%    zd = f0'(nu') at quadrature points;
%
%   Written by Tae Eun Kim, 12/03/2017
    n = length(zz);
    nq = length(z);
    a0 = zeros(size(zz));
    for i = 1:n
        a0(i) = 1/nq*...
                sum( (abs(zz(i))^2 - abs(z).^2)./(zz(i)^2 - z.^2).*zd );
    end
    q0 = 1./(u + a0);
    % % spectral differentiation of a0
    % a0ft = real( fft(a0)/n );               % real coeffs
    % dfac = 1i*fftshift([0, -n/2+1:n/2-1]'); % column vector
    % a0d = ifft( dfac.*a0ft )*n;
end
