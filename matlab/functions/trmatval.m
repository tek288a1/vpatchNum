function M = trmatval(n, rho, beta)
%% transfer matrix calculation
% input:
%   n = dimension of matrix
%   rho
%   beta
% output:
%   M = transfer matrix
%
% This m-file calculates the "transfer matrix" M which relates the
% positive Fourier coefficients of functions in S^2 to the negative
% modes. The matrix entry is calculated by the formula
%
%   M(k,j) = 1/(2\pi) \int_{0}^{2\pi} \hat{\eta}^j \eta^{k-1} d\nu
%
% for 1 \le j,k \le N. Since eta = exp(i\nu), we note that M(k,j) =
% F[\hat{\eta}^j](1-k). So we can use FFT to calculate the entries
% efficiently.
%
% 03/01/2018
  M = zeros(n, n);
  nu = linspace(0, 2*pi, 2*n+1)'; % need at least twice as many points
  nu(end) = [];
  eta = exp(1i*nu);
  zeta = (eta-beta)./(1-beta*eta);
  etah = (rho^2./zeta+beta)./(1+beta*rho^2./zeta);
  P = etah.^(1:n);
  P = fft(P, [], 1)/(2*n);
  M = real([P(1, :); P(end:-1:end-n+2, :)]);
end
