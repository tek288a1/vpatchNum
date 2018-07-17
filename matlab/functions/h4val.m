function [H, B] = h4val(z, zd)
  npt = length(z);  % npt needs to be an even number
  nt = npt/2;       % number of transform modes
  io = 1:2:npt;
  ie = 2:2:npt;
  znu = z(io);
  zth = z(ie);
  zdth = zd(ie);
  Znu = repmat(znu, 1, nt);
  Zth = repmat(zth.', nt, 1);
  Zdth = repmat(zdth.', nt, 1);
  % point values
  H = Zdth .* ( abs(Znu).^2 - abs(Zth).^2 ) ./ ( Znu.^2 - Zth.^2 );

  % double Fourier coefficients 
  B = fft2(H)/nt^2;                 % one 2-D FFT
  v = fftshift(-nt/2+(0:nt-1))/nt;  % phase shift along rows (theta)
  psh = exp(-1i*pi*v);
  B = B .* psh;
  B(nt/2+1, :) = 0;
  B(:, nt/2+1) = 0;
end 
