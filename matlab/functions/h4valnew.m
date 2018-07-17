function [H, B] = h4valnew(z, zd)
%% [H, B] = h4valnew(z, zd)
  z = z(:);
  zd = zd(:);
  n = length(z);
  Znu = repmat( z, 1, n );
  Zth = repmat( z.', n, 1 );
  Zdth = repmat( zd.', n, 1 );
  H = Zdth .* ( abs(Znu).^2 - abs(Zth).^2 )./( Znu.^2 - Zth.^2 );
  id_diag = find( eye(n)==1 );
  H(id_diag) = zd/2.*...
      ( conj(zd)./zd + conj(z)./z );
  B = fft2(H)/n^2;
end
