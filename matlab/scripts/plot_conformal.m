FS = 'FontSize';
fs = 20;

rho = 0.8;
beta = (1 - sqrt(1-rho^2))/rho;
nr = 3;
nth = 361;
r = linspace(rho, 1, nr)';
theta = linspace(0, 2*pi, nth);

%%
clf
eta = r.*exp(1i*theta);
zeta = (eta - beta)./(1 - beta*eta);
etah = (rho^2./zeta + beta)./(1 + beta*rho^2./zeta);

figure(1);
subplot(1, 2, 1)
hold on
for i = 1:nr
    plot(real(eta(i,:)), imag(eta(i,:)), 'b');
end
for j = 1:15:nth
    plot(real(eta(:,j)), imag(eta(:,j)), 'b');
end
plot(beta, 0, 'ro')
title('\eta-plane', FS, fs)
axis image
axis([-1.01 1.01 -1.01 1.01])

subplot(1, 2, 2)
hold on
for i = 1:nr
    plot(real(etah(i,:)), imag(etah(i,:)), 'b');
end
for j = 1:15:nth
    plot(real(etah(:,j)), imag(etah(:,j)), 'b');
end
plot(0, 0, 'ro')
title('\eta (\rho^2/\zeta)-plane', FS, fs)
axis image





%%
hold on
for i = 1:nr
    plot(real(zeta(i,:)), imag(zeta(i,:)), 'b');
end
for j = 1:15:nth
    plot(real(zeta(:,j)), imag(zeta(:,j)), 'b');
end
plot(0, 0, 'ro')
title('\zeta-plane', FS, fs)
axis off, axis image

