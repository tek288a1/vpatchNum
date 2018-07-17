%% load data
clear;
dir = '~/Google Drive/Vpatch1/matlab1/data.matlab/';
beta_opt = 0;
n = 128;
i_rng = 1:n/2;

rho = 0.5;
% rho = input(' Input rho = ');
if beta_opt == 0
    beta = 0;
elseif beta_opt == 1
    beta = (1 - sqrt(1-rho^2))/rho;
elseif beta_opt == 2
    beta = (1 - 2*sqrt(1-rho^2))/rho;
end
file = sprintf('amat_b%1d_n%02d_r%04d.bin',...
               beta_opt, log2(n), round(rho*1e4));
filename = strcat(dir, file);
fid = fopen(filename);
A = fread(fid, [n/2, n/2], 'double');
fclose(fid);

%% inspecting (super-/sub-)diagonal terms
imagesc(A);
axis image, colorbar
hline = repmat('-', 1, 70);
fprintf('\n');
fprintf('%s \n', hline);
fprintf(' Checking validity of approximate inverse\n');
fprintf('%s \n', hline);
fprintf(' * %8s %6.4f\n', 'rho = ', rho);
fprintf(' * %8s %6.4f\n', 'beta = ', beta);
fprintf('\n');
fprintf(' * Mean of k-th diagonal terms:\n');
fprintf('\n');
fprintf('    The 2nd super-diagonal terms are dominant.\n');
fprintf('    (Check MatLab figure window for colored \n');
fprintf('     image of the matrix A.)\n');
fprintf('\n');
fprintf(' %5s %13s \n', 'k', 'mean');
fprintf(' %19s\n', repmat('-', 1, 17));
for k = 2+(-3:3)
    fprintf(' %5d %13.4e\n', k, mean(diag(A, k)));
end

%% Lsinv, approximate inverse, without block structure
nb2 = length(i_rng);
I = eye(nb2);
D = diag( -i_rng );
Dinv = inv(D);

L = D + A;
Linv = inv(L);                          % not a recommended practice
sh = 1;
d(0+sh) = mean(diag(A, 0));
d(1+sh) = mean(diag(A, 1));
d(2+sh) = mean(diag(A, 2));

N0 = d(0+sh)*Dinv*diag(ones(nb2-0,1), 0); % diagonal
N1 = d(1+sh)*Dinv*diag(ones(nb2-1,1), 1); % 1-sup. diagonal; nilpotent
N2 = d(2+sh)*Dinv*diag(ones(nb2-2,1), 2); % 2-sup. diagonal; nilpotent
Lsinv = (I - N0 - N1 - N2)*Dinv;        % Neumann series, linear
fprintf('\n');
fprintf(' * Calculating the norm || I - Ls^{I}*L ||_2: \n');
fprintf('   ( dim(L) = %d x %d )\n', size(I));
fprintf('\n')
fprintf('   ** Without block structure: %8.4f\n',...
        norm(I - Lsinv*L, 'fro'));
fprintf('      Here Ls^{I} = (I - As)*D^{-1} is a banded matrix\n');
fprintf('      of width 3 - a diagonal and 2 super-diagonals.\n');
fprintf('\n');
fprintf('      d_0 = %8.4e\n', d(0+sh));
fprintf('      d_1 = %8.4e\n', d(1+sh));
fprintf('      d_2 = %8.4e\n', d(2+sh));


%% (super-)diagonal elements 
fileOut = sprintf('diag_terms.txt');
fileOutName = strcat(dir, fileOut);
fidOut = fopen(fileOutName, 'w+');
for rho = 0.01:0.01:0.5;
    if beta_opt == 0
        beta = 0;
    elseif beta_opt == 1
        beta = (1 - sqrt(1-rho^2))/rho;
    elseif beta_opt == 2
        beta = (1 - 2*sqrt(1-rho^2))/rho;
    end
    fileIn = sprintf('amat_b%1d_n%02d_r%04d.bin',...
                     beta_opt, log2(n), round(rho*1e4));
    fileInName = strcat(dir, fileIn);
    fid = fopen(fileInName);
    A = fread(fid, [n/2, n/2], 'double');
    fclose(fid);
    
    d(0+sh) = mean(diag(A, 0));
    d(1+sh) = mean(diag(A, 1));
    d(2+sh) = mean(diag(A, 2));
    fprintf(fidOut, '%6.4f %22.16e %22.16e %22.16e\n', rho, d);
end
fclose(fidOut);

%% chebyshev interpolation
data = load(fileOutName);
rhodp = data(:, 1);
d0dp = data(:, 2);

xdp = (2*rhodp - rhodp(1) - rhodp(end))/(rhodp(end) - rhodp(1));
C = cos( acos(xdp) * (0:length(rhodp)-1) );
y = C \ d0dp;
fid = fopen('../data.matlab/cheb_coeff.txt', 'w+');
fprintf(fid, '%23.16e\n', y);
fclose(fid);
norm(C*y - d0dp, 'fro')

figure(1)
subplot(131); plot(data(:,1), data(:,2), '.-')
subplot(132); plot(data(:,1), data(:,3), '.-')
subplot(133); plot(data(:,1), data(:,4), '.-')

%% calculating chebyshev nodes
ndp = 50;
a = 0.01;
b = 0.5;
th = linspace(0, pi, ndp)';
xdp = cos(pi-th);
rhoch = (b-a)/2*xdp + (a+b)/2;
fid = fopen('../data.fortran/rho_cheb.txt', 'w+');
fprintf(fid, '%23.16e\n', rhoch);
fclose(fid);

C = cos( acos(xdp) * (0:length(xdp)-1) );
