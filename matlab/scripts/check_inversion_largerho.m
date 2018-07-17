% Add path (at beginning of script)
clear;
% Add path (at beginning of script)
added_path = ['../functions']; 
addpath(added_path);

%% load data

rho = 0.90;
kk = 16;
ld = 0;  % lower bound of diagonal index
ud = 5;  % uppder bound of diagonal index
trunc = 2; % truncation
n = 512;

% fprintf('\n Specify rho:\n')
% rho = input(' rho = ');
% fprintf('\n For construction of a preconditioner, specify:\n')
% kk = input(' K (size of (1,1)-block) = ');
% ld = input(' lower index of diag = ');
% ud = input(' upper index of diag = ');
% trunc = input(' truncate Neumann series at = ');

% dir = '~/Google Drive/Vpatch1/matlab1/data.matlab/';
dir = '../../data/data.matlab/';
beta_opt = 1;
nb2 = n/2;
i_rng = 1:nb2;
sh = 1;
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

%% visually inspecting (super-/sub-)diagonal terms
figure(1);
imagesc(A);
axis image, colorbar
hline = repmat('-', 1, 70);
fprintf('\n');
fprintf('%s \n', hline);
fprintf(' Checking validity of approximate inverse:');
fprintf(' case of large rho\n');
fprintf('%s \n', hline);
fprintf(' * %8s %6.4f\n', 'rho = ', rho);
fprintf(' * %8s %6.4f\n', 'beta = ', beta);
fprintf(' * %8s %6d %s\n', 'N = ', nb2, '(size of matrix)');
fprintf(' * %8s %6d %s\n', 'K = ', kk, '(truncation)');
fprintf(' * Mean of k-th diagonal terms:\n');
fprintf('\n');
fprintf('    Super-diagonal terms are dominant.\n');
fprintf('    (Check MatLab figure window for colored \n');
fprintf('     image of the matrix A.)\n');
fprintf('\n');
fprintf(' %5s %13s \n', 'k', 'mean');
fprintf(' %19s\n', repmat('-', 1, 17));
for k = 2+(-6:6)
    fprintf(' %5d %13.4e\n', k, mean(diag(A, k)));
end

%% Construction of preconditioner

ind_f = 1:kk;
ind_c = kk+1:nb2;
nf = length(ind_f);
nc = length(ind_c);

I = eye(nb2);
D = diag( -i_rng );
Dinv = inv(D);
L = D + A;

As = zeros(size(A));
Atmp = A(:, ind_c);
sh = -ld+1;
dd = zeros(ud-ld+1,1);
for k = ld:ud                        % As12 and As22: banded, diagonal
    dd(k+sh) = mean(diag(Atmp, -kk+k));
    As = As + dd(k+sh)*diag( ones(nb2-abs(k), 1), k);
end
As(ind_f, ind_f) = triu( A(ind_f, ind_f) );% As11: upper tri.
As(ind_c, ind_f) = 0;                   % As21: zeros
N = Dinv*As;

Lsinv = I;
for k = 1:trunc
    Lsinv = Lsinv + (-N)^k;
end
Lsinv = Lsinv*Dinv;             % Neumann series expansion

figure(2)
imagesc(As), axis image, colorbar

%%
fprintf('\n');
fprintf(' * || I - Ls^{I}*L ||_2 =  %8.4f \n', norm(I - Lsinv*L));
fprintf('\n')
fprintf('   Here Ls^{I} = (I - As)*D^{-1} where:\n');
fprintf('    As_11: upper-triangular part of A;\n');
fprintf('    As_22: zero matrix;\n');
fprintf('    As_12: banded matrix with constant diagonals;\n');
fprintf('    As_22: banded matrix with constant diagonals;\n');
fprintf('\n');
fprintf(' * Diagonal terms used in construction of As\n')
fprintf('\n');
fprintf(' %5s %13s \n', 'k', 'mean');
fprintf(' %19s\n', repmat('-', 1, 17));
for k = ld:ud
    fprintf(' %5d %13.4f\n', k, dd(k+sh));
end


%% Remove path (at end of script/script clean-up)
rmpath(added_path);
