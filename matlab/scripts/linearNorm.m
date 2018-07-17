%% Linear norm calculation
% This script is meant to calculate/estimate the operator norm of 
%    L_s^I P_{>} (T - T_s) : h^2 --> h^2
%
% A simplified analytical matrix A_s is constructed by noting the
% dominance of near constant diagonal and some superdiagonal terms of
% A. A small number of elements on top left corner are preserved.

clear;
%% Initialize: adding path (at beginning of script)
added_path = ['../functions']; 
addpath(added_path);

%% Setting parameters
rho = 0.90;
kk = 16; % size of (1,1)-block
ld = 0;  % lower bound of diagonal index
ud = 5;  % uppder bound of diagonal index
trunc = 2;
n = 512; % size of matrix A

%% for interactive data selection
% fprintf('\n Specify rho:\n')
% rho = input(' rho = ');
% fprintf('\n For construction of a preconditioner, specify:\n')
% kk = input(' K (size of (1,1)-block) = ');
% ld = input(' lower index of diag = ');
% ud = input(' upper index of diag = ');
% trunc = input(' truncate Neumann series at = ');

%% Loading matrix data
% Matrix A calculation for large K (in this script, n) is time
% consuming and so A for different values of rho are calculated and
% stored as binary files.
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

%% Construction of preconditioner

ind_f = 1:kk;                           % indices for finite block
ind_c = kk+1:nb2;                       % complementary indices
nf = length(ind_f);
nc = length(ind_c);

I = eye(nb2);
D = diag( i_rng );
Dinv = diag( 1./i_rng );
L = D - A;

As = zeros(size(A));
Atmp = A(:, ind_c);
sh = -ld+1;
dd = zeros(ud-ld+1,1);
for k = ld:ud                        % As12 and As22: banded, diagonal
    dd(k+sh) = mean(diag(Atmp, -kk+k));
    As = As + dd(k+sh)*diag( ones(nb2-abs(k), 1), k);
end
As(ind_f, ind_f) = triu( A(ind_f, ind_f) );  % As11: upper tri.
As(ind_c, ind_f) = 0;                        % As21: zeros
N = Dinv*As;

LsI = I;
for k = 1:trunc
    LsI = LsI + N^k;
end
LsI = LsI*Dinv;             % Neumann series expansion

%% 
K = nb2/2;
Es = N + N^2;
diffA = D*(A-As)*Dinv^2;
KK = ( I + Es ) * diffA(:,1:K);
mu11 = norm( KK, 'fro' )
mu02 = 1 + norm( Es, 'fro' )


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


%%
fprintf('\n');
fprintf(' * || I - Ls^{I}*L ||_2 =  %8.4f \n', norm(I - LsI*L));
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
