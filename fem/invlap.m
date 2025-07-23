function [inv, d, V] = invlap(mesh)
%% INVLAP Inverse of Laplacian matrix
%
%   Copyright (C) Hailong Guo, Chengrun Jiang
%   16/02/2025

%% Initalize
nt = mesh.nt;
ht = 1/(nt-1);
%t = linspace(0,1,nt)';
K = mesh.K;
M = mesh.M;

%% Compute the eigenvalue of time discretization matrix
e = ones(nt,1);
E = spdiags([-e 2*e -e],-1:1,nt,nt);
E(1,1) = 1; % we modify the diagonal part to make the matrix spd
E(nt, nt) = 1;
E = E/ht^2;
[V, D] = eig(full(E));
[~,ind] = sort(diag(D));
Ds = D(ind,ind);
V = V(:,ind);
d = abs(diag(Ds));


%% Precompute the factorization of  the matrix  using eigen decomposition
%eps = 1e-6;
inv = cell(nt,1);
%n = size(K, 1);
%I = speye(n);
for i = 1:nt
    A = d(i)*M + K + M;
    %A = d(i)*M + K + eps*I;
    dA = decomposition(A,'lu');
    inv{i} = dA;
end



