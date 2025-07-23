function [iv, d, V, IVTV] = invlap2(mesh)
%% INVLAP Inverse of Laplacian matrix
%
%   Copyright (C) Hailong Guo, Chengrun Jiang
%   16/02/2025

%% Initialize
nt = mesh.nt;
ht = 1/(nt-1);
%t = linspace(0,1,nt)';
K = mesh.K;
M = mesh.M;

%% Compute the eigenvalue of time discretization matrix
m = 0:nt-1;
t = (m*ht)';
d = 2-2*cos(m*pi*ht); % eigenvalue
d = d/ht^2; % scale the eigenvalue for the finite difference matrix
V = cos(m.*pi.*t); % eigenvector
tmp = sqrt(sum(V.^2,1));
V = V./tmp;
VT = transpose(V);
VTV = VT*V;
IVTV = inv(VTV);

%% Precompute the factorization of  the matrix  using eigen decomposition
%eps = 0*1e-6;
iv = cell(nt,1);
%n = size(K, 1);
%I = speye(n);
for i = 1:nt
    %A = d(i)*M + K + eps*I + M;
    A = d(i)*M + K + M;
    dA = decomposition(A,'lu');
    iv{i} = dA;
end



