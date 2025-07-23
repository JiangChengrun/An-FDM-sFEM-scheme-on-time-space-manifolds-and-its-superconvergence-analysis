function u = laplacesolver2(mesh, f)
%% LAPLACESOLVER Fast solver for time-space laplacian
%
%   Copyright (C) Hailong Guo, Chengrun Jiang
%   16/02/2025

%% Initialize
n = size(mesh.M,1);
nt = mesh.nt;
inv = mesh.inv;
V = mesh.V;
VT = transpose(V);
IVTVT = transpose(mesh.IVTV);

%% Reshape the right hand side and time the mass matrix
F = reshape(f, n, nt);
F = mesh.M*F;
F = F*V; % transform to eigenspace
F = F*IVTVT;
%F(:,1) = F(:,1)/2;
%F(:,nt) = F(:, nt)/2;

%% Solve the laplace in space
U = zeros(n, nt);
for i = 1:nt
    U(:, i) = inv{i}\F(:,i);
end

%% Transform back in the origal space
U = U*VT;

%% Reshape for output
u = U(:);
