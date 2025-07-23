function [u, err] = TSLaplace2(mesh, pde)
%% TSLAPLACE Fast solver for Time-space laplace equation
%
%   Copyright (C) Hailong Guo, Chengrun Jiang
%   16/02/2025

%% Construct the matrix
[K, M] = LBMatrix(mesh);
mesh.K = K;
mesh.M = M;
[inv, d, V, IVTV] = invlap2(mesh);
mesh.inv = inv;
mesh.d = d;
mesh.V = V;
mesh.IVTV = IVTV;

%% Preprocess for TS
node = mesh.node;
nxyz = size(node,1);
nt = mesh.nt;
ht = 1/(nt-1);
x = node(:,1);
y = node(:,2);
z = node(:,3);
t = linspace(0,1,nt);
tt = repmat(t, nxyz, 1);
xx = repmat(x, nt, 1);
yy = repmat(y, nt, 1);
zz = repmat(z, nt, 1);
xx = xx(:);
yy = yy(:);
zz = zz(:);
tt = tt(:);
point = [xx, yy, zz];

%% Construct the right hand side function 
rhs = pde.f(point, tt);

%% Initial the initial conditions
idxt0 = find(abs(tt)<5*eps);
idxt1 = find(abs(tt-1)<5*eps);
bct0 = pde.g_N0(node);
bct1 = pde.g_N1(node);
%rhs(idxt0) = rhs(idxt0) - 2*bct0/ht;  % change for ghost point
%rhs(idxt1) = rhs(idxt1) + 2*bct1/ht;
rhs(idxt0) = rhs(idxt0) - 2*bct0/ht; 
rhs(idxt1) = rhs(idxt1) + 2*bct1/ht;

%% Solve the discrete system using fast solver
u = laplacesolver2(mesh, rhs);

%% Compute the error
uI = pde.exactu(point, tt);
u = u - u(1) + uI(1);
err = max(abs(u-uI));
