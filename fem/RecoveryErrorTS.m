function h1rerr = RecoveryErrorTS(mesh, rDuh, pde)
%% RECOVERYERRORTS Recovered H1 Error for time space
%   [L2ERR, H1ERR, H1RERR] = TSERROR(MESH, UH, PDE) compute the numerical
%   error of time-space Laplace equation using FD-FE method.
%
%   Copyright (C) Hailong Guo, Chengrun Jiang
%   16/02/2025

%% Initialize
node = mesh.node;
%elem = mesh.elem;
N = size(node,1);
%NT = size(elem,1);
NS = mesh.n + 1;
dt = mesh.T/mesh.n;
t = 0:dt:mesh.T;


%% differentiation matrix in time



%% Compute the error at each time slice
%l2err = zeros(NS,1);
%h1err = zeros(NS,1);
h1rerr = zeros(NS,1);
h0err = zeros(NS,1);
for j = 1:NS
    index = ((j-1)*N+1):j*N;
    %exactu = @(x) pde.exactu(x, t(j));
    Du = @(x) pde.Du(x, t(j));
    DuDt = @(x) pde.DuDt(x, t(j));
    %rDuh2 = PPPRLB2(mesh, uh(:,j));
    %h1rerr(j) = H1ExactErrorLB(mesh, rDuh2, Du);
    h1rerr(j) = H1ExactErrorLB(mesh, rDuh(index,1:3),Du);
    h0err(j) = max(abs(rDuh(index,4) - DuDt(node)));
end

%% Compute the maximum error
h1rerr = max(h1rerr);
h0err = max(h0err);
h1rerr = max(h1rerr, h0err);