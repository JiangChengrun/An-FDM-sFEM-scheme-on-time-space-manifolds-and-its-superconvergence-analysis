function [l2err, h1err, h1err_t, h1terr, h1rerr, h00terr, h00rerr] = TSError(mesh, uh, pde)
%% TSERROR Error function for time-space Laplace equation
%   [L2ERR, H1ERR, H1RERR] = TSERROR(MESH, UH, PDE) compute the numerical
%   error of time-space Laplace equation using FD-FE method.
%
%   Copyright (C) Hailong Guo, Chengrun Jiang
%   16/02/2025

%% Initialize
node = mesh.node;
%elem = mesh.elem;
nxyz = size(node,1);
%NT = size(elem,1);
nt = mesh.nt;

%% Construct the recovery matrix in time and space
ht = 1/(nt-1);
e = ones(nt,1);
t = 0:ht:1;
if strcmp(mesh.rm, 'PPPR')
    %Differentiation matrix in time
    D = spdiags([-e 0*e e],-1:1,nt,nt);
    D(1,[1 2 3]) = [-3, 4, -1]; % recovery gradient approximation
    D(nt,[nt-2 nt-1 nt]) = [1, -4, 3];% recovery gradient approximation
    Bt = D/2/ht;
    % Differentiation matrix in space
    [Bx, By, Bz] = PPPRLB(mesh);
elseif strcmp(mesh.rm, 'SA')
    %Differentiation matrix in time
    Bt = 0.5*spdiags([-e 0*e e],-1:1,nt,nt)/ht;
    Bt(1, 1) = -1/ht;
    Bt(1, 2) = 1/ht;
    Bt(nt, nt-1) = -1/ht;
    Bt(nt, nt) = 1/ht;
    % Differentiation matrix in space
    [Bx, By, Bz] = SALB(mesh);
else
    error('The recovery method is not implemented')
end



%% Compute the error at each time slice
uh = reshape(uh, nxyz, nt);


rDut = uh*Bt'; % compute the differentiation in time
rDuhx = Bx*uh;% recovery the gradient
rDuhy = By*uh;% recovery the gradient
rDuhz = Bz*uh; % recovery the gradient
%rDuh = SALB2(mesh, uh);

D1=spdiags([-e,e],[-1,0],nt,nt);
D1(1,1:2) = [-1, 1];  
Bt1 = D1 / ht;
Du1t= uh*Bt1';
%backward difference

l2err = zeros(nt,1);
h1err_s = zeros(nt,1);
h1err_t = zeros(nt,1);
h1err = zeros(nt,1);
h1terr = zeros(nt,1);
h1rerr = zeros(nt,1);
h0err = zeros(nt,1);
for j = 1:nt
    exactu = @(x) pde.exactu(x, t(j));
    Du = @(x) pde.Du(x, t(j));
    DuDt = @(x) pde.DuDt(x, t(j));
    l2err(j) = L2ErrorLB(mesh,uh(:,j),exactu); %e=l2err
    h1err_s(j) = H1ErrorLB(mesh,uh(:,j),Du);     %De_2^G=h1err
    h1err_t(j) = L2ErrorLB(mesh,Du1t(:,j),DuDt); 
    h1err(j) = sqrt(h1err_s(j)^2+h1err_t(j)^2);
    %h1err(j) = H1ExactErrorLB(mesh,rDut(:,j),DuDt);
    %rDuh2 = PPPRLB2(mesh, uh(:,j));
    %h1rerr(j) = H1ExactErrorLB(mesh, rDuh2, Du);
    h1terr(j) = L2ErrorLB(mesh,rDut(:,j),DuDt);
    h1rerr(j) = H1ExactErrorLB(mesh, [rDuhx(:,j), rDuhy(:,j), rDuhz(:,j)],Du);
    % h00rerr(j) = LInfError(mesh, [rDuhx(:,j), rDuhy(:,j), rDuhz(:,j)],Du);
    h00rerr(j) = max(max((abs([rDuhx(:,j), rDuhy(:,j), rDuhz(:,j)] - Du(node)))));
    % h00rerr2(j)= max(abs([rDuhx(:,j), rDuhy(:,j), rDuhz(:,j)] - Du(node)));
    %h00rerr(j) = (abs([rDuhx(:,j), rDuhy(:,j), rDuhz(:,j)] - Du(node)));
    %h00rerr(j) = norm([rDuhx(:,j), rDuhy(:,j), rDuhz(:,j)] - Du(node),inf);
    h00terr(j) = max(abs(rDut(:,j) - DuDt(node)));
    %h0rerr(j) = max(abs([rDuhx(:,j), rDuhy(:,j), rDuhz(:,j)],Du))
end

%% Compute the maximum error
l2err = max(l2err);
h1err = max(h1err);
h1err_t = max(h1err_t);
h1terr = max(h1terr);
h1rerr = max(h1rerr);
h00rerr = max(h00rerr);
h00terr = max(h00terr);
% h1err =sqrt(sum(h1err(j).^2)/nt);
% h1rerr = sqrt(sum(h1rerr(j).^2)/nt);

