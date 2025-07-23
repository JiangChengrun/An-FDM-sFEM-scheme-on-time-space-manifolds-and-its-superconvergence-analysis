function err = L2ErrorLB(mesh, uh, exactu)
%% H1ERROR   H1 Error of finite element solution
%   ERR = H1ERROR(NODE, ELEM, UH) returns H1 error of surface finite
%   element solution UH on the mesh given by NODE and ELME.
%
%   Copyright (C) Hailong Guo, Chengrun Jiang
%   16/02/2025

%% Initialize
node = mesh.node;
elem = mesh.elem;
Nu = size(uh,1);
N = size(node,1);
NT = size(elem,1);
NE = N + NT;
NP2 = N + NE;
if Nu > N
    elem2dof = dofP2(elem);
    NP2 = max(elem2dof(:));
end

%% Default quadrature order
switch Nu
    case N
        quadOrder = 4;
    case NP2
        quadOrder = 6;
end

%% compute area of surface triangle
ve1 = node(elem(:,3),:)-node(elem(:,2),:);
ve2 = node(elem(:,1),:)-node(elem(:,3),:);
%ve3 = node(elem(:,2),:)-node(elem(:,1),:);
cp = cross(ve1, ve2, 2);
area =  sqrt(sum(cp.^2,2))/2;
 
%% basis function at quadrature points
err = zeros(NT, 1);
[lambda,weight] = quadpts(quadOrder);
switch Nu
    case N    % P1 piecewise linear function
        phi = lambda; % linear bases
    case NP2 % P2 piecewise quadratic elements
        phi(:,1) = lambda(:,1).*(2*lambda(:,1)-1);
        phi(:,2) = lambda(:,2).*(2*lambda(:,2)-1);
        phi(:,3) = lambda(:,3).*(2*lambda(:,3)-1);
        phi(:,4) = 4*lambda(:,2).*lambda(:,3);
        phi(:,5) = 4*lambda(:,1).*lambda(:,3);
        phi(:,6) = 4*lambda(:,2).*lambda(:,1);
end

%% compute elementwise rayleigh quotient
nQuad = size(lambda,1);
for p = 1:nQuad
    switch Nu
        case N % piecewise linear function
            uhp = uh(elem(:,1))*phi(p,1) + ...
                uh(elem(:,2))*phi(p,2) + ...
                uh(elem(:,3))*phi(p,3);
        case NP2 % piecewise quadratic function
            uhp = uh(elem2dof(:,1)).*phi(p,1) + ...
                uh(elem2dof(:,2)).*phi(p,2) + ...
                uh(elem2dof(:,3)).*phi(p,3) + ...
                uh(elem2dof(:,4)).*phi(p,4) + ...
                uh(elem2dof(:,5)).*phi(p,5) + ...
                uh(elem2dof(:,6)).*phi(p,6);
    end
    pxy = lambda(p, 1)*node(elem(:, 1), :) + ...
        lambda(p, 2)*node(elem(:, 2), :) + ...
        lambda(p, 3)*node(elem(:, 3), :);
    %pxy = surfacedata.project(pxy);
    err = err + weight(p)*(uhp - exactu(pxy)).^2;
end

%% compute total H1 Error
err = abs(err.*area);
err = sqrt(sum(err));