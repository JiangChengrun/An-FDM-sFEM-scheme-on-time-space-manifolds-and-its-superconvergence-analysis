function [err,vec_err] = H1ExactErrorLB(mesh, rDuh, Du)
%% H1ERROR   H1 Error of finite element solution
%   ERR = H1ERROR(NODE, ELEM, UH) returns H1 error of surface finite 
%   element solution UH on the mesh given by NODE and ELME.
%
%   Copyright (C) Hailong Guo, Chengrun Jiang
%   16/02/2025

%% Note there is one bug in the program to be fixed

%% Initialize
node = mesh.node;
elem = mesh.elem;
Nu = size(rDuh,1);
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

%% compute gradient of finite element basis function
ve1 = node(elem(:,3),:)-node(elem(:,2),:);
ve2 = node(elem(:,1),:)-node(elem(:,3),:);
ve3 = node(elem(:,2),:)-node(elem(:,1),:);
le1 = sqrt(sum(ve1.^2,2));
le2 = sqrt(sum(ve2.^2,2));
le3 = sqrt(sum(ve3.^2,2));
cp = cross(ve1, ve2, 2);
area =  sqrt(sum(cp.^2,2))/2;
%% Hi, Hailong, I made a mistake here. all the sign are reversed at the begining
Dlambda(1:NT,:,1) = (repmat(dot(ve1,ve3,2)./le1.^2,1,3).*ve1-ve3).* ...
    repmat(le1.^2./area.^2/4,1,3);
Dlambda(1:NT,:,2) = (repmat(dot(ve2,ve1,2)./le2.^2,1,3).*ve2-ve1).* ...
    repmat(le2.^2./area.^2/4,1,3);
Dlambda(1:NT,:,3) = (repmat(dot(ve2,ve3,2)./le3.^2,1,3).*ve3-ve2).* ...
    repmat(le3.^2./area.^2/4,1,3);

%% compute elementwise rayleigh quotient
vec_err = zeros(NT, 1);
[lambda,weight] = quadpts(quadOrder);
nQuad = size(lambda,1);
phi = lambda;
for p = 1:nQuad
    if Nu == N
       rDuhp = rDuh(elem(:,1),:)*phi(p,1) + ...
                rDuh(elem(:,2),:)*phi(p,2) + ...
                rDuh(elem(:,3),:)*phi(p,3);
    elseif Nu == NP2 % piecewise quadratic function
        Dphip1 = (4*lambda(p,1)-1).*Dlambda(:,:,1);
        Dphip2 = (4*lambda(p,2)-1).*Dlambda(:,:,2);
        Dphip3 = (4*lambda(p,3)-1).*Dlambda(:,:,3);
        Dphip4 = 4*(lambda(p,2)*Dlambda(:,:,3)+lambda(p,3)*Dlambda(:,:,2));
        Dphip5 = 4*(lambda(p,3)*Dlambda(:,:,1)+lambda(p,1)*Dlambda(:,:,3));
        Dphip6 = 4*(lambda(p,1)*Dlambda(:,:,2)+lambda(p,2)*Dlambda(:,:,1));
        Duh = repmat(uh(elem2dof(:,1)),1,3).*Dphip1 + ...
            repmat(uh(elem2dof(:,2)),1,3).*Dphip2 + ...
            repmat(uh(elem2dof(:,3)),1,3).*Dphip3 + ...
            repmat(uh(elem2dof(:,4)),1,3).*Dphip4 + ...
            repmat(uh(elem2dof(:,5)),1,3).*Dphip5 + ...
            repmat(uh(elem2dof(:,6)),1,3).*Dphip6;
    end
   pxy =lambda(p, 1)*node(elem(:, 1), :) + ...
        lambda(p, 2)*node(elem(:, 2), :) + ...
        lambda(p, 3)*node(elem(:, 3), :);
   %pxy = surfacedata.project(pxy);
   Duhxy = Du(pxy);
   vec_err = vec_err + weight(p)*sum((rDuhp - Duhxy).^2, 2);
end

%% compute total H1 Error
vec_err = vec_err.*area;
err = sqrt(sum(vec_err));