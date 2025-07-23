function [K, M, area] = LBMatrix(mesh)
%% LBMATRIX Stiffness and mass matrix for LB operator
%   [K, M] = LBMATRIX(MESH) returns the stiffness and mass matrix on 
%   general surface by surface finite elemene method.
%
%
%   Copyright (C) Hailong Guo, Chengrun Jiang
%   16/02/2025

%% Initialize
node = mesh.node;
elem = mesh.elem;
N = size(node,1); % number of vertices
NT = size(elem, 1); % number of triangles
NDof = N;

%% Compute geometric quantities and gradient of basis function
ve1 = node(elem(:,3),:)-node(elem(:,2),:);
ve2 = node(elem(:,1),:)-node(elem(:,3),:);
ve3 = node(elem(:,2),:)-node(elem(:,1),:);
le1 = sqrt(sum(ve1.^2,2));
le2 = sqrt(sum(ve2.^2,2));
le3 = sqrt(sum(ve3.^2,2));
cp = cross(ve1, ve2, 2);
area =  sqrt(sum(cp.^2,2))/2;
Dlambda(1:NT,:,1) = (repmat(dot(ve1,ve3,2)./le1.^2,1,3).*ve1-ve3).* ...
    repmat(le1.^2./area.^2/4,1,3);
Dlambda(1:NT,:,2) = (repmat(dot(ve2,ve1,2)./le2.^2,1,3).*ve2-ve1).* ...
    repmat(le2.^2./area.^2/4,1,3);
Dlambda(1:NT,:,3) = (repmat(dot(ve2,ve3,2)./le3.^2,1,3).*ve3-ve2).* ...
    repmat(le3.^2./area.^2/4,1,3);

%% Assemble stiffness matrix
K = sparse(NDof,NDof); % stiffness matrix
M = sparse(NDof,NDof); % mass matrix
for i = 1:3
    for j = i:3
        Kij = dot(Dlambda(:,:,i), Dlambda(:,:,j),2).*area;
        if i==j
            Mij=1/6*area;
            K = K + sparse(elem(:,i),elem(:,j),Kij,NDof,NDof);
            M =  M + sparse(elem(:,i),elem(:,j),Mij,NDof,NDof);
        else
            Mij=1/12*area;
            K = K + sparse([elem(:,i);elem(:,j)],[elem(:,j);elem(:,i)],...
                [Kij; Kij],NDof,NDof);
            M = M + sparse([elem(:,i);elem(:,j)],[elem(:,j);elem(:,i)],...
                [Mij; Mij],NDof,NDof);
        end
    end
end
clear  Kij Mij;