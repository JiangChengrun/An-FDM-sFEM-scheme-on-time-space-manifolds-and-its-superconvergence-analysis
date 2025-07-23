function [err, vec_err] = LInfError(mesh, rDuh, Du)
%% LINFERROR   L^{\infty} Error of finite element solution
%   ERR = LINFERROR(NODE, ELEM, UH) returns L^{\infty} error of surface finite 
%   element solution UH on the mesh given by NODE and ELEM.
%
%   Copyright (C) Hailong Guo, Chengrun Jiang
%   01/04/2022, 15/02/2025

%% Initialize
% node = mesh.node;
% elem = mesh.elem;
% Nu = size(rDuh,1);
% N = size(node,1);
% NT = size(elem,1);
% NE = N + NT;
% NP2 = N + NE;
% if Nu > N
%     elem2dof = dofP2(elem);
%     NP2 = max(elem2dof(:));
% end
% 
% %% Default quadrature order
% switch Nu
%     case N
%         quadOrder = 4;
%     case NP2
%         quadOrder = 6;
% end
% 
% %% compute gradient of finite element basis function
% ve1 = node(elem(:,3),:)-node(elem(:,2),:);
% ve2 = node(elem(:,1),:)-node(elem(:,3),:);
% ve3 = node(elem(:,2),:)-node(elem(:,1),:);
% le1 = sqrt(sum(ve1.^2,2));
% le2 = sqrt(sum(ve2.^2,2));
% le3 = sqrt(sum(ve3.^2,2));
% cp = cross(ve1, ve2, 2);
% area =  sqrt(sum(cp.^2,2))/2;
% 
% %% compute elementwise L^\infty error
% vec_err = zeros(NT, 1);
% [lambda,weight] = quadpts(quadOrder);
% nQuad = size(lambda,1);
% 
% for p = 1:nQuad
%     if Nu == N
%        rDuhp = rDuh(elem(:,1),:)*lambda(p,1) + ...
%                 rDuh(elem(:,2),:)*lambda(p,2) + ...
%                 rDuh(elem(:,3),:)*lambda(p,3);
%     elseif Nu == NP2 % piecewise quadratic function
%         % Compute the gradients of the quadratic basis functions
%         Dphip1 = (4*lambda(p,1)-1).*Dlambda(:,:,1);
%         Dphip2 = (4*lambda(p,2)-1).*Dlambda(:,:,2);
%         Dphip3 = (4*lambda(p,3)-1).*Dlambda(:,:,3);
%         Dphip4 = 4*(lambda(p,2)*Dlambda(:,:,3)+lambda(p,3)*Dlambda(:,:,2));
%         Dphip5 = 4*(lambda(p,3)*Dlambda(:,:,1)+lambda(p,1)*Dlambda(:,:,3));
%         Dphip6 = 4*(lambda(p,1)*Dlambda(:,:,2)+lambda(p,2)*Dlambda(:,:,1));
%         rDuhp = repmat(rDuh(elem2dof(:,1)),1,3).*Dphip1 + ...
%             repmat(rDuh(elem2dof(:,2)),1,3).*Dphip2 + ...
%             repmat(rDuh(elem2dof(:,3)),1,3).*Dphip3 + ...
%             repmat(rDuh(elem2dof(:,4)),1,3).*Dphip4 + ...
%             repmat(rDuh(elem2dof(:,5)),1,3).*Dphip5 + ...
%             repmat(rDuh(elem2dof(:,6)),1,3).*Dphip6;
%     end
%    pxy =lambda(p, 1)*node(elem(:, 1), :) + ...
%         lambda(p, 2)*node(elem(:, 2), :) + ...
%         lambda(p, 3)*node(elem(:, 3), :);
%    %pxy = surfacedata.project(pxy);
%    Duhxy = Du(pxy);
% 
%    % Calculate the L^\infty error for the current quadrature point
%    elem_err = max(abs(rDuhp - Duhxy), [], 2);
% 
%    % Update the elementwise error vector
%    vec_err = vec_err + weight(p) * elem_err;
% end
% 
% %% compute total L^\infty Error
% err = max(vec_err);
% end

%% H1ERROR_INF Norm of finite element solution
%   ERR = H1ERROR_INF(NODE, ELEM, UH) returns L^\infty error of surface finite 
%   element solution UH on the mesh given by NODE and ELEM.
%
%   Copyright (C) Hailong Guo
%   01/04/2022

% Initialize
% node = mesh.node; 
% elem = mesh.elem;
% Nu = size(rDuh,1);
% N = size(node,1);
% NT = size(elem,1);
% NE = N + NT;  
% NP2 = N + NE;
% if Nu > N
%     elem2dof = dofP2(elem);
%     NP2 = max(elem2dof(:));
% end
% 
% %% Default quadrature order
% switch Nu
%     case N
%         quadOrder = 4;
%     case NP2
%         quadOrder = 6;
% end
% 
% %% compute gradient of finite element basis function
% ve1 = node(elem(:,3),:)-node(elem(:,2),:);
% ve2 = node(elem(:,1),:)-node(elem(:,3),:);
% ve3 = node(elem(:,2),:)-node(elem(:,1),:);
% le1 = sqrt(sum(ve1.^2,2));
% le2 = sqrt(sum(ve2.^2,2));
% le3 = sqrt(sum(ve3.^2,2));
% cp = cross(ve1, ve2, 2);
% area =  sqrt(sum(cp.^2,2))/2;
% 
% Dlambda(1:NT,:,1) = (repmat(dot(ve1,ve3,2)./le1.^2,1,3).*ve1-ve3).* ...
%     repmat(le1.^2./area.^2/4,1,3);
% Dlambda(1:NT,:,2) = (repmat(dot(ve2,ve1,2)./le2.^2,1,3).*ve2-ve1).* ...
%     repmat(le2.^2./area.^2/4,1,3);
% Dlambda(1:NT,:,3) = (repmat(dot(ve2,ve3,2)./le3.^2,1,3).*ve3-ve2).* ...
%     repmat(le3.^2./area.^2/4,1,3);
% 
% %% compute elementwise L^\infty norm
% vec_err = zeros(NT, 1);
% [lambda,weight] = quadpts(quadOrder);
% nQuad = size(lambda,1);
% phi = lambda;
% max_err = -Inf;  % Start with negative infinity to find the maximum error
% for p = 1:nQuad
%     if Nu == N
%         rDuhp = rDuh(elem(:,1),:)*phi(p,1) + ...
%                  rDuh(elem(:,2),:)*phi(p,2) + ...
%                  rDuh(elem(:,3),:)*phi(p,3);
%     elseif Nu == NP2 % piecewise quadratic function
%         Dphip1 = (4*lambda(p,1)-1).*Dlambda(:,:,1);
%         Dphip2 = (4*lambda(p,2)-1).*Dlambda(:,:,2);
%         Dphip3 = (4*lambda(p,3)-1).*Dlambda(:,:,3);
%         Dphip4 = 4*(lambda(p,2)*Dlambda(:,:,3)+lambda(p,3)*Dlambda(:,:,2));
%         Dphip5 = 4*(lambda(p,3)*Dlambda(:,:,1)+lambda(p,1)*Dlambda(:,:,3));
%         Dphip6 = 4*(lambda(p,1)*Dlambda(:,:,2)+lambda(p,2)*Dlambda(:,:,1));
%         Duh = repmat(uh(elem2dof(:,1)),1,3).*Dphip1 + ...
%             repmat(uh(elem2dof(:,2)),1,3).*Dphip2 + ...
%             repmat(uh(elem2dof(:,3)),1,3).*Dphip3 + ...
%             repmat(uh(elem2dof(:,4)),1,3).*Dphip4 + ...
%             repmat(uh(elem2dof(:,5)),1,3).*Dphip5 + ...
%             repmat(uh(elem2dof(:,6)),1,3).*Dphip6;
%     end
%     pxy = lambda(p, 1)*node(elem(:, 1), :) + ...
%           lambda(p, 2)*node(elem(:, 2), :) + ...
%           lambda(p, 3)*node(elem(:, 3), :);
% 
%     Duhxy = Du(pxy);
%     current_err = max(abs(rDuhp - Duhxy), [], 2);  % Maximum error at this quadrature point
%     max_err = max(max_err, current_err);  % Update the global maximum error
% end
% 
% %% compute total L^\infty error
% err = max_err;
% end
% function [err, vec_err] = H1ExactErrorLB_Linfty(mesh, rDuh, Du)
%% LINFERROR   L∞ Error of finite element solution
%   ERR = LINFERROR(NODE, ELEM, UH) returns L∞ error of surface finite 
%   element solution UH on the mesh given by NODE and ELEM.
%
%   Copyright (C) Hailong Guo
%   01/04/2022

%% Initialize
% node = mesh.node; 
% elem = mesh.elem;
% Nu = size(rDuh,1);
% N = size(node,1);
% NT = size(elem,1);
% NE = N + NT;
% NP2 = N + NE;
% if Nu > N
%     elem2dof = dofP2(elem);
%     NP2 = max(elem2dof(:));
% end
% 
% %% Default quadrature order
% switch Nu
%     case N
%         quadOrder = 4;
%     case NP2
%         quadOrder = 6;
% end
% 
% %% Compute gradient of finite element basis function
% ve1 = node(elem(:,3),:) - node(elem(:,2),:);
% ve2 = node(elem(:,1),:) - node(elem(:,3),:);
% ve3 = node(elem(:,2),:) - node(elem(:,1),:);
% le1 = sqrt(sum(ve1.^2,2));
% le2 = sqrt(sum(ve2.^2,2));
% le3 = sqrt(sum(ve3.^2,2));
% cp = cross(ve1, ve2, 2);
% area = sqrt(sum(cp.^2,2)) / 2;
% 
% Dlambda(1:NT,:,1) = (repmat(dot(ve1, ve3, 2)./le1.^2, 1, 3).*ve1 - ve3) .* ...
%     repmat(le1.^2./area.^2/4, 1, 3);
% Dlambda(1:NT,:,2) = (repmat(dot(ve2, ve1, 2)./le2.^2, 1, 3).*ve2 - ve1) .* ...
%     repmat(le2.^2./area.^2/4, 1, 3);
% Dlambda(1:NT,:,3) = (repmat(dot(ve2, ve3, 2)./le3.^2, 1, 3).*ve3 - ve2) .* ...
%     repmat(le3.^2./area.^2/4, 1, 3);
% 
% %% Compute elementwise rayleigh quotient
% vec_err = zeros(NT, 1);
% [lambda, weight] = quadpts(quadOrder);
% nQuad = size(lambda, 1);
% phi = lambda;
% 
% max_err = 0; % To store the maximum error
% for p = 1:nQuad
%     if Nu == N
%         rDuhp = rDuh(elem(:,1),:) * phi(p,1) + ...
%                  rDuh(elem(:,2),:) * phi(p,2) + ...
%                  rDuh(elem(:,3),:) * phi(p,3);
%     elseif Nu == NP2 % Piecewise quadratic function
%         Dphip1 = (4*lambda(p,1) - 1) .* Dlambda(:,:,1);
%         Dphip2 = (4*lambda(p,2) - 1) .* Dlambda(:,:,2);
%         Dphip3 = (4*lambda(p,3) - 1) .* Dlambda(:,:,3);
%         Dphip4 = 4 * (lambda(p,2) * Dlambda(:,:,3) + lambda(p,3) * Dlambda(:,:,2));
%         Dphip5 = 4 * (lambda(p,3) * Dlambda(:,:,1) + lambda(p,1) * Dlambda(:,:,3));
%         Dphip6 = 4 * (lambda(p,1) * Dlambda(:,:,2) + lambda(p,2) * Dlambda(:,:,1));
%         Duh = repmat(uh(elem2dof(:,1)),1,3).*Dphip1 + ...
%               repmat(uh(elem2dof(:,2)),1,3).*Dphip2 + ...
%               repmat(uh(elem2dof(:,3)),1,3).*Dphip3 + ...
%               repmat(uh(elem2dof(:,4)),1,3).*Dphip4 + ...
%               repmat(uh(elem2dof(:,5)),1,3).*Dphip5 + ...
%               repmat(uh(elem2dof(:,6)),1,3).*Dphip6;
%     end
% 
%     % Compute the coordinates of the quadrature point
%     pxy = lambda(p, 1) * node(elem(:,1),:) + ...
%           lambda(p, 2) * node(elem(:,2),:) + ...
%           lambda(p, 3) * node(elem(:,3),:);
% 
%     % Evaluate Du at the quadrature points
%     Duhxy = Du(pxy);
% 
%     % Compute the pointwise error and update the maximum error
%     pointwise_error = max(abs(rDuhp - Duhxy), [], 2);
%     vec_err = vec_err + weight(p) * pointwise_error;
%     max_err = max(max_err, max(pointwise_error)); % Update the maximum error
% end
% 
% %% Compute total L∞ Error
% err = max_err; % The L∞ norm is the maximum error
% end
% function [err, vec_err] = H1ExactErrorLB_Linfty(mesh, rDuh, Du, uh)
%% LINFERROR   L∞ Error of finite element solution
%   ERR = LINFERROR(NODE, ELEM, UH) returns L∞ error of surface finite 
%   element solution UH on the mesh given by NODE and ELEM.
%
%   Copyright (C) Hailong Guo
%   01/04/2022

%% Initialize
% node = mesh.node; 
% elem = mesh.elem;
% Nu = size(rDuh,1);
% N = size(node,1);
% NT = size(elem,1);
% NE = N + NT;
% NP2 = N + NE;
% if Nu > N
%     elem2dof = dofP2(elem);
%     NP2 = max(elem2dof(:));
% end
% 
% %% Default quadrature order
% switch Nu
%     case N
%         quadOrder = 4;
%     case NP2
%         quadOrder = 6;
% end
% 
% %% Compute gradient of finite element basis function
% ve1 = node(elem(:,3),:) - node(elem(:,2),:);
% ve2 = node(elem(:,1),:) - node(elem(:,3),:);
% ve3 = node(elem(:,2),:) - node(elem(:,1),:);
% le1 = sqrt(sum(ve1.^2,2));
% le2 = sqrt(sum(ve2.^2,2));
% le3 = sqrt(sum(ve3.^2,2));
% cp = cross(ve1, ve2, 2);
% area = sqrt(sum(cp.^2,2)) / 2;
% 
% Dlambda(1:NT,:,1) = (repmat(dot(ve1, ve3, 2)./le1.^2, 1, 3).*ve1 - ve3) .* ...
%     repmat(le1.^2./area.^2/4, 1, 3);
% Dlambda(1:NT,:,2) = (repmat(dot(ve2, ve1, 2)./le2.^2, 1, 3).*ve2 - ve1) .* ...
%     repmat(le2.^2./area.^2/4, 1, 3);
% Dlambda(1:NT,:,3) = (repmat(dot(ve2, ve3, 2)./le3.^2, 1, 3).*ve3 - ve2) .* ...
%     repmat(le3.^2./area.^2/4, 1, 3);
% 
% %% Compute elementwise rayleigh quotient
% vec_err = zeros(NT, 1);
% [lambda, weight] = quadpts(quadOrder);
% nQuad = size(lambda, 1);
% phi = lambda;
% 
% max_err = 0; % To store the maximum error
% for p = 1:nQuad
%     if Nu == N
%         rDuhp = rDuh(elem(:,1),:) * phi(p,1) + ...
%                  rDuh(elem(:,2),:) * phi(p,2) + ...
%                  rDuh(elem(:,3),:) * phi(p,3);
%     elseif Nu == NP2 % Piecewise quadratic function
%         Dphip1 = (4*lambda(p,1) - 1) .* Dlambda(:,:,1);
%         Dphip2 = (4*lambda(p,2) - 1) .* Dlambda(:,:,2);
%         Dphip3 = (4*lambda(p,3) - 1) .* Dlambda(:,:,3);
%         Dphip4 = 4 * (lambda(p,2) * Dlambda(:,:,3) + lambda(p,3) * Dlambda(:,:,2));
%         Dphip5 = 4 * (lambda(p,3) * Dlambda(:,:,1) + lambda(p,1) * Dlambda(:,:,3));
%         Dphip6 = 4 * (lambda(p,1) * Dlambda(:,:,2) + lambda(p,2) * Dlambda(:,:,1));
%         Duh = repmat(uh(elem2dof(:,1)),1,3).*Dphip1 + ...
%               repmat(uh(elem2dof(:,2)),1,3).*Dphip2 + ...
%               repmat(uh(elem2dof(:,3)),1,3).*Dphip3 + ...
%               repmat(uh(elem2dof(:,4)),1,3).*Dphip4 + ...
%               repmat(uh(elem2dof(:,5)),1,3).*Dphip5 + ...
%               repmat(uh(elem2dof(:,6)),1,3).*Dphip6;
%     end
% 
%     % Compute the coordinates of the quadrature point
%     pxy = lambda(p, 1) * node(elem(:,1),:) + ...
%           lambda(p, 2) * node(elem(:,2),:) + ...
%           lambda(p, 3) * node(elem(:,3),:);
% 
%     % Evaluate Du at the quadrature points
%     Duhxy = Du(pxy);
% 
%     % Compute the pointwise error for each element
%     pointwise_error = max(abs(rDuhp - Duhxy), [], 2);
% 
%     % Update the maximum error globally
%     max_err = max(max_err, max(pointwise_error)); 
% 
%     % Store the elementwise error for debugging or further analysis
%     vec_err = vec_err + weight(p) * pointwise_error;
% end
% 
% %% Compute total L∞ Error
% err = max_err; % The L∞ norm is the maximum error
% end
% function [err, vec_err] = H1ExactErrorLB_Linfty(mesh, rDuh, Du, uh)
%% LINFERROR   L∞ Error of finite element solution
%   ERR = LINFERROR(NODE, ELEM, UH) returns L∞ error of surface finite 
%   element solution UH on the mesh given by NODE and ELEM.
%
%   Copyright (C) Hailong Guo
%   01/04/2022

%% Initialize
node = mesh.node; 
elem = mesh.elem;
Nu = size(rDuh,1);  % Number of basis functions or degrees of freedom
N = size(node,1);   % Number of nodes
NT = size(elem,1);  % Number of elements (triangles)
NE = N + NT;        % Extended number of nodes/elements
NP2 = N + NE;
if Nu > N
    elem2dof = dofP2(elem);  % Degrees of freedom for quadratic elements
    NP2 = max(elem2dof(:));
end

%% Default quadrature order
switch Nu
    case N
        quadOrder = 4;
    case NP2
        quadOrder = 6;
end

%% Compute gradients of finite element basis functions
ve1 = node(elem(:,3),:) - node(elem(:,2),:);
ve2 = node(elem(:,1),:) - node(elem(:,3),:);
ve3 = node(elem(:,2),:) - node(elem(:,1),:);
le1 = sqrt(sum(ve1.^2,2));
le2 = sqrt(sum(ve2.^2,2));
le3 = sqrt(sum(ve3.^2,2));
cp = cross(ve1, ve2, 2);
area = sqrt(sum(cp.^2,2)) / 2;

Dlambda(1:NT,:,1) = (repmat(dot(ve1, ve3, 2)./le1.^2, 1, 3).*ve1 - ve3) .* ...
    repmat(le1.^2./area.^2/4, 1, 3);
Dlambda(1:NT,:,2) = (repmat(dot(ve2, ve1, 2)./le2.^2, 1, 3).*ve2 - ve1) .* ...
    repmat(le2.^2./area.^2/4, 1, 3);
Dlambda(1:NT,:,3) = (repmat(dot(ve2, ve3, 2)./le3.^2, 1, 3).*ve3 - ve2) .* ...
    repmat(le3.^2./area.^2/4, 1, 3);

%% Compute elementwise Rayleigh quotient
vec_err = zeros(NT, 1);  % This can be used for debugging
[lambda, weight] = quadpts(quadOrder);  % Quadrature points and weights
nQuad = size(lambda, 1);  % Number of quadrature points
phi = lambda;

max_err = 0;  % To store the maximum error over all elements
for p = 1:nQuad
    % Evaluate rDuh at the quadrature points based on the basis functions
    if Nu == N
        rDuhp = rDuh(elem(:,1),:) * phi(p,1) + ...
                 rDuh(elem(:,2),:) * phi(p,2) + ...
                 rDuh(elem(:,3),:) * phi(p,3);
    elseif Nu == NP2  % Piecewise quadratic case
        Dphip1 = (4*lambda(p,1) - 1) .* Dlambda(:,:,1);
        Dphip2 = (4*lambda(p,2) - 1) .* Dlambda(:,:,2);
        Dphip3 = (4*lambda(p,3) - 1) .* Dlambda(:,:,3);
        Dphip4 = 4 * (lambda(p,2) * Dlambda(:,:,3) + lambda(p,3) * Dlambda(:,:,2));
        Dphip5 = 4 * (lambda(p,3) * Dlambda(:,:,1) + lambda(p,1) * Dlambda(:,:,3));
        Dphip6 = 4 * (lambda(p,1) * Dlambda(:,:,2) + lambda(p,2) * Dlambda(:,:,1));
        Duh = repmat(uh(elem2dof(:,1)),1,3).*Dphip1 + ...
              repmat(uh(elem2dof(:,2)),1,3).*Dphip2 + ...
              repmat(uh(elem2dof(:,3)),1,3).*Dphip3 + ...
              repmat(uh(elem2dof(:,4)),1,3).*Dphip4 + ...
              repmat(uh(elem2dof(:,5)),1,3).*Dphip5 + ...
              repmat(uh(elem2dof(:,6)),1,3).*Dphip6;
    end
    
    % Compute the coordinates of the quadrature point in physical space
    pxy = lambda(p, 1) * node(elem(:,1),:) + ...
          lambda(p, 2) * node(elem(:,2),:) + ...
          lambda(p, 3) * node(elem(:,3),:);

    % Evaluate Du at the quadrature points in physical space
    Duhxy = Du(pxy);
    
    % Compute the pointwise error at each quadrature point for the element
    pointwise_error = max(abs(rDuhp - Duhxy), [], 2);
    
    % Update the global maximum error across all elements
    max_err = max(max_err, max(pointwise_error));

    % Accumulate errors for debugging purposes (optional)
    vec_err = vec_err + weight(p) * pointwise_error;
end

%% Compute total L∞ Error
err = max_err;  % The L∞ norm is the maximum error over all elements
end

